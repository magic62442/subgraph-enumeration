//
// Created by Qiyan LI on 25-11-3.
//

#include <chrono>
#include "candidate_space.h"
#include "join.h"
#include "command.h"
#include "optimizer.h"
#include "estimator.h"
#include "all_td.h"
#include "parallel.h"
#include <iomanip>
#include <omp.h>
#include "parallel.h"

int main(int argc, char **argv) {
    Command cmd(argc, argv);
    std::string dataGraphPath = cmd.getDataGraphPath();
    std::string resultPath = cmd.getResultPath();
    DataGraph dataGraph;
    dataGraph.loadDataGraph(dataGraphPath);
    auto start = std::chrono::steady_clock::now();
    size_t count = 0;
    int numThreads = cmd.getThreadNumber();
    if (numThreads == 0) numThreads = 1;
    omp_set_num_threads(numThreads);
    size_t numEdges = dataGraph.getNumEdges();
    ui n = dataGraph.getNumVertices();
    int threshold1 = 100000000, threshold2 = threshold1 * 10;
    omp_set_num_threads(numThreads);
    std::ofstream outFile;
    if (!resultPath.empty()) outFile.open(resultPath);
    std::ostream &outStream = resultPath.empty() ? std::cout : outFile;
    ui maxNum2Hop = 0;
    size_t perThreadTrieStorage;
    size_t totalMaxLeafSize = 0;
    std::vector<std::vector<VertexID>> threadTrieStorage(numThreads);
    std::vector<std::vector<VertexID*>> threadTries(numThreads);
    std::vector<std::vector<ui>> threadTrieSizes(numThreads);
    std::vector<std::vector<ui>> threadLeafSizes(numThreads);
    std::vector<std::vector<std::pair<VertexID, VertexID>>> threadTuples(numThreads);
    std::vector<ui> threadLengths(numThreads);
    std::vector<VertexID *> threadCandidates(numThreads);
    std::vector<ui> maxLeafSize(n, 0);
    for (int i = 0; i < numThreads; ++i) {
        threadCandidates[i] = new VertexID [dataGraph.getMaxDegree()];
        threadTrieSizes[i] = std::vector<ui>(n, 0);
    }
    if (numEdges < threshold1) {
        maxNum2Hop = numEdges;
    }
    else {
#pragma omp parallel reduction(max : maxNum2Hop)
        {
#pragma omp for schedule(dynamic) nowait
            for (VertexID v = 0; v < dataGraph.getNumVertices(); ++v) {
                ui num1;
                ui num2Hop = 0;
                const VertexID *vNeighbors = dataGraph.getNeighbors(v, num1);
                for (int i = 0; i < num1; ++i) {
                    VertexID v1 = vNeighbors[i];
                    ui num2;
                    const VertexID *v1Neighbors = dataGraph.getNeighbors(v1, num2);
                    for (int j = 0; j < num2; ++j) {
                        VertexID v2 = v1Neighbors[j];
                        if (v2 == v) continue;
                        ++num2Hop;
                    }
                }
                if (num2Hop > maxNum2Hop) maxNum2Hop = num2Hop;
            }
        }
    }
    // Start timing memory allocation
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        // Allocate threadTuples in parallel
        threadTuples[tid].resize(maxNum2Hop);
    }
    if (numEdges < threshold1) {
        totalMaxLeafSize = numEdges;
    }
    else if (numEdges < threshold2) {
#pragma omp parallel
        {
#pragma omp for schedule(dynamic) nowait
            for (VertexID v = 0; v < dataGraph.getNumVertices(); ++v) {
                int tid = omp_get_thread_num();
                ui num1;
                const VertexID *vNeighbors = dataGraph.getNeighbors(v, num1);
                std::vector<std::pair<VertexID, VertexID>> &tuples = threadTuples[tid];
                ui length = 0;
                for (int i = 0; i < num1; ++i) {
                    VertexID v1 = vNeighbors[i];
                    ui num2;
                    const VertexID *v1Neighbors = dataGraph.getNeighbors(v1, num2);
                    for (int j = 0; j < num2; ++j) {
                        VertexID v2 = v1Neighbors[j];
                        if (v2 == v) continue;
                        tuples[length].first = v2;
                        tuples[length].second = v1;
                        length++;
                    }
                }
                std::sort(tuples.begin(), tuples.begin() + length);
                ui i = 0, j = 0;
                while (i < length) {
                    VertexID v2 = tuples[i].first;
                    j = i + 1;
                    while (j < length && tuples[j].first == v2) {
                        j++;
                    }
                    ui cnt = j - i;
                    if (cnt > maxLeafSize[v2])
                        maxLeafSize[v2] = cnt;
                    i = j;
                }
            }
        }
        for (VertexID v = 0; v < n; ++v) {
            totalMaxLeafSize += maxLeafSize[v];
        }
    }
    perThreadTrieStorage = totalMaxLeafSize;
    // Calculate total memory to be allocated
    size_t trieDataSize = perThreadTrieStorage * numThreads * sizeof(VertexID);  // trie data storage
    size_t tuplesSize = maxNum2Hop * numThreads * sizeof(std::pair<VertexID, VertexID>);  // tuples storage
    size_t triePointersSize = n * numThreads * sizeof(VertexID*);  // trie pointers storage
    size_t totalMemoryBytes = trieDataSize + tuplesSize + triePointersSize;
    double totalMemoryGB = totalMemoryBytes / (1024.0 * 1024.0 * 1024.0);
    outStream << "Estimated memory allocation: " << std::fixed << std::setprecision(2) << totalMemoryGB << " GB" << std::endl;
    outStream << "  Tuples: " << (tuplesSize / (1024.0 * 1024.0 * 1024.0)) << " GB" << std::endl;
    outStream << "total leaf size: " << totalMaxLeafSize << ", total number of edges: " << numEdges << std::endl;
    outStream << "  Trie data: " << (trieDataSize / (1024.0 * 1024.0 * 1024.0)) << " GB" << std::endl;
    if (numEdges < threshold2) {
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            // Allocate independent storage for this thread
            threadTrieStorage[tid].resize(perThreadTrieStorage);

            // Resize threadTries for this thread
            threadTries[tid].resize(n);

            // Assign pointers to this thread's memory
            size_t offset = 0;
            for (VertexID v = 0; v < n; ++v) {
                ui allocSize = (numEdges < threshold1) ? dataGraph.getDegree(v) : maxLeafSize[v];
                threadTries[tid][v] = &threadTrieStorage[tid][offset];
                offset += allocSize;
            }
        }
    }
    outStream << "start enumeration" << std::endl;
    if (numEdges <= threshold2) {
#pragma omp parallel reduction(+ : count)
        {
#pragma omp for schedule(dynamic) nowait
            for (VertexID v = 0; v < dataGraph.getNumVertices(); ++v) {
                size_t localCount = 0;
                int tid = omp_get_thread_num();
                VertexID *candidates2 = threadCandidates[tid], *candidates3;
                ui candCount2 = 0, candCount3 = 0;
                ui num1;
                const VertexID *vNeighbors = dataGraph.getNeighbors(v, num1);
                std::vector<std::pair<VertexID, VertexID>> &tuples = threadTuples[tid];
                std::vector<VertexID*> &trie = threadTries[tid];
                std::vector<ui> &trieSize = threadTrieSizes[tid];
                ui length = 0;
                for (int i = 0; i < num1; ++i) {
                    VertexID v3 = vNeighbors[i];
                    ui num2;
                    const VertexID *v1Neighbors = dataGraph.getNeighbors(v3, num2);
                    for (int j = 0; j < num2; ++j) {
                        VertexID v4 = v1Neighbors[j];
                        if (v4 == v) continue;
                        tuples[length].first = v4;
                        tuples[length].second = v3;
                        length++;
                    }
                }
                std::sort(tuples.begin(), tuples.begin() + length);
                ui pos1 = 0, pos2 = 0;
                while (pos1 < length) {
                    VertexID v4 = tuples[pos1].first;
                    pos2 = pos1 + 1;
                    while (pos2 < length && tuples[pos2].first == v4) {
                        pos2++;
                    }
                    for (ui k = pos1; k < pos2; ++k) {
                        trie[v4][k - pos1] = tuples[k].second;
                    }
                    ui cnt = pos2 - pos1;
                    trieSize[v4] = cnt;
                    pos1 = pos2;
                }
                for (int i = 0; i < num1; ++i) {
                    VertexID v1 = vNeighbors[i];
                    if (v1 <= v) continue;
                    ui num2;
                    const VertexID *v1Neighbors = dataGraph.getNeighbors(v1, num2);
                    ComputeSetIntersection::ComputeCandidates(vNeighbors, num1, v1Neighbors, num2, candidates2, candCount2);
                    if (candCount2 < 1) continue;
                    for (int j = 0; j < num2; ++j) {
                        VertexID v4 = v1Neighbors[j];
                        if (v4 == v) continue;
                        candidates3 = trie[v4];
                        candCount3 = trieSize[v4];
                        if (candCount3 < 2) continue;
                        ui tmpCount = 0;
                        ComputeSetIntersection::ComputeCandidates(candidates2, candCount2, candidates3, candCount3, tmpCount);
                        ui tmp = candCount2;
                        if (ComputeSetIntersection::Contain(candidates2, 0, candCount2, v4)) --tmp;
                        ui numResult = (candCount3 - 2) * tmpCount + (candCount3 - 1) * (tmp - tmpCount);
                        localCount += numResult;
                    }
                }
                count += localCount;
            }
        }
    }
    else {
#pragma omp parallel reduction(+ : count)
        {
#pragma omp for schedule(dynamic) nowait
            for (VertexID v = 0; v < dataGraph.getNumVertices(); ++v) {
                size_t localCount = 0;
                int tid = omp_get_thread_num();
                VertexID *candidates2 = threadCandidates[tid], *candidates3;
                ui candCount2 = 0, candCount3 = 0;
                ui num1;
                const VertexID *vNeighbors = dataGraph.getNeighbors(v, num1);
                std::vector<std::pair<VertexID, VertexID>> &tuples = threadTuples[tid];
                std::map<VertexID, std::vector<VertexID>> trie;
                std::map<VertexID, ui> trieSize;
                ui length = 0;
                for (int i = 0; i < num1; ++i) {
                    VertexID v3 = vNeighbors[i];
                    ui num2;
                    const VertexID *v3Neighbors = dataGraph.getNeighbors(v3, num2);
                    for (int j = 0; j < num2; ++j) {
                        VertexID v4 = v3Neighbors[j];
                        if (v4 == v) continue;
                        tuples[length].first = v4;
                        tuples[length].second = v3;
                        length++;
                    }
                }
                std::sort(tuples.begin(), tuples.begin() + length);
                ui pos1 = 0, pos2 = 0;
                while (pos1 < length) {
                    VertexID v4 = tuples[pos1].first;
                    pos2 = pos1 + 1;
                    while (pos2 < length && tuples[pos2].first == v4) {
                        pos2++;
                    }
                    ui cnt = pos2 - pos1;
                    trie[v4] = std::vector<VertexID>(cnt);
                    for (ui k = pos1; k < pos2; ++k) {
                        trie[v4][k - pos1] = tuples[k].second;
                    }
                    trieSize[v4] = cnt;
                    pos1 = pos2;
                }
                for (int i = 0; i < num1; ++i) {
                    VertexID v1 = vNeighbors[i];
                    if (v1 <= v) continue;
                    ui num2;
                    const VertexID *v1Neighbors = dataGraph.getNeighbors(v1, num2);
                    ComputeSetIntersection::ComputeCandidates(vNeighbors, num1, v1Neighbors, num2, candidates2, candCount2);
                    if (candCount2 < 1) continue;
                    for (int j = 0; j < num2; ++j) {
                        VertexID v4 = v1Neighbors[j];
                        if (v4 == v) continue;
                        candidates3 = trie[v4].data();
                        candCount3 = trieSize[v4];
                        if (candCount3 < 2) continue;
                        ui tmpCount = 0;
                        ComputeSetIntersection::ComputeCandidates(candidates2, candCount2, candidates3, candCount3, tmpCount);
                        ui tmp = candCount2;
                        if (ComputeSetIntersection::Contain(candidates2, 0, candCount2, v4)) --tmp;
                        ui numResult = (candCount3 - 2) * tmpCount + (candCount3 - 1) * (tmp - tmpCount);
                        localCount += numResult;
                    }
                }
                count += localCount;
            }
        }
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    outStream << "Number of matches: " << count << std::endl;
    outStream << "Execution Time: " << elapsedSeconds.count() << std::endl;
    return 0;
}