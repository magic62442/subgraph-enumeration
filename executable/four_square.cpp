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
#include <unordered_map>
#include <atomic>
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
    std::ofstream outFile;
    if (!resultPath.empty()) outFile.open(resultPath);
    std::ostream &outStream = resultPath.empty() ? std::cout : outFile;
    // Calculate total storage needed for threadTries
    ui n = dataGraph.getNumVertices();
    std::vector<ui> maxLeafSize(n, 0);
    // Allocate continuous memory for all threadTries
    std::vector<std::vector<VertexID*>> threadTries(numThreads);
    std::vector<std::vector<ui>> threadLeafSizes(numThreads);
    std::vector<std::vector<std::pair<VertexID, VertexID>>> threadTuples(numThreads);
    std::vector<ui> threadLengths(numThreads);
    ui maxNum2Hop = 0;
    size_t numEdges = dataGraph.getNumEdges();
    int threshold = 100000000;
    if (numEdges < threshold) {
        // For smaller graphs, use simple allocation
        maxNum2Hop = numEdges;
    } else {
        // For large graphs, calculate precise maxLeafSize and maxNum2Hop
#pragma omp parallel reduction(max : maxNum2Hop)
        {
#pragma omp for schedule(dynamic) nowait
            for (VertexID v = 0; v < dataGraph.getNumVertices(); ++v) {
                ui num1;
                ui num2Hop = 0;
                const VertexID *vNeighbors = dataGraph.getNeighborsLargerID(v, num1);
                for (int i = 0; i < num1; ++i) {
                    VertexID v1 = vNeighbors[i];
                    ui num2;
                    const VertexID *v1Neighbors = dataGraph.getNeighbors(v1, num2);
                    for (int j = 0; j < num2; ++j) {
                        VertexID v2 = v1Neighbors[j];
                        if (v2 <= v) continue;
                        ++num2Hop;
                    }
                }
                if (num2Hop > maxNum2Hop) maxNum2Hop = num2Hop;
            }
        }
    }

    // Start timing memory allocation
    auto memAllocStart = std::chrono::steady_clock::now();

    // Calculate storage needed per thread for threadTries
    size_t perThreadTrieStorage;
    size_t totalMaxLeafSize = 0;
    // Create separate storage for each thread
    std::vector<std::vector<VertexID>> threadTrieStorage(numThreads);

    // Parallel allocation of all memory structures
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        // Allocate threadTuples in parallel
        threadTuples[tid].resize(maxNum2Hop);
    }

    if (numEdges >= threshold) {
#pragma omp parallel
        {
#pragma omp for schedule(dynamic) nowait
            for (VertexID v = 0; v < dataGraph.getNumVertices(); ++v) {
                int tid = omp_get_thread_num();
                ui num1;
                const VertexID *vNeighbors = dataGraph.getNeighborsLargerID(v, num1);
                std::vector<std::pair<VertexID, VertexID>> &tuples = threadTuples[tid];
                std::vector<VertexID*> &trie = threadTries[tid];
                ui length = 0;
                for (int i = 0; i < num1; ++i) {
                    VertexID v1 = vNeighbors[i];
                    ui num2;
                    const VertexID *v1Neighbors = dataGraph.getNeighbors(v1, num2);
                    for (int j = 0; j < num2; ++j) {
                        VertexID v2 = v1Neighbors[j];
                        if (v2 <= v) continue;
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
    }
    if (numEdges < threshold) {
        // For smaller graphs, use degree-based allocation
        perThreadTrieStorage = numEdges;
    } else {
        // For large graphs, use maxLeafSize-based allocation
        for (VertexID v = 0; v < n; ++v) {
            totalMaxLeafSize += maxLeafSize[v];
        }
        perThreadTrieStorage = totalMaxLeafSize;
    }
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
            ui allocSize = (numEdges < threshold) ? dataGraph.getDegree(v) : maxLeafSize[v];
            threadTries[tid][v] = &threadTrieStorage[tid][offset];
            offset += allocSize;
        }
    }
    outStream << "start enumeration" << std::endl;
    // End memory allocation timing, start parallel computation timing
    auto memAllocEnd = std::chrono::steady_clock::now();
    auto parallelStart = std::chrono::steady_clock::now();

#pragma omp parallel reduction(+ : count)
    {
#pragma omp for schedule(dynamic) nowait
        for (VertexID v = 0; v < dataGraph.getNumVertices(); ++v) {
            size_t localCount = 0;
            int tid = omp_get_thread_num();
            ui num1;
            const VertexID *vNeighbors = dataGraph.getNeighborsLargerID(v, num1);
            std::vector<std::pair<VertexID, VertexID>> &tuples = threadTuples[tid];
            std::vector<VertexID*> &trie = threadTries[tid];
            ui length = 0;
            for (int i = 0; i < num1; ++i) {
                VertexID v1 = vNeighbors[i];
                ui num2;
                const VertexID *v1Neighbors = dataGraph.getNeighbors(v1, num2);
                for (int j = 0; j < num2; ++j) {
                    VertexID v2 = v1Neighbors[j];
                    if (v2 <= v) continue;
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
                for (ui k = i; k < j; ++k) {
                    trie[v2][k-i] = tuples[k].second;
                }
                ui cnt = j - i;
                localCount += cnt * (cnt - 1) / 2;
                i = j;
            }

            count += localCount;
        }
    }
    auto parallelEnd = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();

    // Calculate timing durations
    std::chrono::duration<double> memAllocTime = memAllocEnd - memAllocStart;
    std::chrono::duration<double> parallelTime = parallelEnd - parallelStart;
    std::chrono::duration<double> elapsedSeconds = end - start;

    outStream << "Number of matches: " << count << std::endl;
    outStream << "Execution Time: " << elapsedSeconds.count() << std::endl;
    outStream << "Memory Allocation Time: " << memAllocTime.count() << std::endl;
    outStream << "Parallel Computation Time: " << parallelTime.count() << std::endl;
    return 0;
}