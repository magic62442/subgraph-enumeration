//
// Created by Qiyan LI on 25-11-3.
// Five Cycle: u_0 -- u_1 -- u_2 -- u_3 -- u_4 -- u_0
// Rules: v_0 < v_1, v_0 < v_2, v_0 < v_3, v_0 < v_4, v_3 < v_2
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
#include "compute_set_intersection.h"

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

    auto memAllocStart = std::chrono::steady_clock::now();

    size_t perThreadTrieStorage;
    size_t totalMaxLeafSize = 0;
    std::vector<std::vector<VertexID>> threadTrieStorage(numThreads);
    // Storage for leaf sizes per thread
    std::vector<std::vector<ui>> threadTrieLeafSizes(numThreads);

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
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
        perThreadTrieStorage = numEdges;
    } else {
        for (VertexID v = 0; v < n; ++v) {
            totalMaxLeafSize += maxLeafSize[v];
        }
        perThreadTrieStorage = totalMaxLeafSize;
    }

    size_t trieDataSize = perThreadTrieStorage * numThreads * sizeof(VertexID);
    size_t tuplesSize = maxNum2Hop * numThreads * sizeof(std::pair<VertexID, VertexID>);
    size_t triePointersSize = n * numThreads * sizeof(VertexID*);
    size_t leafSizesSize = n * numThreads * sizeof(ui);
    size_t totalMemoryBytes = trieDataSize + tuplesSize + triePointersSize + leafSizesSize;
    double totalMemoryGB = totalMemoryBytes / (1024.0 * 1024.0 * 1024.0);
    outStream << "Estimated memory allocation: " << std::fixed << std::setprecision(2) << totalMemoryGB << " GB" << std::endl;
    outStream << "  Tuples: " << (tuplesSize / (1024.0 * 1024.0 * 1024.0)) << " GB" << std::endl;
    outStream << "total leaf size: " << totalMaxLeafSize << ", total number of edges: " << numEdges << std::endl;
    outStream << "  Trie data: " << (trieDataSize / (1024.0 * 1024.0 * 1024.0)) << " GB" << std::endl;

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        threadTrieStorage[tid].resize(perThreadTrieStorage);
        threadTries[tid].resize(n);
        threadTrieLeafSizes[tid].resize(n, 0);

        size_t offset = 0;
        for (VertexID v = 0; v < n; ++v) {
            ui allocSize = (numEdges < threshold) ? dataGraph.getDegree(v) : maxLeafSize[v];
            threadTries[tid][v] = &threadTrieStorage[tid][offset];
            offset += allocSize;
        }
    }
    outStream << "start enumeration" << std::endl;
    auto memAllocEnd = std::chrono::steady_clock::now();
    auto parallelStart = std::chrono::steady_clock::now();

#pragma omp parallel reduction(+ : count)
    {
#pragma omp for schedule(dynamic) nowait
        for (VertexID v0 = 0; v0 < dataGraph.getNumVertices(); ++v0) {
            size_t localCount = 0;
            int tid = omp_get_thread_num();
            ui num1;
            const VertexID *v0Neighbors = dataGraph.getNeighborsLargerID(v0, num1);
            std::vector<std::pair<VertexID, VertexID>> &tuples = threadTuples[tid];
            std::vector<VertexID*> &trie = threadTries[tid];
            std::vector<ui> &leafSizes = threadTrieLeafSizes[tid];

            // Step 1: Build trie for wedge u_0 -- u_1 -- u_2
            // Condition: v_0 < v_1, v_0 < v_2
            // trie[v2] stores all v1 values (sorted)
            ui length = 0;
            for (ui i = 0; i < num1; ++i) {
                VertexID v1 = v0Neighbors[i];  // v1 > v0 (from getNeighborsLargerID)
                ui num2;
                const VertexID *v1Neighbors = dataGraph.getNeighbors(v1, num2);
                for (ui j = 0; j < num2; ++j) {
                    VertexID v2 = v1Neighbors[j];
                    if (v2 <= v0) continue;  // v2 > v0
                    tuples[length].first = v2;
                    tuples[length].second = v1;
                    length++;
                }
            }
            std::sort(tuples.begin(), tuples.begin() + length);

            // Build trie: for each v2, store sorted list of v1 values
            ui i = 0, j = 0;
            while (i < length) {
                VertexID v2 = tuples[i].first;
                j = i + 1;
                while (j < length && tuples[j].first == v2) {
                    j++;
                }
                // Store v1 values for this v2
                for (ui k = i; k < j; ++k) {
                    trie[v2][k - i] = tuples[k].second;
                }
                leafSizes[v2] = j - i;
                i = j;
            }

            // Step 2: Enumerate u_0 -- u_4 -- u_3 -- u_2 path
            // Conditions: v_0 < v_4, v_0 < v_3, v_0 < v_2, v_3 < v_2
            for (ui i4 = 0; i4 < num1; ++i4) {
                VertexID v4 = v0Neighbors[i4];  // v4 > v0
                ui numV4Neighbors;
                const VertexID *v4Neighbors = dataGraph.getNeighborsLargerID(v0, numV4Neighbors);
                v4Neighbors = dataGraph.getNeighbors(v4, numV4Neighbors);

                for (ui i3 = 0; i3 < numV4Neighbors; ++i3) {
                    VertexID v3 = v4Neighbors[i3];
                    if (v3 <= v0) continue;  // v3 > v0

                    ui numV3Neighbors;
                    const VertexID *v3Neighbors = dataGraph.getNeighborsLargerID(v3, numV3Neighbors);

                    for (ui i2 = 0; i2 < numV3Neighbors; ++i2) {
                        VertexID v2 = v3Neighbors[i2];  // v2 > v3 > v0
                        if (v2 == v4) continue;
                        // Check trie[v2] for counting
                        ui leafSize = leafSizes[v2];
                        if (leafSize == 0) continue;

                        // Count valid v1 values: all v1 in trie[v2] except v3 and v4
                        ui validCount = leafSize;

                        // Check if v3 is in trie[v2] (sorted array)
                        VertexID *leaf = trie[v2];
                        ui pos3 = ComputeSetIntersection::BinarySearch(leaf, 0, leafSize, v3);
                        if (pos3 < leafSize && leaf[pos3] == v3) {
                            validCount--;
                        }
                        // Check if v4 is in trie[v2] (only if v4 != v3)
                        ui pos4 = ComputeSetIntersection::BinarySearch(leaf, 0, leafSize, v4);
                        if (pos4 < leafSize && leaf[pos4] == v4) {
                            validCount--;
                        }
                        localCount += validCount;
                    }
                }
            }

            // Reset leaf sizes for next iteration
            i = 0;
            while (i < length) {
                VertexID v2 = tuples[i].first;
                j = i + 1;
                while (j < length && tuples[j].first == v2) {
                    j++;
                }
                leafSizes[v2] = 0;
                i = j;
            }

            count += localCount;
        }
    }
    auto parallelEnd = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double> memAllocTime = memAllocEnd - memAllocStart;
    std::chrono::duration<double> parallelTime = parallelEnd - parallelStart;
    std::chrono::duration<double> elapsedSeconds = end - start;

    outStream << "Number of matches: " << count << std::endl;
    outStream << "Execution Time: " << elapsedSeconds.count() << std::endl;
    outStream << "Memory Allocation Time: " << memAllocTime.count() << std::endl;
    outStream << "Parallel Computation Time: " << parallelTime.count() << std::endl;
    return 0;
}
