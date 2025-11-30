#include <chrono>
#include "candidate_space.h"
#include "command.h"
#include <iomanip>
#include <omp.h>

int main(int argc, char **argv) {
    Command cmd(argc, argv);
    std::string dataGraphPath = cmd.getDataGraphPath();
    std::string resultPath = cmd.getResultPath();
    DataGraph dataGraph;
    dataGraph.loadDataGraph(dataGraphPath);

    // Calculate max_degree_sum
    ui max_degree_sum = 0;
    for (VertexID i = 0; i < dataGraph.getNumVertices(); ++i) {
        ui neighbors_count;
        const VertexID *neighbors = dataGraph.getNeighbors(i, neighbors_count);
        ui temp_sum = 0;
        for (ui j = 0; j < neighbors_count; ++j) {
            VertexID v = neighbors[j];
            temp_sum += dataGraph.getDegree(v);
        }
        if (temp_sum > max_degree_sum)
            max_degree_sum = temp_sum;
    }

    auto start = std::chrono::steady_clock::now();
    size_t count = 0;
    int numThreads = cmd.getThreadNumber();
    if (numThreads == 0) numThreads = 1;
    omp_set_num_threads(numThreads);
    std::vector<VertexID *> threadCandidates(numThreads);
    std::vector<VertexID *> threadCache(numThreads);
    std::vector<ui *> threadOffset(numThreads);
    for (int i = 0; i < numThreads; ++i) {
        threadCandidates[i] = new VertexID[dataGraph.getMaxDegree()];
        threadCache[i] = new VertexID[max_degree_sum];
        threadOffset[i] = new ui[dataGraph.getMaxDegree() + 1];
    }
#pragma omp parallel reduction(+ : count)
    {
#pragma omp for schedule(dynamic) nowait
        for (VertexID v0 = 0; v0 < dataGraph.getNumVertices(); ++v0) {
            int tid = omp_get_thread_num();
            VertexID *cache = threadCache[tid];
            ui *offset = threadOffset[tid];
            size_t localCount = 0;

            if (dataGraph.getDegree(v0) >= 4) {
                ui v0_neighbors_count;
                const VertexID *v0_neighbors = dataGraph.getNeighbors(v0, v0_neighbors_count);

                // Compute common neighbors for all edges from v0
                offset[0] = 0;
                for (ui j = 0; j < v0_neighbors_count; ++j) {
                    VertexID v1 = v0_neighbors[j];
                    offset[j + 1] = offset[j];

                    if (dataGraph.getDegree(v1) >= 3) {
                        ui v1_neighbors_count;
                        const VertexID *v1_neighbors = dataGraph.getNeighbors(v1, v1_neighbors_count);
                        ui temp_count;
                        ComputeSetIntersection::ComputeCandidates(v0_neighbors, v0_neighbors_count,
                                                                v1_neighbors, v1_neighbors_count,
                                                                cache + offset[j], temp_count);
                        offset[j + 1] = offset[j] + temp_count;
                    }
                }

                // Start enumeration
                for (ui j = 0; j < v0_neighbors_count; ++j) {
                    VertexID v1 = v0_neighbors[j];

                    if (dataGraph.getDegree(v1) >= 3) {
                        ui v0_v1_cn_count = offset[j + 1] - offset[j];
                        VertexID *v0_v1_cn = cache + offset[j];

                        if (v0_v1_cn_count >= 2) {
                            for (ui k = 0; k < v0_v1_cn_count; ++k) {
                                VertexID v2 = v0_v1_cn[k];
                                if (v2 <= v1) continue;

                                if (dataGraph.getDegree(v2) >= 3) {
                                    // Find v2's index in v0_neighbors
                                    ui v2_index = ComputeSetIntersection::BinarySearch(v0_neighbors, 0, v0_neighbors_count, v2);
                                    ui v0_v2_cn_count = offset[v2_index + 1] - offset[v2_index];
                                    VertexID *v0_v2_cn = cache + offset[v2_index];

                                    if (v0_v2_cn_count >= 2) {
                                        ui temp_count;
                                        ComputeSetIntersection::ComputeCandidates(v0_v1_cn, v0_v1_cn_count,
                                                                                v0_v2_cn, v0_v2_cn_count,
                                                                                temp_count);

                                        ui result = (v0_v1_cn_count - temp_count - 1) * (v0_v2_cn_count - 1) +
                                                   temp_count * (v0_v2_cn_count - 2);
                                        localCount += result;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            count += localCount;
        }
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    std::ofstream outFile;
    if (!resultPath.empty()) outFile.open(resultPath);
    std::ostream &outStream = resultPath.empty() ? std::cout : outFile;
    outStream << "Number of matches: " << count << std::endl;
    outStream << "Execution Time: " << elapsedSeconds.count() << std::endl;

    // Cleanup
    for (int i = 0; i < numThreads; ++i) {
        delete[] threadCandidates[i];
        delete[] threadCache[i];
        delete[] threadOffset[i];
    }
    return 0;
}