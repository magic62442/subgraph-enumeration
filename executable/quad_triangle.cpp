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
    std::vector<VertexID *> threadCache(numThreads);
    std::vector<ui *> threadOffset(numThreads);
    for (int i = 0; i < numThreads; ++i) {
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

            if (dataGraph.getDegree(v0) >= 5) {
                ui v0_nbrs_cnt;
                const VertexID *v0_nbrs = dataGraph.getNeighbors(v0, v0_nbrs_cnt);

                // Compute common neighbors for all edges from v0
                offset[0] = 0;
                for (ui j = 0; j < v0_nbrs_cnt; ++j) {
                    VertexID v1 = v0_nbrs[j];
                    offset[j + 1] = offset[j];

                    if (dataGraph.getDegree(v1) >= 3) {
                        ui v1_nbr_cnt;
                        const VertexID *v1_nbrs = dataGraph.getNeighbors(v1, v1_nbr_cnt);
                        ui temp_count;
                        ComputeSetIntersection::ComputeCandidates(v0_nbrs, v0_nbrs_cnt,
                                                                v1_nbrs, v1_nbr_cnt,
                                                                cache + offset[j], temp_count);
                        offset[j + 1] = offset[j] + temp_count;
                    }
                }

                // Start enumeration
                for (ui j = 0; j < v0_nbrs_cnt; ++j) {
                    VertexID v1 = v0_nbrs[j];
                    ui v0_v1_cn_count = offset[j + 1] - offset[j];
                    VertexID *v0_v1_cn = cache + offset[j];

                    if (v0_v1_cn_count >= 2) {
                        for (ui k = 0; k < v0_v1_cn_count - 1; ++k) {
                            VertexID v2 = v0_v1_cn[k];
                            if (dataGraph.getDegree(v2) < 3) continue;

                            ui v2_index = ComputeSetIntersection::BinarySearch(v0_nbrs, 0, v0_nbrs_cnt, v2);
                            ui v0_v2_cn_count = offset[v2_index + 1] - offset[v2_index];
                            VertexID *v0_v2_cn = cache + offset[v2_index];

                            if (v0_v2_cn_count >= 2) {
                                for (ui l = k + 1; l < v0_v1_cn_count; ++l) {
                                    VertexID v3 = v0_v1_cn[l];
                                    if (dataGraph.getDegree(v3) < 3 || v3 <= v2) continue;

                                    ui v3_index = ComputeSetIntersection::BinarySearch(v0_nbrs, 0, v0_nbrs_cnt, v3);
                                    ui v0_v3_cn_count = offset[v3_index + 1] - offset[v3_index];
                                    VertexID *v0_v3_cn = cache + offset[v3_index];

                                    if (v0_v3_cn_count >= 2) {
                                        ui v4_temp_count = v0_v2_cn_count - 1;
                                        ui v5_temp_count = v0_v3_cn_count - 1;
                                        if (ComputeSetIntersection::Contain(v0_v2_cn, 0, v0_v2_cn_count, v3)) {
                                            v4_temp_count -= 1;
                                            v5_temp_count -= 1;
                                        }

                                        ui temp_count;
                                        ComputeSetIntersection::ComputeCandidates(v0_v2_cn, v0_v2_cn_count,
                                                                                v0_v3_cn, v0_v3_cn_count,
                                                                                temp_count);
                                        temp_count -= 1;

                                        if (v4_temp_count >= temp_count && v5_temp_count > 0) {
                                            ui temp_sum = (v4_temp_count - temp_count) * v5_temp_count +
                                                         temp_count * (v5_temp_count - 1);
                                            localCount += temp_sum;
                                        }
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
        delete[] threadCache[i];
        delete[] threadOffset[i];
    }
    return 0;
}