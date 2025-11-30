//
// Created by Àîì÷Ñå on 25-11-3.
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
    std::vector<std::pair<VertexID, VertexID>> edges;
    edges.reserve(dataGraph.getNumEdges());
    int numThreads = cmd.getThreadNumber();
    if (numThreads == 0) numThreads = 1;
    omp_set_num_threads(numThreads);
    for (VertexID u = 0; u < dataGraph.getNumVertices(); ++u) {
        ui degree;
        const VertexID *neighbors = dataGraph.getNeighbors(u, degree);
        for (int i = 0; i < degree; ++i) {
            VertexID v = neighbors[i];
            if (v > u) edges.emplace_back(u, v);
        }
    }
    std::vector<VertexID *> candidates;
    for (int tid = 0; tid < numThreads; ++tid) {
        candidates.emplace_back(new VertexID[dataGraph.getMaxDegree()]);
    }
#pragma omp parallel reduction(+ : count)
    {
#pragma omp for schedule(dynamic) nowait
        for (int i = 0; i < edges.size(); ++i) {
            VertexID v1 = edges[i].first;
            VertexID v2 = edges[i].second;
            ui degree1, degree2;
            const VertexID *neighbors1 = dataGraph.getNeighbors(v1, degree1);
            const VertexID *neighbors2 = dataGraph.getNeighbors(v2, degree2);
            ui candCount = 0;
            int tid = omp_get_thread_num();
            ComputeSetIntersection::ComputeCandidates(neighbors1, degree1, neighbors2, degree2, candidates[tid], candCount);
            if (candCount > 1)
                count += candCount * (candCount - 1) / 2;
        }
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    std::ofstream outFile;
    if (!resultPath.empty()) outFile.open(resultPath);
    std::ostream &outStream = resultPath.empty() ? std::cout : outFile ;
    outStream << "Number of matches: " << count << std::endl;
    outStream << "Execution Time: " << elapsedSeconds.count() << std::endl;
    return 0;
}