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

void singleThreadMatch(int argc, char **argv) {
    Command cmd(argc, argv);
    std::string queryGraphPath = cmd.getQueryGraphPath();
    std::string dataGraphPath = cmd.getDataGraphPath();
    std::string resultPath = cmd.getResultPath();
    bool iep = cmd.getUseIEP();
    bool intersectType = !cmd.getIntersectType();
    PatternGraph patternGraph;
    DataGraph dataGraph;
    std::vector<std::vector<VertexID>> result;
    patternGraph.loadPatternGraph(queryGraphPath);
    dataGraph.loadDataGraph(dataGraphPath);
    CandidateSpace cs(patternGraph, dataGraph);
    std::ofstream outFile;
    std::ostream &outStream = resultPath.empty() ? std::cout : outFile ;
    if (!resultPath.empty()) outFile.open(resultPath);
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    if (!cs.buildCandCFL()) {
        return;
    }
    cs.setQueryGraphWeights(patternGraph);
    bool *visited = new bool[dataGraph.getNumVertices()];
    memset(visited, false, sizeof(bool) * dataGraph.getNumVertices());
    VertexID *partMatch = new VertexID[patternGraph.getNumVertices()];
    VertexID **candidates = new VertexID * [patternGraph.getNumVertices()];
    VertexID **totalCandidates = new VertexID * [patternGraph.getNumVertices()];
    ui maxSize = cs.getMaxSize();
    for (int i = 0; i < patternGraph.getNumVertices(); ++i) {
        candidates[i] = new VertexID[maxSize];
        totalCandidates[i] = new VertexID[maxSize];
    }
    ui *candCount = new ui[patternGraph.getNumVertices()];
    ui *totalCandCount = new ui[patternGraph.getNumVertices()];
    memset(candCount, 0, sizeof(ui) * patternGraph.getNumVertices());
    memset(totalCandCount, 0, sizeof(ui) * patternGraph.getNumVertices());
    std::vector<ui> poses(patternGraph.getNumVertices(), 0);
    std::vector<VertexID> tmpCand(maxSize);
    size_t count = 0;
    bool traverse = false;
#ifdef ALL_LEVEL
    traverse = true;
#endif
    start = std::chrono::steady_clock::now();
    HyperTree t;
    PrefixNode *pt;
    double intersectCost = 0.0, materializeCost = 0.0;
    optCostPlan(patternGraph, dataGraph, cs, visited, partMatch, candidates, candCount, totalCandidates, totalCandCount,
                poses, tmpCand, t, pt);
    t.writeToStream(outStream);
    t.selectSymmetry(patternGraph);
    t.buildTraverseUnlabeled(patternGraph, cs, iep);
    refineIntersectionInfo(patternGraph, t, pt, intersectType);
    std::vector<TrieLevel> levels(t.numNodes);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    outStream << "Planning time (s): " << elapsedSeconds.count() << std::endl;
    outStream << "traverse all level: " << traverse << std::endl;
    outStream << "use iep: " << iep << std::endl;
    if (iep && !(t.symmLastLevel.empty() && t.subsetLastLevel.empty())) {
        outStream << "optimized traversal" << std::endl;
    }
    for (int i = 0; i < levels.size(); ++i) {
        if (t.trieOrder[i].size() == 1) levels[i].oneLevel = true;
        else levels[i].oneLevel = false;
    }
    std::vector<ui> beginPoses(t.numAttributes, 0);
    std::vector<ui> endPoses(t.numAttributes, 0);
    start = std::chrono::steady_clock::now();
    sharedJoin(t, pt, patternGraph, cs, levels, visited, result, count, beginPoses, endPoses, traverse);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    count /= t.divideFactor;
#ifndef ALL_LEVEL
    gNumResult = count;
#endif
#ifdef COLLECT_RESULT
    gNumResult = result.size();
#endif
    if (gNumResult == 0) gNumResult = count;
    outStream << "Number of matches: " << gNumResult << std::endl;
    outStream << "Number of intersections: " << gNumInterSection << std::endl;
    outStream << "Execution Time: " << elapsedSeconds.count() << std::endl;

}

void parallelMatch(int argc, char **argv) {
    Command cmd(argc, argv);
    std::string queryGraphPath = cmd.getQueryGraphPath();
    std::string dataGraphPath = cmd.getDataGraphPath();
    std::string resultPath = cmd.getResultPath();
    bool iep = cmd.getUseIEP();
    bool intersectType = !cmd.getIntersectType();
    int numThreads = cmd.getThreadNumber();
    if (numThreads == 0) numThreads = 1;
    if (numThreads > 0) {
        omp_set_num_threads(numThreads);
    }

    PatternGraph patternGraph;
    DataGraph dataGraph;
    std::vector<std::vector<VertexID>> result;
    patternGraph.loadPatternGraph(queryGraphPath);
    dataGraph.loadDataGraph(dataGraphPath);
    CandidateSpace cs(patternGraph, dataGraph);
    std::ofstream outFile;
    std::ostream &outStream = resultPath.empty() ? std::cout : outFile ;
    if (!resultPath.empty()) outFile.open(resultPath);
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    if (!cs.buildCandCFL()) {
        return;
    }
    cs.setQueryGraphWeights(patternGraph);
    bool *visited = new bool[dataGraph.getNumVertices()];
    memset(visited, false, sizeof(bool) * dataGraph.getNumVertices());
    VertexID *partMatch = new VertexID[patternGraph.getNumVertices()];
    VertexID **candidates = new VertexID * [patternGraph.getNumVertices()];
    VertexID **totalCandidates = new VertexID * [patternGraph.getNumVertices()];
    ui maxSize = cs.getMaxSize();
    for (int i = 0; i < patternGraph.getNumVertices(); ++i) {
        candidates[i] = new VertexID[maxSize];
        totalCandidates[i] = new VertexID[maxSize];
    }
    ui *candCount = new ui[patternGraph.getNumVertices()];
    ui *totalCandCount = new ui[patternGraph.getNumVertices()];
    memset(candCount, 0, sizeof(ui) * patternGraph.getNumVertices());
    memset(totalCandCount, 0, sizeof(ui) * patternGraph.getNumVertices());
    std::vector<ui> poses(patternGraph.getNumVertices(), 0);
    std::vector<VertexID> tmpCand(maxSize);
    size_t count = 0;
    bool traverse = false;
#ifdef ALL_LEVEL
    traverse = true;
#endif
    start = std::chrono::steady_clock::now();
    HyperTree t;
    PrefixNode *pt;
    double intersectCost = 0.0, materializeCost = 0.0;
    optCostPlan(patternGraph, dataGraph, cs, visited, partMatch, candidates, candCount, totalCandidates, totalCandCount,
                poses, tmpCand, t, pt);
    t.writeToStream(outStream);
    t.selectSymmetry(patternGraph);
    t.buildTraverseUnlabeled(patternGraph, cs, iep);
    refineIntersectionInfo(patternGraph, t, pt, intersectType);
    std::vector<TrieLevel> levels(t.numNodes);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    outStream << "Planning time (s): " << elapsedSeconds.count() << std::endl;
    outStream << "traverse all level: " << traverse << std::endl;
    outStream << "use iep: " << iep << std::endl;
    outStream << "Number of threads: " << (numThreads > 0 ? numThreads : omp_get_max_threads()) << std::endl;
    if (iep && !(t.symmLastLevel.empty() && t.subsetLastLevel.empty())) {
        outStream << "optimized traversal" << std::endl;
    }
    for (int i = 0; i < levels.size(); ++i) {
        if (t.trieOrder[i].size() == 1) levels[i].oneLevel = true;
        else levels[i].oneLevel = false;
    }
    start = std::chrono::steady_clock::now();
    if (intersectType) parSharedJoin(t, pt, patternGraph, cs, levels, visited, result, count, traverse);
    else parSharedJoin(t, pt, patternGraph, dataGraph, cs, levels, visited, result, count, traverse);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    count /= t.divideFactor;
#ifndef ALL_LEVEL
    gNumResult = count;
#endif
#ifdef COLLECT_RESULT
    gNumResult = result.size();
#endif
    if (gNumResult == 0) gNumResult = count;
    outStream << "Number of matches: " << gNumResult << std::endl;
    outStream << "Execution Time: " << elapsedSeconds.count() << std::endl;
    outStream << "Number of intersections: " << gNumInterSection << std::endl;
}

int main(int argc, char **argv) {
    Command cmd(argc, argv);
    bool useParallel = false;
    int numThreads = cmd.getThreadNumber();
    if (numThreads > 1) useParallel = true;
    if (useParallel) {
        parallelMatch(argc, argv);
    } else {
        singleThreadMatch(argc, argv);
    }
    return 0;
}
