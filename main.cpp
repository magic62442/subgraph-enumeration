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
    patternGraph.loadPatternGraph(queryGraphPath);
    std::vector<std::vector<VertexID>> result;
    LabeledGraph q, g;
    patternGraph.loadPatternGraph(queryGraphPath);
    dataGraph.loadDataGraph(dataGraphPath);
    q.loadFromGraph(patternGraph);
    g.loadFromGraph(dataGraph);
    CandidateSpace cs(q, g, false);
    std::ofstream outFile;
    std::ostream &outStream = resultPath.empty() ? std::cout : outFile ;
    if (!resultPath.empty()) outFile.open(resultPath);
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    if (!cs.buildCandCFL()) {
        return;
    }
    cs.setQueryGraphWeights(q);
    bool *visited = new bool[g.getNumVertices()];
    memset(visited, false, sizeof(bool) * g.getNumVertices());
    VertexID *partMatch = new VertexID[q.getNumVertices()];
    VertexID **candidates = new VertexID * [q.getNumVertices()];
    VertexID **totalCandidates = new VertexID * [q.getNumVertices()];
    ui maxSize = cs.getMaxSize();
    for (int i = 0; i < q.getNumVertices(); ++i) {
        candidates[i] = new VertexID[maxSize];
        totalCandidates[i] = new VertexID[maxSize];
    }
    ui *candCount = new ui[q.getNumVertices()];
    ui *totalCandCount = new ui[q.getNumVertices()];
    memset(candCount, 0, sizeof(ui) * q.getNumVertices());
    memset(totalCandCount, 0, sizeof(ui) * q.getNumVertices());
    std::vector<ui> poses(q.getNumVertices(), 0);
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
    optCostPlan(patternGraph, g, cs, visited, partMatch, candidates, candCount, totalCandidates, totalCandCount,
                poses, tmpCand, t, pt);
    t.writeToStream(outStream);
    t.selectSymmetry(patternGraph);
    t.buildTraverseUnlabeled(q, cs, iep);
    refineIntersectionInfo(q, t, pt, intersectType);
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
    sharedJoin(t, pt, q, cs, levels, visited, result, count, beginPoses, endPoses, traverse);
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
    #ifdef PROFILE_PARALLEL
    outStream << "Par region enter: " << gParRegionEnter << std::endl;
    outStream << "Par thread enter: " << gParThreadEnter << std::endl;
    outStream << "OMP for regions: " << gOmpForRegions << std::endl;
    outStream << "OMP for iterations: " << gOmpForIters << std::endl;
    outStream << "subtreeSharedJoin calls: " << gSubtreeJoinCalls << std::endl;
    outStream << "subtreeSharedJoin time (s): " << gSubtreeJoinTime << std::endl;
    #endif
    outStream << "Number of intersections: " << gNumInterSection << std::endl;
    outStream << "Execution Time: " << elapsedSeconds.count() << std::endl;
#ifdef COLLET_GLOBAL_TIME
    outStream << "global join time: " << gTraverseTime << std::endl;
#endif

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
    LabeledGraph q, g;
    patternGraph.loadPatternGraph(queryGraphPath);
    dataGraph.loadDataGraph(dataGraphPath);
    q.loadFromGraph(patternGraph);
    g.loadFromGraph(dataGraph);
    CandidateSpace cs(q, g, false);
    std::ofstream outFile;
    std::ostream &outStream = resultPath.empty() ? std::cout : outFile ;
    if (!resultPath.empty()) outFile.open(resultPath);
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    if (!cs.buildCandCFL()) {
        return;
    }
    cs.setQueryGraphWeights(q);
    bool *visited = new bool[g.getNumVertices()];
    memset(visited, false, sizeof(bool) * g.getNumVertices());
    VertexID *partMatch = new VertexID[q.getNumVertices()];
    VertexID **candidates = new VertexID * [q.getNumVertices()];
    VertexID **totalCandidates = new VertexID * [q.getNumVertices()];
    ui maxSize = cs.getMaxSize();
    for (int i = 0; i < q.getNumVertices(); ++i) {
        candidates[i] = new VertexID[maxSize];
        totalCandidates[i] = new VertexID[maxSize];
    }
    ui *candCount = new ui[q.getNumVertices()];
    ui *totalCandCount = new ui[q.getNumVertices()];
    memset(candCount, 0, sizeof(ui) * q.getNumVertices());
    memset(totalCandCount, 0, sizeof(ui) * q.getNumVertices());
    std::vector<ui> poses(q.getNumVertices(), 0);
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
    optCostPlan(patternGraph, g, cs, visited, partMatch, candidates, candCount, totalCandidates, totalCandCount,
                poses, tmpCand, t, pt);
    t.writeToStream(outStream);
    t.selectSymmetry(patternGraph);
    t.buildTraverseUnlabeled(q, cs, iep);
    refineIntersectionInfo(q, t, pt, intersectType);
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
    if (intersectType) parSharedJoin(t, pt, q, cs, levels, visited, result, count, traverse);
    else parSharedJoin(t, pt, q, dataGraph, cs, levels, visited, result, count, traverse);
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
    #ifdef PROFILE_PARALLEL
    outStream << "Par region enter: " << gParRegionEnter << std::endl;
    outStream << "Par thread enter: " << gParThreadEnter << std::endl;
    outStream << "OMP for regions: " << gOmpForRegions << std::endl;
    outStream << "OMP for iterations: " << gOmpForIters << std::endl;
    outStream << "subtreeSharedJoin calls: " << gSubtreeJoinCalls << std::endl;
    outStream << "subtreeSharedJoin time (s): " << gSubtreeJoinTime << std::endl;
    #endif
    outStream << "Execution Time: " << elapsedSeconds.count() << std::endl;
    outStream << "Number of intersections: " << gNumInterSection << std::endl;
#ifdef COLLET_GLOBAL_TIME
    outStream << "global join time: " << gTraverseTime << std::endl;
#endif
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