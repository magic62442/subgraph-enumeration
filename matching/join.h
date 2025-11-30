//
// Created by anonymous authors on 2024/3/8.
//

#ifndef ASDMATCH_JOIN_H
#define ASDMATCH_JOIN_H

#include "estimator.h"
#include <functional>
#include <atomic>

extern size_t gNumResult;

void twoIntersection(const VertexID *lArray, ui lBegin, ui lEnd, const VertexID *rArray, ui rBegin, ui rEnd, ui **iters, VertexID *cn, ui &cn_count);
void multiIntersection(const VertexID **values, ui num, ui **iters, ui &iterSize,
                       std::vector<ui> &beginPoses, std::vector<ui> &endPoses, VertexID *cn, ui &cn_count, ui numBag);
void
generateCandidates(const HyperTree &t, VertexID nID, int depth, const VertexID **arrays, ui *counts, ui num,
                   VertexID *cn, ui &cn_count, VertexID *partMatch, CandidateSpace &cs, VertexID current,
                   VertexID cartesianParent, std::vector<ui> &beginPoses, std::vector<ui> &endPoses);
void setIters(const VertexID *values, ui begin, ui end, ui **iters, int id, VertexID *cn, ui cn_count);
void
generateCandidates(const HyperTree &t, VertexID nID, int depth, const VertexID **arrays, ui *counts, ui num,
                   ui **iters, ui &iterSize, VertexID *cn, ui &cn_count, ui numBags,
                   VertexID *partMatch, CandidateSpace &cs, VertexID current, VertexID cartesianParent,
                   std::vector<ui> &beginPoses, std::vector<ui> &endPoses);
bool generateCandidates(const DataGraph &g, VertexID v, VertexID **allCandidates, ui *allCandCount, const std::vector<int> &verticesAfter,
                        const std::vector<std::vector<VertexID>> &verticesLarger, const std::vector<std::vector<VertexID>> &verticesSmaller, const std::vector<int> &copyAfter,
                        const std::vector<int> &copyTypes, VertexID *partMatch);
void nodeJoin(const HyperTree &t, VertexID nID, CandidateSpace &cs, bool *visited, VertexID *partMatch,
              int mappingSize, VertexID **allCandidates, ui *allCandCount, ui prefixSum, const VertexID **neighbors, ui *neighborCount, ui *allPoses,
              std::vector<std::vector<VertexID>> &tuples, TrieLevel &level, ui &length,
              std::vector<ui> &beginPoses, std::vector<ui> &endPoses);
void nodeJoin(const HyperTree &t, VertexID nID, const DataGraph &g, bool *visited, VertexID *partMatch,
              int mappingSize, VertexID **nodeCandidates, ui *nodeCandCount,
              VertexID **allCandidates, ui *allCandCount, ui prefixSum,
              const VertexID **neighbors, ui *neighborCount, ui *nodePoses,
              std::vector<std::vector<VertexID>> &tuples, TrieLevel &level, ui &length);
void sharedJoin(const HyperTree &t, const PrefixNode *pt, const Graph &query, CandidateSpace &cs,
                std::vector<TrieLevel> &trieLevels, bool *visited, std::vector<std::vector<VertexID>> &result,
                size_t &count, std::vector<ui> &beginPoses, std::vector<ui> &endPoses, bool traverse);
void globalJoin(std::vector<std::vector<VertexID>> &result, size_t &count, const Graph &query, const HyperTree &t,
                CandidateSpace &cs, VertexID *partMatch, int mappingSize, int pathLength, VertexID **allCandidates, ui *allCandCount,
                ui prefixSum, bool *visited, const std::vector<std::vector<VertexID>> &nIDs,
                const std::vector<std::vector<VertexID>> &vertexParents, const std::vector<VertexID> &cartesianParent,
                ui ***iters, ui *iterSize, bool traversal, const VertexID **neighbors, ui *neighborCount,
                std::vector<std::vector<TrieLevel *>> &traversedLevels, std::vector<TrieLevel *>& lastLevels,
                ui *allPoses, std::vector<ui> &traversePoses, std::vector<ui> &numBranches,
                TrieLevel &level, std::vector<std::vector<VertexID>> &tuples, ui &length,
                std::vector<ui> &beginPoses, std::vector<ui> &endPoses);
void globalJoin(std::vector<std::vector<VertexID>> &result, size_t &count, const Graph &query, const HyperTree &t,
                const DataGraph &g, CandidateSpace &cs, VertexID *partMatch, int mappingSize, int pathLength,
                VertexID **nodeCandidates, ui *nodeCandCount, VertexID **allCandidates, ui *allCandCount,
                ui prefixSum, bool *visited, const std::vector<std::vector<VertexID>> &nIDs,
                const std::vector<std::vector<VertexID>> &vertexParents, const std::vector<VertexID> &cartesianParent,
                ui ***iters, ui *iterSize, bool traversal, const VertexID **neighbors, ui *neighborCount,
                std::vector<std::vector<TrieLevel *>> &traversedLevels, std::vector<TrieLevel *>& lastLevels,
                ui *nodePoses, std::vector<ui> &traversePoses, std::vector<ui> &numBranches,
                TrieLevel &level, std::vector<std::vector<VertexID>> &tuples, ui &length,
                std::vector<ui> &beginPoses, std::vector<ui> &endPoses);
void enumerate(std::vector<std::vector<VertexID>> &result, const HyperTree &t, const std::vector<VertexID> &order,
               VertexID *partMatch, int mappingSize, std::vector<std::vector<TrieLevel *>> &traversedLevels, bool *visited);
void
traverse(size_t &count, const Graph &query, const HyperTree &t, const std::vector<VertexID> &order, VertexID *partMatch,
         int mappingSize, const std::vector<TrieLevel *> &levels, std::vector<std::vector<TrieLevel *>> &traversedNodes,
         std::vector<ui> &poses, std::vector<ui> &numBranches, bool *visited, int extendLevel, bool labeled);
void
produceResult(std::vector<std::vector<VertexID>> &result, VertexID *partMatch, ui size);
void storeMatches(const std::vector<std::vector<VertexID>> &result, std::ofstream &outStream);
void buildTrie(std::vector<std::vector<VertexID>> &tuples, TrieLevel &root, ui &length);

#endif //ASDMATCH_JOIN_H
