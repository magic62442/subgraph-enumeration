//
// Created by lqy on 2025/1/28.
//

#ifndef GRAPH_MINING_SYSTEM_PARALLEL_H
#define GRAPH_MINING_SYSTEM_PARALLEL_H

#include "join.h"
#include <omp.h>
#include <atomic>
#include <vector>
#include <map>

#ifdef PROFILE_PARALLEL
// 诊断统计
extern long long gParRegionEnter;        // 并行区进入次数（single）
extern long long gParThreadEnter;        // 并行区线程进入次数（所有线程总和）
extern long long gOmpForRegions;         // omp for 区域次数
extern long long gOmpForIters;           // omp for 总迭代数
extern long long gSubtreeJoinCalls;      // 并行区内 subtreeSharedJoin 调用次数
extern double gSubtreeJoinTime;          // 并行区内 subtreeSharedJoin 累计时间（秒）
#endif


// subtreeSharedJoin 的上下文，用于在调用前分配、调用后统一释放内部使用的临时内存
struct alignas(64) SubtreeSharedJoinContext {
    // 按节点展平后的候选与计数
    VertexID **nodeCandidates = nullptr;
    ui *nodeCandCount = nullptr;
    ui *nodePoses = nullptr;

    // 所有中间结果的候选与计数（每个位置(nID,i)分配一段连续空间用于存储中间结果）
    VertexID **allCandidates = nullptr;
    ui *allCandCount = nullptr;

    // 迭代器缓冲
    ui ***iters = nullptr;
    ui *iterSizes = nullptr;

    // 路径候选（前缀树遍历）
    VertexID **pCandidates = nullptr;
    ui *pCandCount = nullptr;

    // 邻接缓存
    const VertexID **neighbors = nullptr;
    ui *neighborCount = nullptr;

    // 运行期辅助向量/状态
    std::vector<ui> prefixSum;           // 大小 t.numNodes + 1
    std::vector<ui> traversePoses;       // 大小 t.numAttributes
    std::vector<ui> numBranches;         // 大小 t.numAttributes
    std::vector<ui> beginPoses;         // 大小 t.numAttributes
    std::vector<ui> endPoses;           // 大小 t.numAttributes
    std::vector<ui> pPoses;              // 大小 height
    std::vector<ui> childPoses;          // 大小 height
    std::vector<const PrefixNode *> nodes; // 大小 height
    ui height = 0; // subtreeRoot->getHeight() + mappingSize
};

// 上下文初始化与释放接口
void initSubtreeSharedJoinContext(SubtreeSharedJoinContext &ctx, const HyperTree &t, const CandidateSpace &cs,
                                  const PrefixNode *subtreeRoot, int mappingSize);
void freeSubtreeSharedJoinContext(SubtreeSharedJoinContext &ctx, const HyperTree &t, const CandidateSpace &cs, const PrefixNode *subtreeRoot);

// 子树版本的 sharedJoin 函数（使用外部分配的上下文内存）
// 用于处理一个rooted subtree，假设从root开始的路径已经匹配了mappingSize个属性
// 这个版本不是多线程的，等待外部函数并行调用
void subtreeSharedJoin(const HyperTree &t, const PrefixNode *subtreeRoot, int mappingSize,
                       const std::vector<int> &mappingSizes, const std::map<const PrefixNode *, std::vector<VertexID>> &bagsBelow,
                       const Graph &query, CandidateSpace &cs,
                       std::vector<TrieLevel> &trieLevels, std::vector<std::vector<TrieLevel *>> &traversedLevels,
                       std::vector<TrieLevel *> &lastLevels, bool *visited, VertexID *partMatch,
                       std::vector<std::vector<VertexID>> &result, size_t &count, bool traverse,
                       std::vector<std::vector<std::vector<VertexID>>> &localTuples, std::vector<ui> &localLengths,
                       SubtreeSharedJoinContext &ctx);
void subtreeSharedJoin(const HyperTree &t, const PrefixNode *subtreeRoot, int mappingSize,
                       const std::vector<int> &mappingSizes, const std::map<const PrefixNode *, std::vector<VertexID>> &bagsBelow,
                       const Graph &query, const DataGraph &g, CandidateSpace &cs,
                       std::vector<TrieLevel> &trieLevels, std::vector<std::vector<TrieLevel *>> &traversedLevels,
                       std::vector<TrieLevel *> &lastLevels, bool *visited, VertexID *partMatch,
                       std::vector<std::vector<VertexID>> &result, size_t &count, bool traverse,
                       std::vector<std::vector<std::vector<VertexID>>> &localTuples, std::vector<ui> &localLengths,
                       SubtreeSharedJoinContext &ctx);

void parNodeJoin(const HyperTree &t, VertexID nID, CandidateSpace &cs, VertexID **allCandidates, ui *allCandCount,
                 ui prefixSum, std::vector<std::vector<VertexID>> &tuples, ui &length);

void parSharedJoin(const HyperTree &t, const PrefixNode *pt, const Graph &query, CandidateSpace &cs,
    std::vector<TrieLevel> &trieLevels, bool *visited, std::vector<std::vector<VertexID>> &result,
    size_t &count, bool traverse);

void parSharedJoin(const HyperTree &t, const PrefixNode *pt, const Graph &query, const DataGraph &g, CandidateSpace &cs,
                   std::vector<TrieLevel> &trieLevels, bool *visited, std::vector<std::vector<VertexID>> &result,
                   size_t &count, bool traverse);

#endif //GRAPH_MINING_SYSTEM_PARALLEL_H