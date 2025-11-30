//
// Created by lqy on 2025/1/28.
//

#include "parallel.h"
#include <chrono>

void initSubtreeSharedJoinContext(SubtreeSharedJoinContext &ctx, const HyperTree &t, const CandidateSpace &cs,
                                  const PrefixNode *subtreeRoot, int mappingSize) {
    ctx.prefixSum.assign(t.numNodes + 1, 0);
    for (VertexID nID = 1; nID <= t.numNodes; ++nID) {
        ctx.prefixSum[nID] = ctx.prefixSum[nID - 1] + t.nodes[nID - 1].numAttributes;
    }
    std::vector<VertexID> below = subtreeRoot->getBagsBelow();
    ui totalFlat = ctx.prefixSum.back();
    ui maxSize = cs.getMaxSize();

    ctx.nodeCandidates = new VertexID *[totalFlat];
    ctx.nodeCandCount = new ui[totalFlat];
    ctx.nodePoses = new ui[totalFlat];
    for (ui i = 0; i < totalFlat; ++i) {
        ctx.nodeCandidates[i] = nullptr;
        ctx.nodeCandCount[i] = 0;
        ctx.nodePoses[i] = 0;
    }
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        ui prefixSum = ctx.prefixSum[nID];
        for (int i = 1; i < t.nodes[nID].numAttributes; ++i) {
            ctx.nodeCandidates[prefixSum + i] = new VertexID [maxSize];
        }
    }

    // 计算每个位置(nID, i)需要的段长度，并分配allCandidates和allCandCount
    ui totalAllCandSize = 0;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            ui segmentSize = (t.nodes[nID].candidatesBefore[i] != 0) + t.nodes[nID].attributesBefore[i].size();
            totalAllCandSize += segmentSize;
        }
    }

    ctx.allCandidates = new VertexID *[totalAllCandSize];
    ctx.allCandCount = new ui[totalAllCandSize];
    for (ui i = 0; i < totalAllCandSize; ++i) {
        ctx.allCandidates[i] = nullptr;
        ctx.allCandCount[i] = 0;
    }

    // 为每段分配内存并设置最后一个位置指向nodeCandidates
    ui allCandIndex = 0;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        ui nodePrefixSum = ctx.prefixSum[nID];
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            ui segmentSize = (t.nodes[nID].candidatesBefore[i] != 0) + t.nodes[nID].attributesBefore[i].size();
            if (segmentSize == 0) continue;
            for (ui j = 0; j < segmentSize - 1; ++j) {
                ctx.allCandidates[allCandIndex + j] = new VertexID[maxSize];
            }
            // 最后一个位置指向nodeCandidates（如果存在的话）
            if (i > 0) { // nodeCandidates[0]通常为nullptr，跳过
                ctx.allCandidates[allCandIndex + segmentSize - 1] = ctx.nodeCandidates[nodePrefixSum + i];
            }
            allCandIndex += segmentSize;
        }
    }
    if (subtreeRoot->pathToGlobal) {
        ui size = t.nodes[t.numNodes - 1].numAttributes;
        ctx.iters = new ui **[size];
        ctx.iterSizes = new ui[size];
        for (int i = 0; i < size; ++i) { ctx.iters[i] = nullptr; ctx.iterSizes[i] = 0; }
        for (int i = mappingSize; i < size; ++i) {
            const auto &nIDs_i = t.nodes[t.numNodes - 1].nIDs[i];
            if (!nIDs_i.empty()) {
                ctx.iters[i] = new ui *[maxSize];
                ui cols = nIDs_i.size();
                ui *block = new ui[maxSize * cols];
                for (ui j = 0; j < maxSize; ++j) {
                    ctx.iters[i][j] = block + j * cols;
                }
            }
        }
    }
    ctx.height = subtreeRoot->getHeight() + mappingSize;
    ctx.pCandidates = new VertexID *[ctx.height];
    ctx.pCandCount = new ui[ctx.height];
    for (ui i = 0; i < ctx.height; ++i) {
        ctx.pCandidates[i] = nullptr; // 非拥有指针，运行期赋值为 nodeCandidates[..] 或者不使用
        ctx.pCandCount[i] = 0;
    }
    ctx.neighbors = new const VertexID *[t.numAttributes];
    ctx.neighborCount = new ui[t.numAttributes];

    ctx.traversePoses.assign(t.numAttributes, 0);
    ctx.numBranches.assign(t.numAttributes, 0);
    ctx.beginPoses.assign(t.numAttributes, 0);
    ctx.endPoses.assign(t.numAttributes, 0);
    ctx.pPoses.assign(ctx.height, 0);
    ctx.childPoses.assign(ctx.height, 0);
    ctx.nodes.assign(ctx.height, nullptr);
}

void freeSubtreeSharedJoinContext(SubtreeSharedJoinContext &ctx, const HyperTree &t, const CandidateSpace &cs, const PrefixNode *subtreeRoot) {
    if (ctx.nodeCandidates) {
        for (ui i = 0; i < ctx.prefixSum.back(); ++i) {
            if (ctx.nodeCandidates[i]) {
                delete[] ctx.nodeCandidates[i];
                ctx.nodeCandidates[i] = nullptr;
            }
        }
        delete[] ctx.nodeCandidates;
        ctx.nodeCandidates = nullptr;
    }
    delete[] ctx.nodeCandCount; ctx.nodeCandCount = nullptr;
    delete[] ctx.nodePoses; ctx.nodePoses = nullptr;

    // 释放allCandidates中的中间结果内存
    if (ctx.allCandidates) {
        ui allCandIndex = 0;
        for (VertexID nID = 0; nID < t.numNodes; ++nID) {
            for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
                ui segmentSize = (t.nodes[nID].candidatesBefore[i] != 0) + t.nodes[nID].attributesBefore[i].size();
                if (segmentSize == 0) continue;
                // 只释放中间结果的内存，不释放最后一个位置（它指向nodeCandidates）
                for (ui j = 0; j < segmentSize - 1; ++j) {
                    if (ctx.allCandidates[allCandIndex + j]) {
                        delete[] ctx.allCandidates[allCandIndex + j];
                        ctx.allCandidates[allCandIndex + j] = nullptr;
                    }
                }
                allCandIndex += segmentSize;
            }
        }
        delete[] ctx.allCandidates; ctx.allCandidates = nullptr;
    }
    if (ctx.allCandCount) {
        delete[] ctx.allCandCount; ctx.allCandCount = nullptr;
    }

    if (ctx.iters) {
        int size = t.nodes[t.numNodes - 1].numAttributes;
        for (int i = 0; i < size; ++i) {
            if (ctx.iters[i]) {
                // 仅释放一次：整块连续内存的基地址在 iters[i][0]
                delete[] ctx.iters[i][0];
                delete[] ctx.iters[i];
            }
        }
        delete[] ctx.iters;
        ctx.iters = nullptr;
    }
    delete[] ctx.iterSizes; ctx.iterSizes = nullptr;
    delete[] ctx.iterSizes; ctx.iterSizes = nullptr;

    delete[] ctx.pCandidates; ctx.pCandidates = nullptr;
    delete[] ctx.pCandCount; ctx.pCandCount = nullptr;

    delete[] ctx.neighbors; ctx.neighbors = nullptr;
    delete[] ctx.neighborCount; ctx.neighborCount = nullptr;
    ctx.prefixSum.clear();
    ctx.traversePoses.clear();
    ctx.numBranches.clear();
    ctx.pPoses.clear();
    ctx.childPoses.clear();
    ctx.nodes.clear();
}

void subtreeSharedJoin(const HyperTree &t, const PrefixNode *subtreeRoot, int mappingSize,
                       const std::vector<int> &mappingSizes, const std::map<const PrefixNode *, std::vector<VertexID>> &bagsBelow,
                       const Graph &query, CandidateSpace &cs,
                       std::vector<TrieLevel> &trieLevels, std::vector<std::vector<TrieLevel *>> &traversedLevels,
                       std::vector<TrieLevel *> &lastLevels, bool *visited, VertexID *partMatch,
                       std::vector<std::vector<VertexID>> &result, size_t &count, bool traverse,
                       std::vector<std::vector<std::vector<VertexID>>> &localTuples, std::vector<ui> &localLengths,
                       SubtreeSharedJoinContext &ctx) {
    std::vector<ui> &prefixSum = ctx.prefixSum;
    VertexID **nodeCandidates = ctx.nodeCandidates;
    ui *nodeCandCount = ctx.nodeCandCount;
    ui *nodePoses = ctx.nodePoses;
    const HyperNode &globalNode = t.nodes[t.numNodes - 1];
    const PrefixNode *pn = subtreeRoot;
    ui maxSize = cs.getMaxSize();
    ui ***iters = ctx.iters;
    ui *iterSizes = ctx.iterSizes;
    std::vector<ui> &traversePoses = ctx.traversePoses;
    std::vector<ui> &numBranches = ctx.numBranches;
    ui height = ctx.height;
    VertexID **pCandidates = ctx.pCandidates;
    ui *pCandCount = ctx.pCandCount;
    std::fill(pCandCount, pCandCount + height, 0);
    std::vector<ui> &pPoses = ctx.pPoses;
    std::vector<ui> &childPoses = ctx.childPoses;
    std::fill(pPoses.begin(), pPoses.end(), 0);
    std::fill(childPoses.begin(), childPoses.end(), 0);
    const VertexID **neighbors = ctx.neighbors;
    ui *neighborCount = ctx.neighborCount;
    ui maxEdgeSize = maxNumBackWard(t.globalOrder, query);
    ui maxNum = height + globalNode.numAttributes - globalNode.prefixSize;
    std::vector<ui> &beginPoses = ctx.beginPoses;
    std::vector<ui> &endPoses = ctx.endPoses;

    pn = subtreeRoot;
    std::vector<const PrefixNode *> &nodes = ctx.nodes;
    nodes.assign(height, nullptr);
    if (mappingSize > 0) nodes[mappingSize - 1] = subtreeRoot;
    int depth = mappingSize;

    // 处理子树根节点的nIDsToCall
    for (VertexID nID: pn -> nIDsToCall) {
        if (nID != t.numNodes - 1) {
            nodeJoin(t, nID, cs, visited, partMatch, mappingSizes[nID], nodeCandidates, nodeCandCount,
                     prefixSum[nID], neighbors, neighborCount, nodePoses, localTuples[nID], trieLevels[nID], localLengths[nID],
                     beginPoses, endPoses);
        }
    }

    // 如果子树根有子节点，为第一个子节点创建候选集
    // 由于已经有了匹配过的点，需要利用已匹配的点来计算候选集
    if (!pn -> children.empty()) {
        const PrefixNode *current = pn -> children[0];
        VertexID u = current -> u;
        if (current -> pathToGlobal) {
            bool flag = true;
            for (VertexID nID : pn -> nIDsToBuild) {
                if (localLengths[nID] == 0) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                for (VertexID nID : pn -> nIDsToBuild)
                    if (!trieLevels[nID].oneLevel) buildTrie(localTuples[nID], trieLevels[nID], localLengths[nID]);
            }
            else {
                return;
            }
            VertexID firstBag = bagsBelow.at(current)[0];
            if (current->nIDsToJoin.empty()) {
                const std::vector<VertexID> &parents = t.nodes[firstBag].attributesBefore[depth];
                int reuse = 0;
                if (t.nodes[firstBag].candidatesBefore[depth] != 0) {
                    reuse = 1;
                    neighbors[0] = nodeCandidates[t.nodes[firstBag].candidatesBefore[depth]];
                    neighborCount[0] = nodeCandCount[t.nodes[firstBag].candidatesBefore[depth]];
                }
                for (int i = 0; i < parents.size(); ++i) {
                    VertexID pU = parents[i];
                    neighbors[i + reuse] = cs.getNeighbors(pU, partMatch[pU], u, neighborCount[i + reuse]);
                }
                generateCandidates(t, firstBag, depth, neighbors, neighborCount, parents.size() + reuse,
                                   nodeCandidates[prefixSum[firstBag] + depth], nodeCandCount[prefixSum[firstBag] + depth],
                                   partMatch, cs, u, current->cartesianParent, beginPoses, endPoses);
                pCandidates[depth] = nodeCandidates[prefixSum[firstBag] + depth];
                pCandCount[depth] = nodeCandCount[prefixSum[firstBag] + depth];
                iterSizes[depth] = pCandCount[depth];
            }
            else {
                for (int i = 0; i < current->nIDsToJoin.size(); ++i) {
                    VertexID nID = current->nIDsToJoin[i];
                    neighbors[i] = lastLevels[nID]->values;
                    neighborCount[i] = lastLevels[nID]->length;
                }
                for (int i = 0; i < current->attributesBefore.size(); ++i) {
                    VertexID pU = current->attributesBefore[i];
                    VertexID pV = partMatch[pU];
                    ui numNeighbors;
                    const VertexID *pVNeighbors = cs.getNeighbors(pU, pV, u, numNeighbors);
                    neighbors[current->nIDsToJoin.size() + i] = pVNeighbors;
                    neighborCount[current->nIDsToJoin.size() + i] = numNeighbors;
                }
                generateCandidates(t, t.numNodes - 1, depth, neighbors, neighborCount, current->nIDsToJoin.size() + current->attributesBefore.size(),
                                   iters[depth], iterSizes[depth], nodeCandidates[prefixSum[firstBag] + depth], nodeCandCount[prefixSum[firstBag] + depth], current->nIDsToJoin.size(),
                                   partMatch, cs, u, current->cartesianParent, beginPoses, endPoses);
            }
        }
        else {
            if (depth != 0) {
                VertexID firstBag = bagsBelow.at(current)[0];
                const std::vector<VertexID> &parents = t.nodes[firstBag].attributesBefore[depth];
                int reuse = 0;
                if (t.nodes[firstBag].candidatesBefore[depth] != 0) {
                    reuse = 1;
                    neighbors[0] = nodeCandidates[t.nodes[firstBag].candidatesBefore[depth]];
                    neighborCount[0] = nodeCandCount[t.nodes[firstBag].candidatesBefore[depth]];
                }
                for (int i = 0; i < parents.size(); ++i) {
                    VertexID pU = parents[i];
                    neighbors[i + reuse] = cs.getNeighbors(pU, partMatch[pU], u, neighborCount[i + reuse]);
                }
                generateCandidates(t, firstBag, depth, neighbors, neighborCount, parents.size() + reuse,
                                   nodeCandidates[prefixSum[firstBag] + depth], nodeCandCount[prefixSum[firstBag] + depth],
                                   partMatch, cs, u, current->cartesianParent, beginPoses, endPoses);
                pCandidates[depth] = nodeCandidates[prefixSum[firstBag] + depth];
                pCandCount[depth] = nodeCandCount[prefixSum[firstBag] + depth];
            }
            else {
                pCandCount[depth] = cs.candSize(u);
                memcpy(pCandidates[depth], cs.candSet(u), pCandCount[depth] * sizeof(VertexID));
            }
        }
    }

    // 主循环：遍历子树（简化版本，主要处理子树结构）
    while (depth >= mappingSize) {
        while (childPoses[depth] < pn -> children.size()) {
            const PrefixNode *current = pn -> children[childPoses[depth]];
            nodes[depth] = current;
            VertexID u = current -> u;
            bool nextLevel = false;

            if (!current -> pathToGlobal) {
                while (pPoses[depth] < pCandCount[depth]) {
                    VertexID v = pCandidates[depth][pPoses[depth]];
                    ++pPoses[depth];
                    if (visited[v]) continue;
                    visited[v] = true;
                    partMatch[u] = v;

                    for (VertexID nID: current -> nIDsToCall) {
                        nodeJoin(t, nID, cs, visited, partMatch, mappingSizes[nID], nodeCandidates,
                                 nodeCandCount, prefixSum[nID], neighbors, neighborCount, nodePoses, localTuples[nID],
                                 trieLevels[nID], localLengths[nID], beginPoses, endPoses);
                    }

                    if (!current -> children.empty()) {
                        nextLevel = true;
                        break;
                    }
                    else visited[v] = false;
                }
            }
            else {
                while (pPoses[depth] < iterSizes[depth]) {
                    VertexID v;
                    if (current->nIDsToJoin.empty()) {
                        v = pCandidates[depth][pPoses[depth]];
                    }
                    else {
                        ui childPos = iters[depth][pPoses[depth]][0];
                        VertexID tmp = current -> nIDsToJoin[0];
                        v = traversedLevels[tmp].back()->values[childPos];

                        if (!visited[v]) {
                            for (int i = 0; i < current -> nIDsToJoin.size(); ++i) {
                                VertexID nID = current -> nIDsToJoin[i];
                                childPos = iters[depth][pPoses[depth]][i];
                                if (traversedLevels[nID].back()->children) {
                                    traversedLevels[nID].push_back(&traversedLevels[nID].back()->children[childPos]);
                                    lastLevels[nID] = traversedLevels[nID].back();
                                }
                                else {
                                    traversedLevels[nID].push_back(nullptr);
                                    lastLevels[nID] = nullptr;
                                }
                            }
                        }
                    }
                    ++pPoses[depth];
                    if (visited[v]) continue;
                    visited[v] = true;
                    partMatch[u] = v;

                    for (VertexID nID: current -> nIDsToCall) {
                        if (nID != t.numNodes - 1) {
                            nodeJoin(t, nID, cs, visited, partMatch, mappingSizes[nID], nodeCandidates,
                                     nodeCandCount, prefixSum[nID], neighbors, neighborCount,
                                     nodePoses, localTuples[nID], trieLevels[nID], localLengths[nID],
                                     beginPoses, endPoses);
                        }
                    }

                    if (!current -> children.empty()) {
                        nextLevel = true;
                        break;
                    }
                    else {
                        for (VertexID nID: current -> nIDsToCall) {
                            if (nID == t.numNodes - 1) {
                                bool flag = true;
                                for (VertexID nID2 : current -> nIDsToBuild) {
                                    if (localLengths[nID2] == 0) {
                                        flag = false;
                                        break;
                                    }
                                }
                                if (flag) {
                                    for (VertexID nID2 : current -> nIDsToBuild)
                                        if (!trieLevels[nID2].oneLevel) buildTrie(localTuples[nID2], trieLevels[nID2], localLengths[nID2]);
                                    globalJoin(result, count, query, t, cs, partMatch, mappingSizes[nID],
                                               depth + 1, nodeCandidates, nodeCandCount, prefixSum[nID],
                                               visited, globalNode.nIDs, globalNode.attributesBefore,
                                               globalNode.cartesianParent, iters, iterSizes, traverse, neighbors, neighborCount,
                                               traversedLevels, lastLevels, nodePoses, traversePoses, numBranches,
                                               trieLevels[nID], localTuples[nID], localLengths[nID],
                                               beginPoses, endPoses);
                                }
                                else {
                                    for (VertexID nID2 : current -> nIDsToBuild) {
                                        localLengths[nID2] = 0;
                                    }
                                }
                            }
                        }
                        visited[v] = false;
                        for (VertexID nID : current -> nIDsToJoin) {
                            traversedLevels[nID].pop_back();
                            lastLevels[nID] = traversedLevels[nID].back();
                        }
                    }
                }
            }

            if (!nextLevel) {
                ++childPoses[depth];
                pPoses[depth] = 0;
                if (childPoses[depth] == pn->children.size()) break;
                current = pn->children[childPoses[depth]];
                u = current->u;
            }
            else {
                ++depth;
                pPoses[depth] = 0;
                childPoses[depth] = 0;
                pn = current;
                current = pn -> children[0];
                u = current->u;
            }

            if (current -> pathToGlobal) {
                bool flag = true;
                for (VertexID nID : pn -> nIDsToBuild) {
                    if (localLengths[nID] == 0) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    for (VertexID nID : pn -> nIDsToBuild)
                        if (!trieLevels[nID].oneLevel) buildTrie(localTuples[nID], trieLevels[nID], localLengths[nID]);
                }
                else {
                    for (VertexID nID2 : pn -> nIDsToBuild) {
                        localLengths[nID2] = 0;
                    }
                    visited[partMatch[pn->u]] = false;
                    for (VertexID nID : pn -> nIDsToJoin) {
                        traversedLevels[nID].pop_back();
                        lastLevels[nID] = traversedLevels[nID].back();
                    }
                    --depth;
                    if (depth < mappingSize) {
                        pn = subtreeRoot;
                        break;
                    }
                    else pn = nodes[depth - 1];
                    continue;
                }
                VertexID firstBag = bagsBelow.at(current)[0];
                if (current->nIDsToJoin.empty()) {
                    const std::vector<VertexID> &parents = t.nodes[firstBag].attributesBefore[depth];
                    int reuse = 0;
                    if (t.nodes[firstBag].candidatesBefore[depth] != 0) {
                        reuse = 1;
                        neighbors[0] = nodeCandidates[t.nodes[firstBag].candidatesBefore[depth]];
                        neighborCount[0] = nodeCandCount[t.nodes[firstBag].candidatesBefore[depth]];
                    }
                    for (int i = 0; i < parents.size(); ++i) {
                        VertexID pU = parents[i];
                        neighbors[i + reuse] = cs.getNeighbors(pU, partMatch[pU], u, neighborCount[i + reuse]);
                    }
                    generateCandidates(t, firstBag, depth, neighbors, neighborCount, parents.size() + reuse,
                                       nodeCandidates[prefixSum[firstBag] + depth], nodeCandCount[prefixSum[firstBag] + depth],
                                       partMatch, cs, u, current->cartesianParent,
                                       beginPoses, endPoses);
                    pCandidates[depth] = nodeCandidates[prefixSum[firstBag] + depth];
                    pCandCount[depth] = nodeCandCount[prefixSum[firstBag] + depth];
                    iterSizes[depth] = pCandCount[depth];
                }
                else {
                    for (int i = 0; i < current->nIDsToJoin.size(); ++i) {
                        VertexID nID = current->nIDsToJoin[i];
                        neighbors[i] = lastLevels[nID]->values;
                        neighborCount[i] = lastLevels[nID]->length;
                    }
                    for (int i = 0; i < current->attributesBefore.size(); ++i) {
                        VertexID pU = current->attributesBefore[i];
                        VertexID pV = partMatch[pU];
                        ui numNeighbors;
                        const VertexID *pVNeighbors = cs.getNeighbors(pU, pV, u, numNeighbors);
                        neighbors[current->nIDsToJoin.size() + i] = pVNeighbors;
                        neighborCount[current->nIDsToJoin.size() + i] = numNeighbors;
                    }
                    generateCandidates(t, t.numNodes - 1, depth, neighbors, neighborCount, current->nIDsToJoin.size() + current->attributesBefore.size(),
                                       iters[depth], iterSizes[depth], nodeCandidates[prefixSum[firstBag] + depth], nodeCandCount[prefixSum[firstBag] + depth], current->nIDsToJoin.size(),
                                       partMatch, cs, u, current->cartesianParent,
                                       beginPoses, endPoses);
                }
            }
            else {
                if (depth != 0) {
                    VertexID firstBag = bagsBelow.at(current)[0];
                    const std::vector<VertexID> &parents = t.nodes[firstBag].attributesBefore[depth];
                    int reuse = 0;
                    if (t.nodes[firstBag].candidatesBefore[depth] != 0) {
                        reuse = 1;
                        neighbors[0] = nodeCandidates[t.nodes[firstBag].candidatesBefore[depth]];
                        neighborCount[0] = nodeCandCount[t.nodes[firstBag].candidatesBefore[depth]];
                    }
                    for (int i = 0; i < parents.size(); ++i) {
                        VertexID pU = parents[i];
                        neighbors[i + reuse] = cs.getNeighbors(pU, partMatch[pU], u, neighborCount[i + reuse]);
                    }
                    generateCandidates(t, firstBag, depth, neighbors, neighborCount, parents.size() + reuse,
                                       nodeCandidates[prefixSum[firstBag] + depth], nodeCandCount[prefixSum[firstBag] + depth],
                                       partMatch, cs, u, current->cartesianParent, beginPoses, endPoses);
                    pCandidates[depth] = nodeCandidates[prefixSum[firstBag] + depth];
                    pCandCount[depth] = nodeCandCount[prefixSum[firstBag] + depth];
                }
                else {
                    pCandCount[depth] = cs.candSize(u);
                    memcpy(pCandidates[depth], cs.candSet(u), pCandCount[depth] * sizeof(VertexID));
                }
            }
        }

        --depth;
        if (depth >= mappingSize) {
            const PrefixNode *current = pn;
            if (depth == mappingSize) pn = subtreeRoot;
            else pn = nodes[depth - 1];
            VertexID u = current -> u;
            VertexID v = partMatch[u];
            visited[v] = false;
            for (VertexID nID: current -> nIDsToCall) {
                if (nID == t.numNodes - 1) {
                    bool flag = true;
                    for (VertexID nID2 : current -> nIDsToBuild) {
                        if (localLengths[nID2] == 0) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        for (VertexID nID2 : current -> nIDsToBuild)
                            if (!trieLevels[nID2].oneLevel)
                                buildTrie(localTuples[nID2], trieLevels[nID2], localLengths[nID2]);
                        globalJoin(result, count, query, t, cs, partMatch, mappingSizes[nID], depth + 1, nodeCandidates, nodeCandCount,
                                   prefixSum[nID], visited, globalNode.nIDs, globalNode.attributesBefore,
                                   globalNode.cartesianParent, iters, iterSizes, traverse, neighbors, neighborCount,
                                   traversedLevels, lastLevels, nodePoses, traversePoses, numBranches,
                                   trieLevels[nID], localTuples[nID], localLengths[nID], beginPoses, endPoses);
                    }
                    else {
                        for (VertexID nID2 : current -> nIDsToBuild) {
                            localLengths[nID2] = 0;
                        }
                    }
                }
            }
            if (current -> pathToGlobal) {
                for (VertexID nID : current -> nIDsToJoin) {
                    traversedLevels[nID].pop_back();
                    lastLevels[nID] = traversedLevels[nID].back();
                }
            }
        }
    }
    for (VertexID nID: pn -> nIDsToCall) {
        if (nID == t.numNodes - 1) {
            bool flag = true;
            for (VertexID nID2 : pn -> nIDsToBuild) {
                if (localLengths[nID2] == 0) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                for (VertexID nID2 : pn -> nIDsToBuild)
                    if (!trieLevels[nID2].oneLevel) buildTrie(localTuples[nID2], trieLevels[nID2], localLengths[nID2]);
                globalJoin(result, count, query, t, cs, partMatch, mappingSizes[nID], depth + 1, nodeCandidates, nodeCandCount,
                           prefixSum[nID], visited, globalNode.nIDs, globalNode.attributesBefore,
                           globalNode.cartesianParent, iters, iterSizes, traverse, neighbors, neighborCount,
                           traversedLevels, lastLevels, nodePoses, traversePoses, numBranches,
                           trieLevels[nID], localTuples[nID], localLengths[nID], beginPoses, endPoses);
            }
            else {
                for (VertexID nID2 : pn -> nIDsToBuild) {
//                    tuples[nID2].clear();
                    localLengths[nID2] = 0;
                }
            }
        }
    }
}

void subtreeSharedJoin(const HyperTree &t, const PrefixNode *subtreeRoot, int mappingSize,
                       const std::vector<int> &mappingSizes, const std::map<const PrefixNode *, std::vector<VertexID>> &bagsBelow,
                       const Graph &query, const DataGraph &g, CandidateSpace &cs,
                       std::vector<TrieLevel> &trieLevels, std::vector<std::vector<TrieLevel *>> &traversedLevels,
                       std::vector<TrieLevel *> &lastLevels, bool *visited, VertexID *partMatch,
                       std::vector<std::vector<VertexID>> &result, size_t &count, bool traverse,
                       std::vector<std::vector<std::vector<VertexID>>> &localTuples, std::vector<ui> &localLengths,
                       SubtreeSharedJoinContext &ctx) {
    std::vector<ui> &prefixSum = ctx.prefixSum;
    VertexID **nodeCandidates = ctx.nodeCandidates;
    VertexID **allCandidates = ctx.allCandidates;
    ui *allCandCount = ctx.allCandCount;
    ui *nodeCandCount = ctx.nodeCandCount;
    ui *nodePoses = ctx.nodePoses;
    const HyperNode &globalNode = t.nodes[t.numNodes - 1];
    const PrefixNode *pn = subtreeRoot;
    ui maxSize = cs.getMaxSize();
    ui ***iters = ctx.iters;
    ui *iterSizes = ctx.iterSizes;
    std::vector<ui> &traversePoses = ctx.traversePoses;
    std::vector<ui> &numBranches = ctx.numBranches;
    ui height = ctx.height;
    VertexID **pCandidates = ctx.pCandidates;
    ui *pCandCount = ctx.pCandCount;
    std::fill(pCandCount, pCandCount + height, 0);
    std::vector<ui> &pPoses = ctx.pPoses;
    std::vector<ui> &childPoses = ctx.childPoses;
    std::fill(pPoses.begin(), pPoses.end(), 0);
    std::fill(childPoses.begin(), childPoses.end(), 0);
    const VertexID **neighbors = ctx.neighbors;
    ui *neighborCount = ctx.neighborCount;
    ui maxEdgeSize = maxNumBackWard(t.globalOrder, query);
    ui maxNum = height + globalNode.numAttributes - globalNode.prefixSize;
    std::vector<ui> &beginPoses = ctx.beginPoses;
    std::vector<ui> &endPoses = ctx.endPoses;
    pn = subtreeRoot;
    std::vector<const PrefixNode *> &nodes = ctx.nodes;
    nodes.assign(height, nullptr);
    if (mappingSize > 0) nodes[mappingSize - 1] = subtreeRoot;
    int depth = mappingSize;
    const auto &candIndex = t.candIndex;
    // 处理子树根节点的nIDsToCall
    for (VertexID nID: pn -> nIDsToCall) {
        if (nID != t.numNodes - 1) {
            nodeJoin(t, nID, g, visited, partMatch, mappingSizes[nID], nodeCandidates, nodeCandCount, allCandidates, allCandCount,
                     prefixSum[nID], neighbors, neighborCount, nodePoses, localTuples[nID], trieLevels[nID], localLengths[nID]);
        }
    }

    // 如果子树根有子节点，为第一个子节点创建候选集
    // 由于已经有了匹配过的点，需要利用已匹配的点来计算候选集
    if (!pn -> children.empty()) {
        const PrefixNode *current = pn -> children[0];
        VertexID u = current -> u;
        if (current -> pathToGlobal) {
            bool flag = true;
            for (VertexID nID : pn -> nIDsToBuild) {
                if (localLengths[nID] == 0) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                for (VertexID nID : pn -> nIDsToBuild)
                    if (!trieLevels[nID].oneLevel) buildTrie(localTuples[nID], trieLevels[nID], localLengths[nID]);
            }
            else {
                return;
            }
            VertexID firstBag = bagsBelow.at(current)[0];
            if (!current->nIDsToJoin.empty()) {
                for (int i = 0; i < current->nIDsToJoin.size(); ++i) {
                    VertexID nID = current->nIDsToJoin[i];
                    neighbors[i] = lastLevels[nID]->values;
                    neighborCount[i] = lastLevels[nID]->length;
                }
                for (int i = 0; i < current->attributesBefore.size(); ++i) {
                    VertexID pU = current->attributesBefore[i];
                    VertexID pV = partMatch[pU];
                    ui numNeighbors;
                    const VertexID *pVNeighbors = cs.getNeighbors(pU, pV, u, numNeighbors);
                    neighbors[current->nIDsToJoin.size() + i] = pVNeighbors;
                    neighborCount[current->nIDsToJoin.size() + i] = numNeighbors;
                }
                generateCandidates(t, t.numNodes - 1, depth, neighbors, neighborCount, current->nIDsToJoin.size() + current->attributesBefore.size(),
                                   iters[depth], iterSizes[depth], nodeCandidates[prefixSum[firstBag] + depth], nodeCandCount[prefixSum[firstBag] + depth], current->nIDsToJoin.size(),
                                   partMatch, cs, u, current->cartesianParent, beginPoses, endPoses);
            }
        }
        if (current->nIDsToJoin.empty()) {
            VertexID firstBag = bagsBelow.at(current)[0];
            ui beginPos = 0, endPos = allCandCount[candIndex[firstBag][depth]];
            if (!t.nodes[firstBag].largerAttrs[depth].empty()) {
                VertexID maxCompared = 0;
                for (VertexID u2: t.nodes[firstBag].largerAttrs[depth]) {
                    if (partMatch[u2] > maxCompared)
                        maxCompared = partMatch[u2];
                }
                beginPos = setBeginPos(nodeCandidates[prefixSum[firstBag] + depth], endPos, maxCompared);
            }
            if (!t.nodes[firstBag].smallerAttrs[depth].empty()) {
                VertexID minCompared = (VertexID) - 1;
                for (VertexID u2: t.nodes[firstBag].smallerAttrs[depth]) {
                    if (partMatch[u2] < minCompared)
                        minCompared = partMatch[u2];
                }
                endPos = setEndPos(nodeCandidates[prefixSum[firstBag] + depth], endPos, minCompared);
            }
            pPoses[depth] = beginPos;
            pCandidates[depth] = nodeCandidates[prefixSum[firstBag] + depth];
            pCandCount[depth] = endPos;
            if (current->pathToGlobal) iterSizes[depth] = pCandCount[depth];
            for (int i = 0; i < t.nodes[firstBag].candidatesAfter[depth].size(); ++i) {
                int pos = t.nodes[firstBag].candidatesAfter[depth][i];
                allCandCount[pos] = endPos - beginPos;
                memcpy(allCandidates[pos], nodeCandidates[prefixSum[firstBag] + depth] + beginPos, allCandCount[pos] * sizeof(VertexID));
            }
        }
    }

    while (depth >= mappingSize) {
        while (childPoses[depth] < pn -> children.size()) {
            const PrefixNode *current = pn -> children[childPoses[depth]];
            nodes[depth] = current;
            VertexID u = current -> u;
            bool nextLevel = false;

            if (!current -> pathToGlobal) {
                while (pPoses[depth] < pCandCount[depth]) {
                    VertexID v = pCandidates[depth][pPoses[depth]];
                    ++pPoses[depth];
                    if (visited[v]) continue;
                    bool valid = true;
                    partMatch[u] = v;
                    for (VertexID nID: bagsBelow.at(current)) {
                        if (!generateCandidates(g, v, allCandidates, allCandCount, t.nodes[nID].attributesAfter[depth],
                                                t.nodes[nID].attrAfterLarger[depth], t.nodes[nID].attrAfterSmaller[depth],
                                                t.nodes[nID].copyAfter[depth], t.nodes[nID].copyAfterTypes[depth],
                                                partMatch)) {
                            valid = false;
                            break;
                        }
                    }
                    if (!valid) continue;
                    visited[v] = true;
                    for (VertexID nID: current -> nIDsToCall) {
                        nodeJoin(t, nID, g, visited, partMatch, mappingSizes[nID], nodeCandidates, nodeCandCount, allCandidates, allCandCount,
                                 prefixSum[nID], neighbors, neighborCount, nodePoses, localTuples[nID],
                                 trieLevels[nID], localLengths[nID]);
                    }

                    if (!current -> children.empty()) {
                        nextLevel = true;
                        break;
                    }
                    else visited[v] = false;
                }
            }
            else {
                while (pPoses[depth] < iterSizes[depth]) {
                    VertexID v;
                    if (current->nIDsToJoin.empty()) {
                        v = pCandidates[depth][pPoses[depth]];
                    }
                    else {
                        ui childPos = iters[depth][pPoses[depth]][0];
                        VertexID tmp = current -> nIDsToJoin[0];
                        v = traversedLevels[tmp].back()->values[childPos];
                    }
                    ++pPoses[depth];
                    if (visited[v]) continue;
                    bool valid = true;
                    partMatch[u] = v;
                    for (VertexID nID: bagsBelow.at(current)) {
                        if (!generateCandidates(g, v, allCandidates, allCandCount,
                                                t.nodes[nID].attributesAfter[depth],
                                                t.nodes[nID].attrAfterLarger[depth], t.nodes[nID].attrAfterSmaller[depth],
                                                t.nodes[nID].copyAfter[depth], t.nodes[nID].copyAfterTypes[depth],
                                                partMatch)) {
                            valid = false;
                            break;
                        }
                    }
                    if (!valid) continue;
                    visited[v] = true;
                    for (int i = 0; i < current -> nIDsToJoin.size(); ++i) {
                        VertexID nID = current -> nIDsToJoin[i];
                        ui childPos = iters[depth][pPoses[depth]-1][i];
                        if (traversedLevels[nID].back()->children) {
                            traversedLevels[nID].push_back(&traversedLevels[nID].back()->children[childPos]);
                            lastLevels[nID] = traversedLevels[nID].back();
                        }
                        else {
                            traversedLevels[nID].push_back(nullptr);
                            lastLevels[nID] = nullptr;
                        }
                    }
                    for (VertexID nID: current -> nIDsToCall) {
                        if (nID != t.numNodes - 1) {
                            nodeJoin(t, nID, g, visited, partMatch, mappingSizes[nID], nodeCandidates, nodeCandCount,
                                     allCandidates, allCandCount, prefixSum[nID], neighbors, neighborCount,
                                     nodePoses, localTuples[nID], trieLevels[nID], localLengths[nID]);
                        }
                    }

                    if (!current -> children.empty()) {
                        nextLevel = true;
                        break;
                    }
                    else {
                        for (VertexID nID: current -> nIDsToCall) {
                            if (nID == t.numNodes - 1) {
                                bool flag = true;
                                for (VertexID nID2 : current -> nIDsToBuild) {
                                    if (localLengths[nID2] == 0) {
                                        flag = false;
                                        break;
                                    }
                                }
                                if (flag) {
                                    for (VertexID nID2 : current -> nIDsToBuild)
                                        if (!trieLevels[nID2].oneLevel) buildTrie(localTuples[nID2], trieLevels[nID2], localLengths[nID2]);
                                    globalJoin(result, count, query, t, g, cs, partMatch, mappingSizes[nID],
                                               depth + 1, nodeCandidates, nodeCandCount, allCandidates, allCandCount, prefixSum[nID],
                                               visited, globalNode.nIDs, globalNode.attributesBefore,
                                               globalNode.cartesianParent, iters, iterSizes, traverse, neighbors, neighborCount,
                                               traversedLevels, lastLevels, nodePoses, traversePoses, numBranches,
                                               trieLevels[nID], localTuples[nID], localLengths[nID],
                                               beginPoses, endPoses);
                                }
                                else {
                                    for (VertexID nID2 : current -> nIDsToBuild) {
                                        localLengths[nID2] = 0;
                                    }
                                }
                            }
                        }
                        visited[v] = false;
                        for (VertexID nID : current -> nIDsToJoin) {
                            traversedLevels[nID].pop_back();
                            lastLevels[nID] = traversedLevels[nID].back();
                        }
                    }
                }
            }

            if (!nextLevel) {
                ++childPoses[depth];
                pPoses[depth] = 0;
                if (childPoses[depth] == pn->children.size()) break;
                current = pn->children[childPoses[depth]];
                u = current->u;
            }
            else {
                ++depth;
                pPoses[depth] = 0;
                childPoses[depth] = 0;
                pn = current;
                current = pn -> children[0];
                u = current->u;
            }
            VertexID firstBag = bagsBelow.at(current)[0];
            if (current -> pathToGlobal) {
                bool flag = true;
                for (VertexID nID : pn -> nIDsToBuild) {
                    if (localLengths[nID] == 0) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    for (VertexID nID : pn -> nIDsToBuild)
                        if (!trieLevels[nID].oneLevel) buildTrie(localTuples[nID], trieLevels[nID], localLengths[nID]);
                }
                else {
                    for (VertexID nID2 : pn -> nIDsToBuild) {
                        localLengths[nID2] = 0;
                    }
                    visited[partMatch[pn->u]] = false;
                    for (VertexID nID : pn -> nIDsToJoin) {
                        traversedLevels[nID].pop_back();
                        lastLevels[nID] = traversedLevels[nID].back();
                    }
                    --depth;
                    if (depth < mappingSize) {
                        pn = subtreeRoot;
                        break;
                    }
                    else pn = nodes[depth - 1];
                    continue;
                }
                if (!current->nIDsToJoin.empty()) {
                    for (int i = 0; i < current->nIDsToJoin.size(); ++i) {
                        VertexID nID = current->nIDsToJoin[i];
                        neighbors[i] = lastLevels[nID]->values;
                        neighborCount[i] = lastLevels[nID]->length;
                    }
                    for (int i = 0; i < current->attributesBefore.size(); ++i) {
                        VertexID pU = current->attributesBefore[i];
                        VertexID pV = partMatch[pU];
                        ui numNeighbors;
                        const VertexID *pVNeighbors = cs.getNeighbors(pU, pV, u, numNeighbors);
                        neighbors[current->nIDsToJoin.size() + i] = pVNeighbors;
                        neighborCount[current->nIDsToJoin.size() + i] = numNeighbors;
                    }
                    generateCandidates(t, t.numNodes - 1, depth, neighbors, neighborCount, current->nIDsToJoin.size() + current->attributesBefore.size(),
                                       iters[depth], iterSizes[depth], nodeCandidates[prefixSum[firstBag] + depth], nodeCandCount[prefixSum[firstBag] + depth], current->nIDsToJoin.size(),
                                       partMatch, cs, u, current->cartesianParent,
                                       beginPoses, endPoses);
                }
            }
            if (current->nIDsToJoin.empty()) {
                ui beginPos = 0, endPos = allCandCount[candIndex[firstBag][depth]];
                if (!t.nodes[firstBag].largerAttrs[depth].empty()) {
                    VertexID maxCompared = 0;
                    for (VertexID u2: t.nodes[firstBag].largerAttrs[depth]) {
                        if (partMatch[u2] > maxCompared)
                            maxCompared = partMatch[u2];
                    }
                    beginPos = setBeginPos(nodeCandidates[prefixSum[firstBag] + depth], endPos, maxCompared);
                }
                if (!t.nodes[firstBag].smallerAttrs[depth].empty()) {
                    VertexID minCompared = (VertexID) - 1;
                    for (VertexID u2: t.nodes[firstBag].smallerAttrs[depth]) {
                        if (partMatch[u2] < minCompared)
                            minCompared = partMatch[u2];
                    }
                    endPos = setEndPos(nodeCandidates[prefixSum[firstBag] + depth], endPos, minCompared);
                }
                pPoses[depth] = beginPos;
                pCandidates[depth] = nodeCandidates[prefixSum[firstBag] + depth];
                pCandCount[depth] = endPos;
                if (current->pathToGlobal) iterSizes[depth] = pCandCount[depth];
                for (int i = 0; i < t.nodes[firstBag].candidatesAfter[depth].size(); ++i) {
                    int pos = t.nodes[firstBag].candidatesAfter[depth][i];
                    allCandCount[pos] = endPos - beginPos;
                    memcpy(allCandidates[pos], nodeCandidates[prefixSum[firstBag] + depth] + beginPos, allCandCount[pos] * sizeof(VertexID));
                }
            }
        }

        --depth;
        if (depth >= mappingSize) {
            const PrefixNode *current = pn;
            if (depth == mappingSize) pn = subtreeRoot;
            else pn = nodes[depth - 1];
            VertexID u = current -> u;
            VertexID v = partMatch[u];
            visited[v] = false;
            for (VertexID nID: current -> nIDsToCall) {
                if (nID == t.numNodes - 1) {
                    bool flag = true;
                    for (VertexID nID2 : current -> nIDsToBuild) {
                        if (localLengths[nID2] == 0) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        for (VertexID nID2 : current -> nIDsToBuild)
                            if (!trieLevels[nID2].oneLevel)
                                buildTrie(localTuples[nID2], trieLevels[nID2], localLengths[nID2]);
                        globalJoin(result, count, query, t, g, cs, partMatch, mappingSizes[nID], depth + 1, nodeCandidates, nodeCandCount,
                                   allCandidates, allCandCount,
                                   prefixSum[nID], visited, globalNode.nIDs, globalNode.attributesBefore,
                                   globalNode.cartesianParent, iters, iterSizes, traverse, neighbors, neighborCount,
                                   traversedLevels, lastLevels, nodePoses, traversePoses, numBranches,
                                   trieLevels[nID], localTuples[nID], localLengths[nID], beginPoses, endPoses);
                    }
                    else {
                        for (VertexID nID2 : current -> nIDsToBuild) {
                            localLengths[nID2] = 0;
                        }
                    }
                }
            }
            if (current -> pathToGlobal) {
                for (VertexID nID : current -> nIDsToJoin) {
                    traversedLevels[nID].pop_back();
                    lastLevels[nID] = traversedLevels[nID].back();
                }
            }
        }
    }
    for (VertexID nID: pn -> nIDsToCall) {
        if (nID == t.numNodes - 1) {
            bool flag = true;
            for (VertexID nID2 : pn -> nIDsToBuild) {
                if (localLengths[nID2] == 0) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                for (VertexID nID2 : pn -> nIDsToBuild)
                    if (!trieLevels[nID2].oneLevel) buildTrie(localTuples[nID2], trieLevels[nID2], localLengths[nID2]);
                globalJoin(result, count, query, t, g, cs, partMatch, mappingSizes[nID], depth + 1, nodeCandidates, nodeCandCount,
                           allCandidates, allCandCount,
                           prefixSum[nID], visited, globalNode.nIDs, globalNode.attributesBefore,
                           globalNode.cartesianParent, iters, iterSizes, traverse, neighbors, neighborCount,
                           traversedLevels, lastLevels, nodePoses, traversePoses, numBranches,
                           trieLevels[nID], localTuples[nID], localLengths[nID], beginPoses, endPoses);
            }
            else {
                for (VertexID nID2 : pn -> nIDsToBuild) {
                    localLengths[nID2] = 0;
                }
            }
        }
    }
}

void parNodeJoin(const HyperTree &t, VertexID nID, CandidateSpace &cs, VertexID **allCandidates, ui *allCandCount,
    ui prefixSum, std::vector<std::vector<VertexID>> &tuples, ui &length) {
    const HyperNode &tau = t.nodes[nID];
    VertexID u1 = tau.attributes[0], u2 = tau.attributes[1];
    int type = 0;
    if (std::find(tau.largerAttrs[1].begin(), tau.largerAttrs[1].end(), u1) != tau.largerAttrs[1].end())
        type = 1;
    if (std::find(tau.smallerAttrs[1].begin(), tau.smallerAttrs[1].end(), u1) != tau.smallerAttrs[1].end())
        type = 2;
    const auto &edges = cs.getEdges(type);
    ui maxVertices = cs.getMaxSize();
    ui maxDegree = cs.getMaxDegree();

    // 预计算内存需求
    ui totalNodes = 0;
    for (ui k = 0; k < t.numNodes; ++k) {
        totalNodes += t.nodes[k].numAttributes;
    }

    // 预分配线程内存池
    int maxThreads = omp_get_max_threads();

    // 为每个线程预分配内存
    bool **threadVisited = new bool*[maxThreads];
    VertexID **threadPartMatch = new VertexID*[maxThreads];
    VertexID ***threadAllCandidates = new VertexID**[maxThreads];
    ui **threadAllCandCount = new ui*[maxThreads];
    ui **threadAllPoses = new ui*[maxThreads];
    const VertexID ***threadNeighbors = new const VertexID**[maxThreads];
    ui **threadNeighborCount = new ui*[maxThreads];
    std::vector<std::vector<ui>> threadBeginPoses(maxThreads, std::vector<ui>(t.numAttributes));
    std::vector<std::vector<ui>> threadEndPoses(maxThreads, std::vector<ui>(t.numAttributes));
    std::vector<std::vector<std::vector<VertexID>>> threadTuples(maxThreads);
    std::vector<std::vector<std::vector<VertexID>>> threadResults(maxThreads);
    std::vector<HyperTree> threadTrees(maxThreads);
    if (t.trieOrder[nID].size() > 1) {
        for (int i = 0; i < maxThreads; ++i) {
            threadTuples[i].reserve(1e6);
            threadResults[i].reserve(1e7);
        }
    }
    for (int tid = 0; tid < maxThreads; ++tid) {
        threadVisited[tid] = new bool[maxVertices];
        memset(threadVisited[tid], false, sizeof(bool) * maxVertices);
        threadPartMatch[tid] = new VertexID[t.numAttributes];
        threadAllCandidates[tid] = new VertexID*[totalNodes];
        threadAllCandCount[tid] = new ui[totalNodes];
        threadAllPoses[tid] = new ui[totalNodes];
        threadNeighbors[tid] = new const VertexID*[t.numAttributes];
        threadNeighborCount[tid] = new ui[t.numAttributes];
        for (ui k = 0; k < totalNodes; ++k) {
            threadAllCandidates[tid][k] = allCandidates[k];
            threadAllCandCount[tid][k] = allCandCount[k];
            threadAllPoses[tid][k] = 0;
        }
        // 为第3个及以后的属性预分配空间
        for (int i = 2; i < tau.numAttributes; ++i) {
            threadAllCandidates[tid][prefixSum + i] = new VertexID[maxDegree];
            threadAllCandCount[tid][prefixSum + i] = 0;
        }
        threadTrees[tid] = t;
    }
    TrieLevel localLevel;
    localLevel.oneLevel = true;
    // 为前两个顶点的每个匹配创建一个任务
#pragma omp parallel
    {
        int tid = omp_get_thread_num();

        // 获取当前线程的预分配内存
        bool *localVisited = threadVisited[tid];
        VertexID *localPartMatch = threadPartMatch[tid];
        VertexID **localAllCandidates = threadAllCandidates[tid];
        ui *localAllCandCount = threadAllCandCount[tid];
        ui *localAllPoses = threadAllPoses[tid];
        const VertexID **localNeighbors = threadNeighbors[tid];
        ui *localNeighborCount = threadNeighborCount[tid];
        std::vector<ui> &localBeginPoses = threadBeginPoses[tid];
        std::vector<ui> &localEndPoses = threadEndPoses[tid];
        std::vector<std::vector<VertexID>> &localTuples = threadTuples[tid];
        ui localLength = 0;
        // 线程级别的结果累积容器
        std::vector<std::vector<VertexID>> &localResults = threadResults[tid];
        // 使用动态或静态调度并行处理边
#pragma omp for schedule(dynamic) nowait
        for (int i = 0; i < edges.size(); ++i) {
            VertexID v1 = edges[i].first;
            VertexID v2 = edges[i].second;
            localPartMatch[u1] = v1;
            localPartMatch[u2] = v2;
            localVisited[v1] = true;
            localVisited[v2] = true;
            localLength = 0;
            // 调用nodeJoin处理从第三个顶点开始的匹配
            nodeJoin(threadTrees[tid], nID, cs, localVisited, localPartMatch, 2,
                     localAllCandidates, localAllCandCount, prefixSum,
                     localNeighbors, localNeighborCount, localAllPoses,
                     localTuples, localLevel, localLength, localBeginPoses, localEndPoses);
            localVisited[v1] = false;
            localVisited[v2] = false;

            // 将本次迭代的结果累积到线程级容器中
            for (ui j = 0; j < localLength && j < localTuples.size(); ++j) {
                localResults.push_back(localTuples[j]);
            }
        }


        // 循环结束后，一次性合并所有结果到全局容器
        if (!localResults.empty()) {
            #pragma omp critical(merge_parnodejoin_results)
            {
                for (auto &&result : localResults) {
                    if (length < tuples.size()) {
                        tuples[length] = std::move(result);
                    } else {
                        tuples.push_back(std::move(result));
                    }
                    ++length;
                }
            }
        }
    }

    // 清理预分配的内存池
    for (int tid = 0; tid < maxThreads; ++tid) {
        delete[] threadVisited[tid];
        delete[] threadPartMatch[tid];
        for (int i = 2; i < tau.numAttributes; ++i) {
            delete[] threadAllCandidates[tid][prefixSum + i];
        }
        delete[] threadAllCandidates[tid];
        delete[] threadAllCandCount[tid];
        delete[] threadAllPoses[tid];
        delete[] threadNeighbors[tid];
        delete[] threadNeighborCount[tid];
    }

    delete[] threadVisited;
    delete[] threadPartMatch;
    delete[] threadAllCandidates;
    delete[] threadAllCandCount;
    delete[] threadAllPoses;
    delete[] threadNeighbors;
    delete[] threadNeighborCount;
}

void parSharedJoin(const HyperTree &t, const PrefixNode *pt, const Graph &query, CandidateSpace &cs,
    std::vector<TrieLevel> &trieLevels, bool *visited, std::vector<std::vector<VertexID>> &result,
    size_t &count, bool traverse) {
    std::vector<std::vector<std::vector<VertexID>>> tuples(t.numNodes);
    std::vector<ui> lengths(t.numNodes, 0);
    for (int i = 0; i < t.numNodes; ++i) {
        if (t.nodes[i].numAttributes > t.nodes[i].prefixSize + 1)
            tuples[i].reserve(1e7);
    }
    std::vector<ui> prefixSum(t.numNodes + 1, 0);
    for (VertexID nID = 1; nID <= t.numNodes; ++nID) {
        prefixSum[nID] = prefixSum[nID - 1] + t.nodes[nID - 1].numAttributes;
    }
    VertexID *partMatch = new VertexID[t.numAttributes];
    VertexID **nodeCandidates = new VertexID *[prefixSum.back()];
    ui *nodeCandCount = new ui [prefixSum.back()];
    ui *nodePoses = new ui [prefixSum.back()];
    for (int i = 0; i < prefixSum.back(); ++i) {
        nodeCandidates[i] = new VertexID[cs.getMaxSize()];
        nodeCandCount[i] = 0;
        nodePoses[i] = 0;
    }
    ui height = pt->getHeight();
    VertexID **pCandidates = new VertexID *[height];
    ui *pCandCount = new ui[height];
    pCandidates[0] = new VertexID[cs.getMaxSize()];
    for (int i = 0; i < height; ++i) {
        pCandCount[i] = 0;
    }
    const VertexID **neighbors = new const VertexID *[t.numAttributes];
    ui *neighborCount = new ui[t.numAttributes];
    ui maxVertices = cs.getMaxSize();
    ui maxDegree = cs.getMaxDegree();
    std::vector<std::vector<TrieLevel *>> traversedLevels(trieLevels.size());
    std::vector<TrieLevel *> lastLevels(trieLevels.size());
    for (VertexID nID = 0; nID < trieLevels.size(); ++nID) {
        traversedLevels[nID].push_back(&trieLevels[nID]);
        lastLevels[nID] = &trieLevels[nID];
    }
    // 计算mappingSizes和bagsBelow
    std::vector<int> mappingSizes = getMappingSizes(t, pt);
    std::map<const PrefixNode *, std::vector<VertexID>> bagsBelow;
    pt->addBagsBelow(bagsBelow, t.numNodes - 1);
    // collect nIDs whose prefix size is 0
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (t.nodes[nID].prefixSize == 0) {
            parNodeJoin(t, nID, cs, nodeCandidates, nodeCandCount, prefixSum[nID],
                        tuples[nID], lengths[nID]);
            buildTrie(tuples[nID], trieLevels[nID], lengths[nID]);
        }
    }
    PrefixNode *firstShared = pt->children.back();
    if (!(firstShared->children.size() == 1 && firstShared->nIDsToCall.empty())) {
        VertexID u = firstShared->u;
        std::vector<VertexID> candidates;
        candidates.reserve(maxVertices + 1);
        for (VertexID v = 0; v < maxVertices; ++v) {
            bool valid = true;
            for (VertexID nID: firstShared->nIDsToJoin) {
                const VertexID *values = traversedLevels[nID].back()->values;
                ui pos = ComputeSetIntersection::BinarySearch(values, 0, traversedLevels[nID].back()->length, v);
                if (pos == traversedLevels[nID].back()->length || values[pos] != v) {
                    valid = false;
                    break;
                }
            }
            if (valid) {
                candidates.push_back(v);
            }
        }
        int maxThreads = omp_get_max_threads();
        bool **threadVisited = new bool *[maxThreads];
        VertexID **threadPartMatch = new VertexID*[maxThreads];
        std::vector<SubtreeSharedJoinContext> threadContexts(maxThreads);
        std::vector<std::vector<TrieLevel>> threadTrieLevels(maxThreads, trieLevels);
        std::vector<std::vector<std::vector<TrieLevel *>>> threadTraversedLevels(maxThreads);
        std::vector<std::vector<TrieLevel *>> threadLastLevels(maxThreads);
        std::vector<std::vector<std::vector<std::vector<VertexID>>>> threadTuples(maxThreads);
        std::vector<std::vector<ui>> threadLengths(maxThreads);
        std::vector<HyperTree> threadTrees(maxThreads);
        for (int tid = 0; tid < maxThreads; ++tid) {
            threadVisited[tid] = new bool[maxVertices];
            memset(threadVisited[tid], false, sizeof(bool) * maxVertices);
            threadPartMatch[tid] = new VertexID[t.numAttributes];
            threadTuples[tid].resize(t.numNodes);
            threadTraversedLevels[tid].resize(t.numNodes);
            threadLastLevels[tid].resize(t.numNodes);
            for (VertexID nID = 0; nID < t.numNodes; ++nID) {
                threadTraversedLevels[tid].reserve(t.numAttributes);
                threadTraversedLevels[tid][nID].push_back(&threadTrieLevels[tid][nID]);
                threadLastLevels[tid][nID] = &threadTrieLevels[tid][nID];
                threadLengths[tid] = std::vector<ui>(t.numNodes, 0);
            }
            for (VertexID nID: bagsBelow[firstShared]) {
                if (t.trieOrder[nID].size() > 1) threadTuples[tid][nID].reserve(1e6);
            }
            initSubtreeSharedJoinContext(threadContexts[tid], t, cs, firstShared, 1);
            threadTrees[tid] = t;
        }
        // 为每个数据图顶点创建一个任务
#pragma omp parallel reduction(+ : count)
        {
            int tid = omp_get_thread_num();
            bool *localVisited = threadVisited[tid];
            memset(localVisited, false, sizeof(bool) * maxVertices);
            VertexID *localPartMatch = threadPartMatch[tid];
            // 初始化线程私有的tuples和trieLevels
            std::vector<TrieLevel> &localTrieLevels = threadTrieLevels[tid];
            std::vector<std::vector<TrieLevel *>> &localTraversedLevels = threadTraversedLevels[tid];
            std::vector<TrieLevel *> &localLastLevels = threadLastLevels[tid];
            std::vector<std::vector<std::vector<VertexID>>> &localTuples = threadTuples[tid];
            std::vector<ui> &localLengths = threadLengths[tid];
            std::vector<std::vector<VertexID>> localResult;
            // 创建子树共享连接上下文
            SubtreeSharedJoinContext ctx = threadContexts[tid];
            
            // 使用动态或静态调度为每个数据图顶点执行subtreeSharedJoin

#pragma omp for schedule(dynamic) nowait
            for (size_t i = 0; i < candidates.size(); ++i) {
                VertexID v = candidates[i];
                for (VertexID nID: firstShared->nIDsToJoin) {
                    const VertexID *values = traversedLevels[nID].back()->values;
                    ui pos = ComputeSetIntersection::BinarySearch(values, 0, traversedLevels[nID].back()->length, v);
                    localTraversedLevels[nID].push_back(&localTrieLevels[nID].children[pos]);
                    localLastLevels[nID] = &localTrieLevels[nID].children[pos];
                }
                // 设置局部匹配和访问状态
                localPartMatch[u] = v;
                localVisited[v] = true;
                size_t localCount = 0;

                // 执行subtreeSharedJoin
                subtreeSharedJoin(threadTrees[tid], firstShared, 1, mappingSizes, bagsBelow, query, cs, localTrieLevels,
                                  localTraversedLevels, localLastLevels,
                                  localVisited, localPartMatch, localResult, localCount, traverse,
                                  localTuples, localLengths, ctx);

                count += localCount;
                // 重置访问状态
                localVisited[v] = false;
                for (VertexID nID: firstShared->nIDsToJoin) {
                    localTraversedLevels[nID].pop_back();
                    localLastLevels[nID] = localTraversedLevels[nID].back();
                }
            }
            
            // 释放内存
            for (VertexID nID = 0; nID < t.numNodes; ++nID) {
                localTrieLevels[nID].children = nullptr;
                localTrieLevels[nID].values = nullptr;
            }
        }
        for (int tid = 0; tid < maxThreads; ++tid) {
            delete[] threadVisited[tid];
            delete[] threadPartMatch[tid];
            freeSubtreeSharedJoinContext(threadContexts[tid], t, cs, firstShared);
        }
    }
    else {
        VertexID u1 = firstShared->u;
        PrefixNode *secondShared = firstShared->children.back();
        VertexID u2 = secondShared->u;
        int type = 0;
        for (auto &rule: t.symmetryRules) {
            if (rule[0] == u1) {
                for (VertexID u3: rule) {
                    if (u3 == u2) {
                        type = 1;
                        break;
                    }
                }
            }
            if (rule[0] == u2) {
                for (VertexID u3: rule) {
                    if (u3 == u1) {
                        type = 2;
                        break;
                    }
                }
            }
        }
        const auto &edges = cs.getEdges(type);
        std::vector<std::pair<VertexID, VertexID>> candidates;
        candidates.reserve(edges.size());
        for (auto &edge: edges) {
            VertexID v1 = edge.first, v2 = edge.second;
            bool valid = true;
            for (VertexID nID: firstShared->nIDsToJoin) {
                const VertexID *values = traversedLevels[nID].back()->values;
                ui pos = ComputeSetIntersection::BinarySearch(values, 0, traversedLevels[nID].back()->length, v1);
                if (pos == traversedLevels[nID].back()->length || values[pos] != v1) {
                    valid = false;
                    break;
                }
            }
            if (!valid) {
                continue;
            }
            for (VertexID nID: secondShared->nIDsToJoin) {
                const VertexID *values = traversedLevels[nID].back()->values;
                ui pos = ComputeSetIntersection::BinarySearch(values, 0, traversedLevels[nID].back()->length, v2);
                if (pos == traversedLevels[nID].back()->length || values[pos] != v2) {
                    valid = false;
                    break;
                }
            }
            if (!valid) {
                continue;
            }
            candidates.push_back(edge);
        }
        int maxThreads = omp_get_max_threads();
        bool **threadVisited = new bool *[maxThreads];
        VertexID **threadPartMatch = new VertexID*[maxThreads];
        std::vector<SubtreeSharedJoinContext> threadContexts(maxThreads);
        std::vector<std::vector<TrieLevel>> threadTrieLevels(maxThreads, trieLevels);
        std::vector<std::vector<std::vector<TrieLevel *>>> threadTraversedLevels(maxThreads);
        std::vector<std::vector<TrieLevel *>> threadLastLevels(maxThreads);
        std::vector<std::vector<std::vector<std::vector<VertexID>>>> threadTuples(maxThreads);
        std::vector<std::vector<ui>> threadLengths(maxThreads);
        std::vector<HyperTree> threadTrees(maxThreads);
        for (int tid = 0; tid < maxThreads; ++tid) {
            threadVisited[tid] = new bool[maxVertices];
            memset(threadVisited[tid], false, sizeof(bool) * maxVertices);
            threadPartMatch[tid] = new VertexID[t.numAttributes];
            threadTuples[tid].resize(t.numNodes);
            threadTraversedLevels[tid].resize(t.numNodes);
            threadLastLevels[tid].resize(t.numNodes);
            for (VertexID nID = 0; nID < t.numNodes; ++nID) {
                threadTraversedLevels[tid].reserve(t.numAttributes);
                threadTraversedLevels[tid][nID].push_back(&threadTrieLevels[tid][nID]);
                threadLastLevels[tid][nID] = &threadTrieLevels[tid][nID];
                threadLengths[tid] = std::vector<ui>(t.numNodes, 0);
            }
            for (VertexID nID: bagsBelow[secondShared]) {
                if (t.trieOrder[nID].size() > 1) threadTuples[tid][nID].reserve(1e6);
            }
            initSubtreeSharedJoinContext(threadContexts[tid], t, cs, secondShared, 2);
            threadTrees[tid] = t;
        }
// 为每个数据图边创建一个任务
#pragma omp parallel reduction(+ : count)
        {
            int tid = omp_get_thread_num();
            bool *localVisited = threadVisited[tid];
            memset(localVisited, false, sizeof(bool) * maxVertices);
            VertexID *localPartMatch = threadPartMatch[tid];
            // 初始化线程私有的tuples和trieLevels
            std::vector<TrieLevel> &localTrieLevels = threadTrieLevels[tid];
            std::vector<std::vector<TrieLevel *>> &localTraversedLevels = threadTraversedLevels[tid];
            std::vector<TrieLevel *> &localLastLevels = threadLastLevels[tid];
            std::vector<std::vector<std::vector<VertexID>>> &localTuples = threadTuples[tid];
            std::vector<ui> &localLengths = threadLengths[tid];
            std::vector<std::vector<VertexID>> localResult;
            // 创建子树共享连接上下文
            SubtreeSharedJoinContext ctx = threadContexts[tid];
            
            // 使用动态或静态调度为每个数据图边执行subtreeSharedJoin

#pragma omp for schedule(dynamic) nowait
            for (int i = 0; i < candidates.size(); ++i) {
                VertexID v1 = candidates[i].first;
                VertexID v2 = candidates[i].second;
                for (VertexID nID: firstShared->nIDsToJoin) {
                    const VertexID *values = traversedLevels[nID].back()->values;
                    ui pos = ComputeSetIntersection::BinarySearch(values, 0, traversedLevels[nID].back()->length, v1);
                    localTraversedLevels[nID].push_back(&localTrieLevels[nID].children[pos]);
                    localLastLevels[nID] = &localTrieLevels[nID].children[pos];
                }
                for (VertexID nID: secondShared->nIDsToJoin) {
                    const VertexID *values = traversedLevels[nID].back()->values;
                    ui pos = ComputeSetIntersection::BinarySearch(values, 0, traversedLevels[nID].back()->length, v2);
                    localTraversedLevels[nID].push_back(&localTrieLevels[nID].children[pos]);
                    localLastLevels[nID] = &localTrieLevels[nID].children[pos];
                }

                // 设置局部匹配和访问状态
                localPartMatch[u1] = v1;
                localPartMatch[u2] = v2;
                localVisited[v1] = true;
                localVisited[v2] = true;

                size_t localCount = 0;

                // 执行subtreeSharedJoin
                subtreeSharedJoin(threadTrees[tid], secondShared, 2, mappingSizes, bagsBelow, query, cs, localTrieLevels,
                                  localTraversedLevels, localLastLevels,
                                  localVisited, localPartMatch, localResult, localCount, traverse,
                                  localTuples, localLengths, ctx);

                count += localCount;

                // 重置访问状态
                localVisited[v1] = false;
                localVisited[v2] = false;
                for (VertexID nID: firstShared->nIDsToJoin) {
                    localTraversedLevels[nID].pop_back();
                    localLastLevels[nID] = localTraversedLevels[nID].back();
                }
                for (VertexID nID: secondShared->nIDsToJoin) {
                    localTraversedLevels[nID].pop_back();
                    localLastLevels[nID] = localTraversedLevels[nID].back();
                }
            }
            
            // 释放内存
            for (VertexID nID = 0; nID < t.numNodes; ++nID) {
                localTrieLevels[nID].children = nullptr;
                localTrieLevels[nID].values = nullptr;
            }
        }
        for (int tid = 0; tid < maxThreads; ++tid) {
            delete[] threadVisited[tid];
            delete[] threadPartMatch[tid];
            freeSubtreeSharedJoinContext(threadContexts[tid], t, cs, secondShared);
        }
    }
    
    // 释放内存
    delete[] partMatch;
    for (int i = 0; i < prefixSum.back(); ++i) {
        delete[] nodeCandidates[i];
    }
    delete[] nodeCandidates;
    delete[] nodeCandCount;
    delete[] nodePoses;
    delete[] pCandidates[0];
    delete[] pCandidates;
    delete[] pCandCount;
    delete[] neighbors;
    delete[] neighborCount;
}

void parSharedJoin(const HyperTree &t, const PrefixNode *pt, const Graph &query, const DataGraph &g, CandidateSpace &cs,
                   std::vector<TrieLevel> &trieLevels, bool *visited, std::vector<std::vector<VertexID>> &result,
                   size_t &count, bool traverse) {
    std::vector<std::vector<std::vector<VertexID>>> tuples(t.numNodes);
    std::vector<ui> lengths(t.numNodes, 0);
    for (int i = 0; i < t.numNodes; ++i) {
        if (t.nodes[i].numAttributes > t.nodes[i].prefixSize + 1)
            tuples[i].reserve(1e7);
    }
    std::vector<ui> prefixSum(t.numNodes + 1, 0);
    for (VertexID nID = 1; nID <= t.numNodes; ++nID) {
        prefixSum[nID] = prefixSum[nID - 1] + t.nodes[nID - 1].numAttributes;
    }
    VertexID *partMatch = new VertexID[t.numAttributes];
    VertexID **nodeCandidates = new VertexID *[prefixSum.back()];
    ui *nodeCandCount = new ui [prefixSum.back()];
    ui *nodePoses = new ui [prefixSum.back()];
    for (int i = 0; i < prefixSum.back(); ++i) {
        nodeCandidates[i] = new VertexID[cs.getMaxSize()];
        nodeCandCount[i] = 0;
        nodePoses[i] = 0;
    }
    ui height = pt->getHeight();
    VertexID **pCandidates = new VertexID *[height];
    ui *pCandCount = new ui[height];
    pCandidates[0] = new VertexID[cs.getMaxSize()];
    for (int i = 0; i < height; ++i) {
        pCandCount[i] = 0;
    }
    const VertexID **neighbors = new const VertexID *[t.numAttributes];
    ui *neighborCount = new ui[t.numAttributes];
    ui maxVertices = cs.getMaxSize();
    ui maxDegree = cs.getMaxDegree();
    std::vector<std::vector<TrieLevel *>> traversedLevels(trieLevels.size());
    std::vector<TrieLevel *> lastLevels(trieLevels.size());
    for (VertexID nID = 0; nID < trieLevels.size(); ++nID) {
        traversedLevels[nID].push_back(&trieLevels[nID]);
        lastLevels[nID] = &trieLevels[nID];
    }
    // 计算mappingSizes和bagsBelow
    std::vector<int> mappingSizes = getMappingSizes(t, pt);
    std::map<const PrefixNode *, std::vector<VertexID>> bagsBelow;
    pt->addBagsBelow(bagsBelow, t.numNodes - 1);
//    parNodeJoin(t, 0, cs, nodeCandidates, nodeCandCount, prefixSum[0], tuples[0], lengths[0]);
//    return;
    // collect nIDs whose prefix size is 0
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (t.nodes[nID].prefixSize == 0) {
            parNodeJoin(t, nID, cs, nodeCandidates, nodeCandCount, prefixSum[nID],
                        tuples[nID], lengths[nID]);
            buildTrie(tuples[nID], trieLevels[nID], lengths[nID]);
        }
    }
    PrefixNode *firstShared = pt->children.back();
    if (!(firstShared->children.size() == 1 && firstShared->nIDsToCall.empty())) {
        VertexID u = firstShared->u;
        std::vector<VertexID> candidates;
        candidates.reserve(maxVertices + 1);
        for (VertexID v = 0; v < maxVertices; ++v) {
            bool valid = true;
            for (VertexID nID: firstShared->nIDsToJoin) {
                const VertexID *values = traversedLevels[nID].back()->values;
                ui pos = ComputeSetIntersection::BinarySearch(values, 0, traversedLevels[nID].back()->length, v);
                if (pos == traversedLevels[nID].back()->length || values[pos] != v) {
                    valid = false;
                    break;
                }
            }
            if (valid) {
                candidates.push_back(v);
            }
        }
        int maxThreads = omp_get_max_threads();
        bool **threadVisited = new bool *[maxThreads];
        VertexID **threadPartMatch = new VertexID*[maxThreads];
        std::vector<SubtreeSharedJoinContext> threadContexts(maxThreads);
        std::vector<std::vector<TrieLevel>> threadTrieLevels(maxThreads, trieLevels);
        std::vector<std::vector<std::vector<TrieLevel *>>> threadTraversedLevels(maxThreads);
        std::vector<std::vector<TrieLevel *>> threadLastLevels(maxThreads);
        std::vector<std::vector<std::vector<std::vector<VertexID>>>> threadTuples(maxThreads);
        std::vector<std::vector<ui>> threadLengths(maxThreads);
        std::vector<HyperTree> threadTrees(maxThreads);
        for (int tid = 0; tid < maxThreads; ++tid) {
            threadVisited[tid] = new bool[maxVertices];
            memset(threadVisited[tid], false, sizeof(bool) * maxVertices);
            threadPartMatch[tid] = new VertexID[t.numAttributes];
            threadTuples[tid].resize(t.numNodes);
            threadTraversedLevels[tid].resize(t.numNodes);
            threadLastLevels[tid].resize(t.numNodes);
            for (VertexID nID = 0; nID < t.numNodes; ++nID) {
                threadTraversedLevels[tid].reserve(t.numAttributes);
                threadTraversedLevels[tid][nID].push_back(&threadTrieLevels[tid][nID]);
                threadLastLevels[tid][nID] = &threadTrieLevels[tid][nID];
                threadLengths[tid] = std::vector<ui>(t.numNodes, 0);
            }
            for (VertexID nID: bagsBelow[firstShared]) {
                if (t.trieOrder[nID].size() > 1) threadTuples[tid][nID].reserve(1e6);
            }
            initSubtreeSharedJoinContext(threadContexts[tid], t, cs, firstShared, 1);
            threadTrees[tid] = t;
        }
        // 为每个数据图顶点创建一个任务
#pragma omp parallel reduction(+ : count)
        {
            int tid = omp_get_thread_num();
            bool *localVisited = threadVisited[tid];
            memset(localVisited, false, sizeof(bool) * maxVertices);
            VertexID *localPartMatch = threadPartMatch[tid];
            // 初始化线程私有的tuples和trieLevels
            std::vector<TrieLevel> &localTrieLevels = threadTrieLevels[tid];
            std::vector<std::vector<TrieLevel *>> &localTraversedLevels = threadTraversedLevels[tid];
            std::vector<TrieLevel *> &localLastLevels = threadLastLevels[tid];
            std::vector<std::vector<std::vector<VertexID>>> &localTuples = threadTuples[tid];
            std::vector<ui> &localLengths = threadLengths[tid];
            std::vector<std::vector<VertexID>> localResult;
            // 创建子树共享连接上下文
            SubtreeSharedJoinContext ctx = threadContexts[tid];

            // 使用动态或静态调度为每个数据图顶点执行subtreeSharedJoin

#pragma omp for schedule(dynamic) nowait
            for (size_t i = 0; i < candidates.size(); ++i) {
                VertexID v = candidates[i];
                localPartMatch[u] = v;
                for (VertexID nID: bagsBelow[firstShared]) {
                    generateCandidates(g, v, ctx.allCandidates, ctx.allCandCount,
                                       t.nodes[nID].attributesAfter[0],
                                       t.nodes[nID].attrAfterLarger[0], t.nodes[nID].attrAfterSmaller[0],
                                       t.nodes[nID].copyAfter[0], t.nodes[nID].copyAfterTypes[0], localPartMatch);
                }
                for (VertexID nID: firstShared->nIDsToJoin) {
                    const VertexID *values = traversedLevels[nID].back()->values;
                    ui pos = ComputeSetIntersection::BinarySearch(values, 0, traversedLevels[nID].back()->length, v);
                    localTraversedLevels[nID].push_back(&localTrieLevels[nID].children[pos]);
                    localLastLevels[nID] = &localTrieLevels[nID].children[pos];
                }
                // 设置局部匹配和访问状态

                localVisited[v] = true;
                size_t localCount = 0;

                // 执行subtreeSharedJoin
                subtreeSharedJoin(threadTrees[tid], firstShared, 1, mappingSizes, bagsBelow, query, g, cs, localTrieLevels,
                                  localTraversedLevels, localLastLevels,
                                  localVisited, localPartMatch, localResult, localCount, traverse,
                                  localTuples, localLengths, ctx);

                count += localCount;
                // 重置访问状态
                localVisited[v] = false;
                for (VertexID nID: firstShared->nIDsToJoin) {
                    localTraversedLevels[nID].pop_back();
                    localLastLevels[nID] = localTraversedLevels[nID].back();
                }
            }

            // 释放内存
            for (VertexID nID = 0; nID < t.numNodes; ++nID) {
                localTrieLevels[nID].children = nullptr;
                localTrieLevels[nID].values = nullptr;
            }
        }
        for (int tid = 0; tid < maxThreads; ++tid) {
            delete[] threadVisited[tid];
            delete[] threadPartMatch[tid];
            freeSubtreeSharedJoinContext(threadContexts[tid], t, cs, firstShared);
        }
    }
    else {
        VertexID u1 = firstShared->u;
        PrefixNode *secondShared = firstShared->children.back();
        VertexID u2 = secondShared->u;
        int type = 0;
        for (auto &rule: t.symmetryRules) {
            if (rule[0] == u1) {
                for (VertexID u3: rule) {
                    if (u3 == u2) {
                        type = 1;
                        break;
                    }
                }
            }
            if (rule[0] == u2) {
                for (VertexID u3: rule) {
                    if (u3 == u1) {
                        type = 2;
                        break;
                    }
                }
            }
        }
        const auto &edges = cs.getEdges(type);
        std::vector<std::pair<VertexID, VertexID>> candidates;
        candidates.reserve(edges.size());
        for (auto &edge: edges) {
            VertexID v1 = edge.first, v2 = edge.second;
            bool valid = true;
            for (VertexID nID: firstShared->nIDsToJoin) {
                const VertexID *values = traversedLevels[nID].back()->values;
                ui pos = ComputeSetIntersection::BinarySearch(values, 0, traversedLevels[nID].back()->length, v1);
                if (pos == traversedLevels[nID].back()->length || values[pos] != v1) {
                    valid = false;
                    break;
                }
            }
            if (!valid) {
                continue;
            }
            for (VertexID nID: secondShared->nIDsToJoin) {
                const VertexID *values = traversedLevels[nID].back()->values;
                ui pos = ComputeSetIntersection::BinarySearch(values, 0, traversedLevels[nID].back()->length, v2);
                if (pos == traversedLevels[nID].back()->length || values[pos] != v2) {
                    valid = false;
                    break;
                }
            }
            if (!valid) {
                continue;
            }
            candidates.push_back(edge);
        }
        int maxThreads = omp_get_max_threads();
        bool **threadVisited = new bool *[maxThreads];
        VertexID **threadPartMatch = new VertexID*[maxThreads];
        std::vector<SubtreeSharedJoinContext> threadContexts(maxThreads);
        std::vector<std::vector<TrieLevel>> threadTrieLevels(maxThreads, trieLevels);
        std::vector<std::vector<std::vector<TrieLevel *>>> threadTraversedLevels(maxThreads);
        std::vector<std::vector<TrieLevel *>> threadLastLevels(maxThreads);
        std::vector<std::vector<std::vector<std::vector<VertexID>>>> threadTuples(maxThreads);
        std::vector<std::vector<ui>> threadLengths(maxThreads);
        std::vector<HyperTree> threadTrees(maxThreads);
        for (int tid = 0; tid < maxThreads; ++tid) {
            threadVisited[tid] = new bool[maxVertices];
            memset(threadVisited[tid], false, sizeof(bool) * maxVertices);
            threadPartMatch[tid] = new VertexID[t.numAttributes];
            threadTuples[tid].resize(t.numNodes);
            threadTraversedLevels[tid].resize(t.numNodes);
            threadLastLevels[tid].resize(t.numNodes);
            for (VertexID nID = 0; nID < t.numNodes; ++nID) {
                threadTraversedLevels[tid].reserve(t.numAttributes);
                threadTraversedLevels[tid][nID].push_back(&threadTrieLevels[tid][nID]);
                threadLastLevels[tid][nID] = &threadTrieLevels[tid][nID];
                threadLengths[tid] = std::vector<ui>(t.numNodes, 0);
            }
            for (VertexID nID: bagsBelow[secondShared]) {
                if (t.trieOrder[nID].size() > 1) threadTuples[tid][nID].reserve(1e6);
            }
            initSubtreeSharedJoinContext(threadContexts[tid], t, cs, secondShared, 2);
            threadTrees[tid] = t;
        }
// 为每个数据图边创建一个任务
#pragma omp parallel reduction(+ : count)
        {
            int tid = omp_get_thread_num();
            bool *localVisited = threadVisited[tid];
            memset(localVisited, false, sizeof(bool) * maxVertices);
            VertexID *localPartMatch = threadPartMatch[tid];
            // 初始化线程私有的tuples和trieLevels
            std::vector<TrieLevel> &localTrieLevels = threadTrieLevels[tid];
            std::vector<std::vector<TrieLevel *>> &localTraversedLevels = threadTraversedLevels[tid];
            std::vector<TrieLevel *> &localLastLevels = threadLastLevels[tid];
            std::vector<std::vector<std::vector<VertexID>>> &localTuples = threadTuples[tid];
            std::vector<ui> &localLengths = threadLengths[tid];
            std::vector<std::vector<VertexID>> localResult;
            // 创建子树共享连接上下文
            SubtreeSharedJoinContext ctx = threadContexts[tid];
            // 使用动态或静态调度为每个数据图边执行subtreeSharedJoin
#pragma omp for schedule(dynamic) nowait
            for (int i = 0; i < candidates.size(); ++i) {
                VertexID v1 = candidates[i].first;
                VertexID v2 = candidates[i].second;
                localPartMatch[u1] = v1;
                localPartMatch[u2] = v2;
                for (VertexID nID: bagsBelow[firstShared]) {
                    generateCandidates(g, v1, ctx.allCandidates, ctx.allCandCount,
                                       t.nodes[nID].attributesAfter[0],
                                       t.nodes[nID].attrAfterLarger[0], t.nodes[nID].attrAfterSmaller[0],
                                       t.nodes[nID].copyAfter[0], t.nodes[nID].copyAfterTypes[0],
                                       localPartMatch);
                }
                bool valid = true;
                for (VertexID nID: bagsBelow[secondShared]) {
                    if (!generateCandidates(g, v2, ctx.allCandidates, ctx.allCandCount,
                                            t.nodes[nID].attributesAfter[1],
                                           t.nodes[nID].attrAfterLarger[1], t.nodes[nID].attrAfterSmaller[1],
                                           t.nodes[nID].copyAfter[1],
                                           t.nodes[nID].copyAfterTypes[1], localPartMatch)) {
                        valid = false;
                        break;
                    }
                }
                if (!valid) continue;
                for (VertexID nID: firstShared->nIDsToJoin) {
                    const VertexID *values = traversedLevels[nID].back()->values;
                    ui pos = ComputeSetIntersection::BinarySearch(values, 0, traversedLevels[nID].back()->length, v1);
                    localTraversedLevels[nID].push_back(&localTrieLevels[nID].children[pos]);
                    localLastLevels[nID] = &localTrieLevels[nID].children[pos];
                }
                for (VertexID nID: secondShared->nIDsToJoin) {
                    const VertexID *values = traversedLevels[nID].back()->values;
                    ui pos = ComputeSetIntersection::BinarySearch(values, 0, traversedLevels[nID].back()->length, v2);
                    localTraversedLevels[nID].push_back(&localTrieLevels[nID].children[pos]);
                    localLastLevels[nID] = &localTrieLevels[nID].children[pos];
                }

                // 设置局部匹配和访问状态
                localVisited[v1] = true;
                localVisited[v2] = true;

                size_t localCount = 0;
                // 执行subtreeSharedJoin
                subtreeSharedJoin(threadTrees[tid], secondShared, 2, mappingSizes, bagsBelow, query, g, cs, localTrieLevels,
                                  localTraversedLevels, localLastLevels,
                                  localVisited, localPartMatch, localResult, localCount, traverse,
                                  localTuples, localLengths, ctx);

                count += localCount;

                // 重置访问状态
                localVisited[v1] = false;
                localVisited[v2] = false;
                for (VertexID nID: firstShared->nIDsToJoin) {
                    localTraversedLevels[nID].pop_back();
                    localLastLevels[nID] = localTraversedLevels[nID].back();
                }
                for (VertexID nID: secondShared->nIDsToJoin) {
                    localTraversedLevels[nID].pop_back();
                    localLastLevels[nID] = localTraversedLevels[nID].back();
                }
            }

            // 释放内存
            for (VertexID nID = 0; nID < t.numNodes; ++nID) {
                localTrieLevels[nID].children = nullptr;
                localTrieLevels[nID].values = nullptr;
            }
        }
        for (int tid = 0; tid < maxThreads; ++tid) {
            delete[] threadVisited[tid];
            delete[] threadPartMatch[tid];
            freeSubtreeSharedJoinContext(threadContexts[tid], t, cs, secondShared);
        }
    }

    // 释放内存
    delete[] partMatch;
    for (int i = 0; i < prefixSum.back(); ++i) {
        delete[] nodeCandidates[i];
    }
    delete[] nodeCandidates;
    delete[] nodeCandCount;
    delete[] nodePoses;
    delete[] pCandidates[0];
    delete[] pCandidates;
    delete[] pCandCount;
    delete[] neighbors;
    delete[] neighborCount;
}