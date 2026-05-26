//
// Created by anonymous authors on 2024/5/22.
//

#include "optimizer.h"

void setTDExtention(HyperTree &t, const Graph &query) {
#ifdef ALL_LEVEL
    t.extendLevel = t.defaultPartition.size();
//    const HyperNode &last = t.nodes[t.numNodes - 1];
//    for (int i = last.prefixSize; i < last.numAttributes; ++i) {
//        ++t.extendLevel;
//        t.defaultPartition.push_back(last.attributes[i]);
//    }
#else
    if (!t.newGlobalNode) {
        VertexID nID = t.numNodes - 1;
        for (int i = t.nodes[nID].prefixSize; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (t.v2n[u].size() > 1)
                t.extendLevel = i + 1;
        }
        for (int i = 0; i < t.extendLevel; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (std::find(t.defaultPartition.begin(), t.defaultPartition.end(), u) == t.defaultPartition.end())
                t.defaultPartition.push_back(u);
        }
    }
    else {
        t.extendLevel = t.nodes[t.numNodes - 1].numAttributes;
        for (int i = 0; i < t.extendLevel; ++i) {
            VertexID u = t.nodes[t.numNodes - 1].attributes[i];
            if (std::find(t.defaultPartition.begin(), t.defaultPartition.end(), u) == t.defaultPartition.end())
                t.defaultPartition.push_back(u);
        }
    }
#endif
    t.buildTraverseStruct(query);
}

std::vector<VertexID> globalCandidates(const Graph &query, HyperTree &t) {
    std::vector<VertexID> globalCand;
    t.numAttributes = query.getNumVertices();
    delete[] t.v2n;
    t.v2n = new std::vector<VertexID> [query.getNumVertices()];
    for (ui i = 0; i < t.numNodes; ++i) {
        for (ui j = 0; j < t.nodes[i].numAttributes; ++j) {
            t.v2n[t.nodes[i].attributes[j]].push_back(i);
        }
    }
    std::vector<std::vector<VertexID>> sharedAttrs(t.numNodes);
    std::vector<VertexID> globalAttr;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (t.v2n[u].size() > 1) {
                if (std::find(globalAttr.begin(), globalAttr.end(), u) == globalAttr.end())
                    globalAttr.push_back(u);
                sharedAttrs[nID].push_back(u);
            }
        }
        std::sort(sharedAttrs[nID].begin(), sharedAttrs[nID].end());
    }
    std::sort(globalAttr.begin(), globalAttr.end());
    bool globalBag = true;
    for (int nID = t.numNodes - 1; nID >= 0; --nID) {
        if (std::includes(sharedAttrs[nID].begin(), sharedAttrs[nID].end(), globalAttr.begin(), globalAttr.end())) {
            globalBag = false;
            bool exists = false;
            for (VertexID nID2: globalCand) {
                if (t.nodes[nID].canonValue == t.nodes[nID2].canonValue && sharedAttrs[nID] == sharedAttrs[nID2])
                    exists = true;
            }
            if (!exists) globalCand.push_back(nID);
        }
    }
    return globalCand;
}

void optCostPlan(const PatternGraph &p, const Graph &g, CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount,
                 VertexID **totalCandidates, ui *totalCandCount, std::vector<ui> &poses, std::vector<VertexID> &tmpCand,
                 HyperTree &t, PrefixNode *&pt) {
    std::vector<std::vector<size_t>> optimalOrders;
    double minWidth;
    std::vector<HyperTree> trees = buildOptimalFHDsByEnumeration(p, optimalOrders, minWidth);
    double minCost = std::numeric_limits<double>::max();
    PrefixNode *bestPT;
    std::vector<std::vector<VertexID>> bestOrders;
    std::vector<std::vector<VertexID>> previousSymmetry;
    for (int i = 0 ; i < trees.size(); ++i) {
        HyperTree &tree = trees[i];
        std::vector<VertexID> globalCands = globalCandidates(p, tree);
        if (globalCands.empty()) {
            std::vector<VertexID> globalAttr;
            for (VertexID u = 0; u < p.getNumVertices(); ++u) {
                if (tree.v2n[u].size() > 1) globalAttr.push_back(u);
            }
            tree.addGlobalNode(globalAttr);
            tree.selectSymmetry(p);
            if (tree.symmetryRules != previousSymmetry)
                subsetToCard.clear();
            optCostOrder(p, g, tree, cs, visited, partMatch, candidates, candCount, totalCandidates, totalCandCount, poses, tmpCand, tree.symmetryRules, bestPT, bestOrders, minCost);
        }
        else {
            tree.selectSymmetry(p);
            if (tree.symmetryRules != previousSymmetry)
                subsetToCard.clear();
            for (VertexID nID: globalCands) {
                swapBags(tree, nID, tree.numNodes - 1);
                optCostOrder(p, g, tree, cs, visited, partMatch, candidates, candCount, totalCandidates, totalCandCount, poses, tmpCand, tree.symmetryRules, bestPT, bestOrders, minCost);
                swapBags(tree, nID, tree.numNodes - 1);
            }
        }
        previousSymmetry = tree.symmetryRules;
    }
    pt = bestPT->clone();
    std::vector<std::vector<VertexID>> sharedAttrs(bestOrders.size());
    for (VertexID nID = 0; nID < bestOrders.size(); ++nID) {
        for (int i = 0; i < bestOrders[nID].size(); ++i) {
            VertexID u = bestOrders[nID][i];
            bool exists = false;
            if (std::find(bestOrders.back().begin(), bestOrders.back().end(), u) != bestOrders.back().end()) {
                exists = true;
            }
            if (exists) sharedAttrs[nID].push_back(u);;
        }
        std::sort(sharedAttrs[nID].begin(), sharedAttrs[nID].end());
    }
//    bestOrders[1][2] = 1;
//    bestOrders[1][3] = 3;
//    bestOrders[2][2] = 0;
//    bestOrders[2][3] = 3;
    t.numAttributes = p.getNumVertices();
    t.numNodes = bestOrders.size();
    t.nodes = new HyperNode[t.numNodes];
    t.v2n = new std::vector<VertexID>[t.numAttributes];
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (VertexID u: bestOrders[nID])
            t.v2n[u].push_back(nID);
    }
    buildFromPrefixTree(pt, bestOrders, t, sharedAttrs, p, cs);
    t.largerAttrs = std::vector<std::vector<VertexID>>(t.numAttributes);
    t.smallerAttrs = std::vector<std::vector<VertexID>>(t.numAttributes);
    t.hasLargerAttrs = std::vector<bool>(t.numAttributes, false);
    t.hasSmallerAttrs = std::vector<bool>(t.numAttributes, false);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        t.nodes[nID].largerAttrs = std::vector<std::vector<VertexID>>(t.nodes[nID].numAttributes);
        t.nodes[nID].smallerAttrs = std::vector<std::vector<VertexID>>(t.nodes[nID].numAttributes);
        t.nodes[nID].hasLargerAttrs = std::vector<bool>(t.nodes[nID].numAttributes, false);
        t.nodes[nID].hasSmallerAttrs = std::vector<bool>(t.nodes[nID].numAttributes, false);
    }
    pt->refineNIDsToBuild(t);
#ifdef ALL_LEVEL
    t.defaultPartition = std::vector<VertexID>(t.nodes[t.numNodes - 1].attributes, t.nodes[t.numNodes - 1].attributes + t.nodes[t.numNodes - 1].numAttributes);
#endif
    setTDExtention(t, p);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            t.nodes[nID].id += 1 << u;
        }
    }
}

double computeOrderCost(const PatternGraph &p, const Graph &g, HyperTree &t, VertexID nID, CandidateSpace &cs, bool *visited,
                      VertexID *partMatch, VertexID **candidates, ui *candCount, const std::vector<VertexID> &order,
                      VertexID **totalCandidates, ui *totalCandCount, std::vector<ui> &poses, std::vector<VertexID> &tmpCand,
                      std::vector<std::vector<VertexID>> &symmetryRules, std::vector<VertexID> &visitedBag,
                      std::vector<std::vector<std::vector<VertexID>>> &attributesBefore,
                      std::vector<std::vector<std::vector<VertexID>>> &smallerAttrs,
                      std::vector<std::vector<std::vector<VertexID>>> &largerAttrs,
                      std::vector<std::vector<int>> &candidatesBefore, std::vector<std::vector<VertexID>> &cartesianParent,
                      std::vector<PrefixNode *> &path, std::set<PrefixNode *> &computed) {
    double cost = 0.0;
    std::vector<std::vector<VertexID>> nIDs;
    if (nID != t.numNodes - 1)
        intersectionInfo(order, attributesBefore[nID], cartesianParent[nID], p, cs, symmetryRules, candidatesBefore[nID],
                     largerAttrs[nID], smallerAttrs[nID], true);
    else
        intersectionInfoLastBag(order, t, attributesBefore[nID], cartesianParent[nID], nIDs, p, cs, symmetryRules,
                                candidatesBefore[nID], largerAttrs[nID], smallerAttrs[nID]);
    uint64_t id = 1 << order[0];
    std::vector<VertexID> vertices;
    vertices.push_back(order[0]);
    for (int i = 1; i < order.size(); ++i) {
        VertexID u = order[i];
        if (path[i] && computed.find(path[i]) != computed.end()) {
            id += 1 << u;
            vertices.push_back(u);
            continue;
        }
        if (subsetToCard.find(id) == subsetToCard.end() && id != 0) {
            cardEstimateWrapper(p, g, cs, vertices, visited, partMatch, candidates, candCount, totalCandidates,
                                totalCandCount, poses, symmetryRules, tmpCand);
        }
        double card = subsetToCard[id];
        id += 1 << u;
        vertices.push_back(u);
        double listSize = 0.0;
        if (i >= 2) {
            int maxCoverSize = 1, maxPos, maxPrevious = -1;
            if (!(nID == t.numNodes - 1 && !nIDs[i].empty())) {
                // previous bags
                bool allInPath = true;
                for (VertexID u2: attributesBefore[nID][i]) {
                    bool inPath = false;
                    for (PrefixNode *node: path) {
                        if (!node) break;
                        if (node->u == u2) inPath = true;
                    }
                    if (!inPath) {
                        allInPath = false;
                        break;
                    }
                }
                if (allInPath) {
                    for (VertexID nID2: visitedBag) {
                        for (int j = 2; j < attributesBefore[nID2].size(); ++j) {
                            if (path[j] != nullptr && !path[j]->nIDsToJoin.empty()) continue;
                            if (orderedSubset(attributesBefore[nID2][j], attributesBefore[nID][i]) &&
                                unorderedSubset(largerAttrs[nID2][j], largerAttrs[nID][i]) &&
                                unorderedSubset(smallerAttrs[nID2][j], smallerAttrs[nID][i])) {
                                if (attributesBefore[nID2][j].size() > maxCoverSize) {
                                    maxCoverSize = attributesBefore[nID2][j].size();
                                    maxPrevious = nID2;
                                    maxPos = j;
                                }
                            }
                        }
                    }
                }
                // current bag
                for (int j = 2; j < i; ++j) {
                    if (path[j] != nullptr && !path[j]->nIDsToJoin.empty()) continue;
                    if (orderedSubset(attributesBefore[nID][j], attributesBefore[nID][i]) &&
                        unorderedSubset(largerAttrs[nID][j], largerAttrs[nID][i]) &&
                        unorderedSubset(smallerAttrs[nID][j], smallerAttrs[nID][i])) {
                        if (attributesBefore[nID][j].size() > maxCoverSize) {
                            maxCoverSize = attributesBefore[nID][j].size();
                            maxPrevious = nID;
                            maxPos = j;
                        }
                    }
                }
                std::vector<VertexID> vertexParents;
                std::set<VertexID> exists;
                if (maxCoverSize > 1) {
                    for (VertexID u2: attributesBefore[maxPrevious][maxPos]) exists.insert(u2);
                }
                for (VertexID u2: attributesBefore[nID][i]) {
                    if (exists.find(u2) == exists.end())
                        vertexParents.push_back(u2);
                }
                for (VertexID u2: vertexParents) {
                    if (std::find(largerAttrs[nID][i].begin(), largerAttrs[nID][i].end(), u2) == largerAttrs[nID][i].end() &&
                        std::find(smallerAttrs[nID][i].begin(), smallerAttrs[nID][i].end(), u2) == smallerAttrs[nID][i].end()) {
                        listSize += (double)g.getNumEdges() / g.getNumVertices();
                    }
                    else listSize += (double)g.getNumEdges() / g.getNumVertices() / 2;
                }
            }
            else {
                listSize += nIDs[i].size() * gTriCnt / g.getNumEdges();
                for (VertexID u2: attributesBefore[nID][i]) {
                    if (std::find(largerAttrs[nID][i].begin(), largerAttrs[nID][i].end(), u2) == largerAttrs[nID][i].end() &&
                        std::find(smallerAttrs[nID][i].begin(), smallerAttrs[nID][i].end(), u2) == smallerAttrs[nID][i].end()) {
                        listSize += (double)g.getNumEdges() / g.getNumVertices();
                    }
                    else listSize += (double)g.getNumEdges() / g.getNumVertices() / 2;
                }
            }
        }
        else {
            if (std::find(largerAttrs[nID][i].begin(), largerAttrs[nID][i].end(), order[0]) == largerAttrs[nID][i].end() &&
                std::find(smallerAttrs[nID][i].begin(), smallerAttrs[nID][i].end(), order[1]) == smallerAttrs[nID][i].end()) {
                listSize += (double)g.getNumEdges() / g.getNumVertices();
            }
            else listSize += (double)g.getNumEdges() / g.getNumVertices() / 2;
        }
        cost += card * listSize;
        if (path[i]) computed.insert(path[i]);
    }
    return cost;
}

// type: 0: optimal, 1: only materialization, 2: only intersection
void optCostOrder(const PatternGraph &p, const Graph &g, HyperTree &t, CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount,
                  VertexID **totalCandidates, ui *totalCandCount, std::vector<ui> &poses, std::vector<VertexID> &tmpCand,
                  std::vector<std::vector<VertexID>> &symmetryRules, PrefixNode *&bestPT, std::vector<std::vector<VertexID>> &bestOrders, double &minCost,
                  int type) {
    if (t.numNodes == 1) {
        bestPT = new PrefixNode(99);
        bestPT->nIDsToCall = {0};
        bestOrders.emplace_back();
        for (VertexID u = 0; u < t.numAttributes; ++u)
            bestOrders[0].push_back(u);
        return;
    }
    std::vector<PrefixNode *> prefixTrees;
    genAllPrefixTree(t, p, cs, prefixTrees);
    std::vector<PrefixNode *> nodes(p.getNumVertices(), nullptr);
    std::vector<ui> childPoses(p.getNumVertices(), 0);
    std::vector<std::vector<VertexID>> orders(t.numNodes);
    for (int i = 0; i < prefixTrees.size(); ++i) {
        std::set<PrefixNode *> computed;
        std::vector<VertexID> visitedBag;
        std::vector<std::vector<std::vector<VertexID>>> attributesBefore(t.numNodes), smallerAttrs(t.numNodes), largerAttrs(t.numNodes);
        std::vector<std::vector<int>> candidatesBefore(t.numNodes);
        std::vector<std::vector<VertexID>> cartesianParent(t.numNodes);
        std::vector<VertexID> attrsInPath;
        PrefixNode *pn = prefixTrees[i];
        std::vector<ui> trieNumLevel(t.numNodes, 0);
        int depth = 0;
        childPoses[0] = 0;
        uint64_t id = 0;
        double cost = 0.0;
        for (VertexID nID: pn->nIDsToCall) {
            if (nID == t.numNodes - 1) continue;
            std::vector<VertexID> vertices(t.nodes[nID].attributes, t.nodes[nID].attributes + t.nodes[nID].numAttributes);
            orders[nID] = RIOrder(p, vertices);
            // local join cost
            if (type != 1)
                cost += computeOrderCost(p, g, t, nID, cs, visited, partMatch, candidates, candCount, orders[nID], totalCandidates,
                                         totalCandCount, poses, tmpCand, symmetryRules, visitedBag, attributesBefore,
                                         smallerAttrs, largerAttrs, candidatesBefore, cartesianParent, nodes, computed);
            // materialization cost
            if (subsetToCard.find(t.nodes[nID].id) == subsetToCard.end()) {
                cardEstimateWrapper(p, g, cs, vertices, visited, partMatch, candidates, candCount, totalCandidates,
                                    totalCandCount, poses, symmetryRules, tmpCand);
            }
            double card = subsetToCard.at(t.nodes[nID].id);
            if (type != 2) cost += card * t.nodes[nID].numAttributes * COEEFICIENT;
        }
        while (depth >= 0) {
            while (childPoses[depth] < pn -> children.size()) {
                PrefixNode *current = pn -> children[childPoses[depth]];
                ++childPoses[depth];
                VertexID u = current -> u;
                nodes[depth] = current;
                id += 1 << u;
                attrsInPath.push_back(u);
                for (VertexID nID: current -> nIDsToCall) {
                    std::vector<VertexID> remaining;
                    for (int j = 0; j < t.nodes[nID].numAttributes; ++j) {
                        if (std::find(attrsInPath.begin(), attrsInPath.end(), t.nodes[nID].attributes[j]) == attrsInPath.end())
                            remaining.push_back(t.nodes[nID].attributes[j]);
                    }
                    if (nID != t.numNodes - 1) {
                        trieNumLevel[nID] = t.nodes[nID].numAttributes;
                        for (int j = 0; j <= depth; ++j) {
                            if (nodes[j]->pathToGlobal) {
                                --trieNumLevel[nID];
                            }
                            else break;
                        }
                        orders[nID] = RIOrder(p, attrsInPath, remaining);
                        // local join cost
                        if (type != 1)
                            cost += computeOrderCost(p, g, t, nID, cs, visited, partMatch, candidates, candCount, orders[nID], totalCandidates,
                                                 totalCandCount, poses, tmpCand, symmetryRules, visitedBag, attributesBefore,
                                                 smallerAttrs, largerAttrs, candidatesBefore, cartesianParent, nodes, computed);
                        // materialization cost
                        if (trieNumLevel[nID] > 1) {
                            if (subsetToCard.find(t.nodes[nID].id) == subsetToCard.end()) {
                                std::vector<VertexID> vertices(t.nodes[nID].attributes, t.nodes[nID].attributes + t.nodes[nID].numAttributes);
                                cardEstimateWrapper(p, g, cs, vertices, visited, partMatch, candidates, candCount, totalCandidates,
                                                    totalCandCount, poses, symmetryRules, tmpCand);
                            }
                            double card = subsetToCard.at(t.nodes[nID].id);
                            if (type != 2)
                                cost += card * trieNumLevel[nID] * COEEFICIENT;
                        }
                    }
                    else {
                        orders[nID] = RIOrder(p, attrsInPath, remaining);
                        // global join cost
                        if (type != 1)
                            cost += computeOrderCost(p, g, t, nID, cs, visited, partMatch, candidates, candCount, orders[nID], totalCandidates,
                                                 totalCandCount, poses, tmpCand, symmetryRules, visitedBag, attributesBefore,
                                                 smallerAttrs, largerAttrs, candidatesBefore, cartesianParent, nodes, computed);
                    }
                    visitedBag.push_back(nID);
                }
                if (!current -> children.empty()) {
                    ++depth;
                    childPoses[depth] = 0;
                }
                else if (childPoses[depth] < pn -> children.size()) {
                    id -= 1 << u;
                    attrsInPath.pop_back();
                }
                if (depth > 0) pn = nodes[depth - 1];
            }
            --depth;
            if (depth >= 0) {
                if (depth == 0) pn = prefixTrees[i];
                else pn = nodes[depth - 1];
                id -= 1 << attrsInPath.back();
                attrsInPath.pop_back();
                if (childPoses[depth] < pn -> children.size()) {
                    id -= 1 << attrsInPath.back();
                    attrsInPath.pop_back();
                }
            }
        }
        for (VertexID nID: pn->nIDsToCall) {
            if (nID != t.numNodes - 1) continue;
            std::vector<VertexID> vertices(t.nodes[nID].attributes, t.nodes[nID].attributes + t.nodes[nID].numAttributes);
            orders[nID] = RIOrder(p, vertices);
            // global join cost
            if (type != 1)
                cost += computeOrderCost(p, g, t, nID, cs, visited, partMatch, candidates, candCount, orders[nID], totalCandidates,
                                     totalCandCount, poses, tmpCand, symmetryRules, visitedBag, attributesBefore,
                                     smallerAttrs, largerAttrs, candidatesBefore, cartesianParent, nodes, computed);
        }
        if (cost < minCost) {
            minCost = cost;
            bestPT = prefixTrees[i];
            bestOrders = orders;
        }
        
    }
}
