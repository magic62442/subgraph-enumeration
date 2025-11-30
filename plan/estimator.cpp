//
// Created by anonymous authors on 2024/2/27.
//

#include "estimator.h"

std::unordered_map<uint64_t, double> subsetToCard = std::unordered_map<uint64_t, double>();
double gTriCnt = 0.0;

void initPoses(const std::vector<VertexID> &order, const Graph &query, std::vector<std::vector<VertexID>> &vertexParents,
               std::vector<VertexID> &cartesianParent, const std::vector<std::vector<size_t>> &dist) {
    ui numAttributes = order.size();
    vertexParents = std::vector<std::vector<VertexID>>(numAttributes);
    cartesianParent = std::vector<VertexID>(numAttributes);
    bool cartesian = false;
    for (int i = 0; i < numAttributes; ++i) {
        VertexID u = order[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = order[j];
            if (query.getEdgeID(u, u2) != -1)
                vertexParents[i].push_back(u2);
        }
        if (i > 0 && vertexParents[i].empty()) cartesian = true;
    }
    if (cartesian) {
        for (ui i = 1; i < numAttributes; ++i) {
            VertexID u = order[i];
            if (vertexParents[i].empty()) {
                size_t minDist = std::numeric_limits<size_t>::max();
                for (VertexID j = 0; j < i; ++j) {
                    VertexID u2 = order[j];
                    if (dist[u][u2] < minDist) {
                        minDist = dist[u][u2];
                        cartesianParent[i] = u2;
                    }
                }
            }
        }
    }
}

// for a cartesian product at attribute u', find a path from a matched attribute u to u'
// use those that are reachable from v instead of all candidates of u'
void
optimizedCartesianProduct(CandidateSpace &cs, VertexID v, const std::vector<VertexID> &path, VertexID *candidate,
                          ui &candCount, VertexID maxCompared, VertexID minCompared, VertexID *partMatch) {
    VertexID u = path[0], uPrime = path.back();
    VertexID *result = new VertexID[cs.candSize(uPrime)];
    ui numResult = 0;
    candCount = 0;
    const VertexID **pathCand = new const VertexID *[path.size()];
    ui *pathCandCount = new ui[path.size()];
    ui *poses = new ui[path.size()];
    memset(pathCandCount, 0, sizeof(ui) * path.size());
    memset(poses, 0, sizeof(ui) * path.size());
    pathCandCount[0] = 1;
    int depth = 0;
    while (depth >= 0) {
        while (poses[depth] < pathCandCount[depth]) {
            VertexID u2 = path[depth];
            VertexID v2;
            if (depth == 0) v2 = v;
            else v2 = pathCand[depth][poses[depth]];
            ++poses[depth];
            VertexID next = path[depth + 1];
            ui nextCount;
            const VertexID *nextCand = cs.getNeighbors(u2, v2, next, nextCount);
            ui beginPos, endPos;
            if (partMatch != nullptr && maxCompared != 0) {
                beginPos = ComputeSetIntersection::BinarySearch(nextCand, 0, nextCount, maxCompared);
                nextCand += beginPos;
                nextCount -= beginPos;
            }
            if (partMatch != nullptr && minCompared != 0) {
                endPos = ComputeSetIntersection::BinarySearch(nextCand, 0, nextCount, minCompared);
                nextCount = endPos;
            }
            if (depth == path.size() - 2) {
                VertexID *temp = new VertexID[cs.candSize(uPrime)];
                VertexID i = 0, j = 0, k = 0;
                while (i < numResult && j < nextCount) {
                    if (result[i] < nextCand[j]) {
                        temp[k++] = result[i++];
                    } else if (result[i] > nextCand[j]) {
                        temp[k++] = nextCand[j++];
                    } else {
                        temp[k++] = result[i++];
                        ++j;
                    }
                }
                while (i < numResult) temp[k++] = result[i++];
                while (j < nextCount) temp[k++] = nextCand[j++];
                for (i = 0; i < k; ++i) result[i] = temp[i];
                numResult = k;
                delete[] temp;
            }
            else {
                ++depth;
                poses[depth] = 0;
                pathCand[depth] = nextCand;
                pathCandCount[depth] = nextCount;
            }
        }
        --depth;
    }

    for (int i = 0; i < numResult; ++i) {
        if (result[i] != v) {
            candidate[candCount] = result[i];
            ++candCount;
        }
    }
    delete[] result;
    delete[] pathCand[0];
    delete[] pathCand;
    delete[] pathCandCount;
    delete[] poses;
}

void
optimizedCartesianProduct(CandidateSpace &cs, VertexID v, const std::vector<VertexID> &path,
                          TrieNode *edgeColumn, ui &candCount, VertexID maxCompared, VertexID minCompared,
                          VertexID *partMatch) {
    VertexID u = path[0], uPrime = path.back();
    VertexID *result = new VertexID[cs.candSize(uPrime)];
    ui numResult = 0;
    candCount = 0;
    const VertexID **pathCand = new const VertexID *[path.size()];
    ui *pathCandCount = new ui[path.size()];
    ui *poses = new ui[path.size()];
    memset(pathCandCount, 0, sizeof(ui) * path.size());
    memset(poses, 0, sizeof(ui) * path.size());
    pathCandCount[0] = 1;
    int depth = 0;
    while (depth >= 0) {
        while (poses[depth] < pathCandCount[depth]) {
            VertexID u2 = path[depth];
            VertexID v2;
            if (depth == 0) v2 = v;
            else v2 = pathCand[depth][poses[depth]];
            ++poses[depth];
            VertexID next = path[depth + 1];
            ui nextCount;
            const VertexID *nextCand = cs.getNeighbors(u2, v2, next, nextCount);
            ui beginPos, endPos;
            if (partMatch != nullptr && maxCompared != 0) {
                beginPos = ComputeSetIntersection::BinarySearch(nextCand, 0, nextCount, maxCompared);
                nextCand += beginPos;
                nextCount -= beginPos;
            }
            if (partMatch != nullptr && minCompared != 0) {
                endPos = ComputeSetIntersection::BinarySearch(nextCand, 0, nextCount, minCompared);
                nextCount = endPos;
            }
            if (depth == path.size() - 2) {
                VertexID *temp = new VertexID[cs.candSize(uPrime)];
                VertexID i = 0, j = 0, k = 0;
                while (i < numResult && j < nextCount) {
                    if (result[i] < nextCand[j]) {
                        temp[k++] = result[i++];
                    } else if (result[i] > nextCand[j]) {
                        temp[k++] = nextCand[j++];
                    } else {
                        temp[k++] = result[i++];
                        ++j;
                    }
                }
                while (i < numResult) temp[k++] = result[i++];
                while (j < nextCount) temp[k++] = nextCand[j++];
                for (i = 0; i < k; ++i) result[i] = temp[i];
                numResult = k;
                delete[] temp;
            }
            else {
                ++depth;
                poses[depth] = 0;
                pathCand[depth] = nextCand;
                pathCandCount[depth] = nextCount;
            }
        }
        --depth;
    }

    for (int i = 0; i < numResult; ++i) {
        if (result[i] != v) {
            edgeColumn[candCount].value = result[i];
            ++candCount;
        }
    }
    delete[] result;
    delete[] pathCand[0];
    delete[] pathCand;
    delete[] pathCandCount;
    delete[] poses;
}

void
cardEstimateLabeled(const std::vector<VertexID> &order, const std::vector<std::vector<VertexID>> &vertexParents,
                    const std::vector<VertexID> &cartesianParent, CandidateSpace &cs, bool *visited,
                    VertexID *partMatch, VertexID **candidates, ui *candCount, std::vector<ui> &poses,
                    std::vector<double> &estimation) {
    if (order.empty()) {
        estimation.push_back(0.0);
        return;
    }
    if (order.size() == 1) {
        estimation.push_back(cs.candSize(order[0]));
        return;
    }
    if (order.size() == 2) {
        VertexID u1 = order[0], u2 = order[1];
        estimation.push_back(cs.candSize(u1));
        bool connected = true;
        double sz = 0.0;
        const VertexID *candSet = cs.candSet(u1);
        for (int i = 0; i < cs.candSize(u1); ++i) {
            VertexID v = candSet[i];
            ui numNeighbors;
            const VertexID *neighbors = cs.getNeighbors(u1, v, u2, numNeighbors);
            sz += numNeighbors;
        }
        if (sz == 0.0) sz = cs.candSize(u1) * cs.candSize(u2);
    }
    int depth = 0;
    VertexID u = order[0];
    VertexID **neighbors = new VertexID *[order.size()];
    ui maxSize = cs.getMaxSize();
    for (int i = 0; i < order.size(); ++i) {
        neighbors[i] = new VertexID[maxSize];
    }
    ui *neighborCount = new ui[order.size()];
    std::vector<ui> numSamples(order.size());
    std::vector<ui> numUsed(order.size(), 0);
    std::vector<std::vector<double>> estimate(order.size());
    for (int i = 0; i < order.size(); ++i) estimate[i] = std::vector<double>(order.size() - i, 0.0);
    numSamples[0] = SAMPLE_SIZE;
    std::vector<ui> totalCandCount(order.size(), 0);
    totalCandCount[0] = cs.candSize(u);
    // randomly sample a subset of candidates
    if (totalCandCount[0] < SAMPLE_SIZE) candCount[0] = totalCandCount[0];
    else candCount[0] = SAMPLE_SIZE;
    memcpy(candidates[0], cs.candSet(u), cs.candSize(u) * sizeof(VertexID));
    if (candCount[0] < cs.candSize(u))
        sampleKElements(candidates[0], cs.candSize(u), candCount[0]);
    ui totalCount = 0;
    numSamples[1] = SAMPLE_SIZE / candCount[0];
    while (depth >= 0) {
        while (poses[depth] < candCount[depth]) {
            VertexID v = candidates[depth][poses[depth]];
            ++poses[depth];
            visited[v] = true;
            partMatch[order[depth]] = v;
            ++depth;
            numUsed[depth] = 0;
            for (int i = 0; i < estimate[depth].size(); ++i) {
                estimate[depth][i] = 0.0;
            }
            poses[depth] = 0;
            u = order[depth];
            const std::vector<VertexID> &parents = vertexParents[depth];
            if (!parents.empty()) {
                for (int i = 0; i < parents.size(); ++i) {
                    VertexID pU = parents[i];
                    neighborCount[i] = 0;
                    ui numNeighbors;
                    const VertexID *vNeighbors = cs.getNeighbors(pU, partMatch[pU], u, numNeighbors);
                    for (int j = 0; j < numNeighbors; ++j) {
                        VertexID v2 = vNeighbors[j];
                        if (visited[v2]) continue;
                        neighbors[i][neighborCount[i]] = v2;
                        ++neighborCount[i];
                    }
                }
                ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), candidates[depth], totalCandCount[depth]);
            }
            else {
                const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent[depth], u);
                optimizedCartesianProduct(cs, partMatch[cartesianParent[depth]], path, candidates[depth],
                                          totalCandCount[depth], 0, 0, nullptr);
            }
            if (depth == order.size() - 1) {
                --depth;
                estimate[depth][0] += 1.0;
                estimate[depth][1] += static_cast<double>(totalCandCount[depth + 1]);
                totalCount += totalCandCount[depth + 1];
                ++numUsed[depth];
                visited[partMatch[order[depth]]] = false;
            }
            else {
                // randomly sample a subset of candidates. set the number of samples
                ui sampleSize;
                sampleSize = totalCandCount[depth] / SAMPLE_PORTION;
                sampleSize = std::min(numSamples[depth], sampleSize);
                if (totalCandCount[depth] != 0 && sampleSize == 0) sampleSize = 1;
                if (sampleSize < totalCandCount[depth])
                    sampleKElements(candidates[depth], totalCandCount[depth], sampleSize);
                candCount[depth] = sampleSize;
                if (sampleSize != 0) numSamples[depth + 1] = numSamples[depth] / sampleSize;
                else numSamples[depth + 1] = 0;
            }
        }
        --depth;
        if (depth >= 0) {
            // update estimations at the current depth
            estimate[depth][0] += 1.0;
            for (int i = 1; i < estimate[depth].size(); ++i) {
                estimate[depth][i] += static_cast<double>(estimate[depth + 1][i - 1]) * totalCandCount[depth] / candCount[depth];
            }
            numUsed[depth] += numUsed[depth + 1];
            if (poses[depth] < candCount[depth])
                numSamples[depth + 1] = (numSamples[depth] - numUsed[depth]) / (candCount[depth] - poses[depth]);
            visited[partMatch[order[depth]]] = false;
        }
    }

    estimation = estimate[0];
    for (int i = 0; i < order.size(); ++i) {
        delete[] neighbors[i];
    }
    delete[] neighbors;
    delete[] neighborCount;
}

void
generateCandidates(const VertexID **arrays, ui *counts, ui num, std::vector<VertexID> &tmpCand, bool *visited,
                   VertexID *cn, ui &cn_count, VertexID *partMatch, CandidateSpace &cs, VertexID current,
                   const std::vector<VertexID> &largerAttrs, const std::vector<VertexID> &smallerAttrs, VertexID cartesianParent) {
    VertexID maxCompared = 0, minCompared = (VertexID) - 1;
    for (VertexID u: largerAttrs) {
        if (partMatch[u] > maxCompared)
            maxCompared = partMatch[u];
    }
    for (VertexID u: smallerAttrs) {
        if (partMatch[u] < minCompared)
            minCompared = partMatch[u];
    }
    if (num != 0) {
        std::vector<ui> beginPoses(num), endPoses(num), oldCounts(num);
        if (!largerAttrs.empty()) {
            for (int i = 0; i < num; ++i) {
                beginPoses[i] = setBeginPos(arrays[i], counts[i], maxCompared);
                arrays[i] += beginPoses[i];
                counts[i] -= beginPoses[i];
            }
        }
        if (!smallerAttrs.empty()) {
            for (int i = 0; i < num; ++i) {
                endPoses[i] = setEndPos(arrays[i], counts[i], minCompared);
                oldCounts[i] = counts[i];
                counts[i] = endPoses[i];
            }
        }
        ComputeSetIntersection::MultiIntersection(arrays, counts, num, cn, cn_count);
    }
    else if (cartesianParent != 99) {
        const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, current);
        optimizedCartesianProduct(cs, partMatch[cartesianParent], path, cn,
                                  cn_count, maxCompared, minCompared, partMatch);
    }
    else {
        const VertexID *candSet = cs.candSet(current);
        cn_count = cs.candSize(current);
        memcpy(cn, candSet, sizeof(VertexID) * cn_count);
    }
    tmpCand.resize(0);
    for (int i = 0; i < cn_count; ++i) {
        if (!visited[cn[i]])
            tmpCand.push_back(cn[i]);
    }
    cn_count = tmpCand.size();
    memcpy(cn, tmpCand.data(), sizeof(VertexID) * cn_count);
}

void intersectionInfo(const std::vector<VertexID> &order, std::vector<std::vector<VertexID>> &vertexParents,
                      std::vector<VertexID> &cartesianParent, const Graph &query, const CandidateSpace &cs,
                      std::vector<std::vector<VertexID>> &symmetryRules, std::vector<int> &candidatesBefore,
                      std::vector<std::vector<VertexID>> &largerAttrs, std::vector<std::vector<VertexID>> &smallerAttrs,
                      bool skipReuse) {
    vertexParents = std::vector<std::vector<VertexID>>(order.size());
    cartesianParent = std::vector<VertexID>(order.size(), 99);
    largerAttrs = std::vector<std::vector<VertexID>>(order.size());
    smallerAttrs = std::vector<std::vector<VertexID>>(order.size());
    candidatesBefore = std::vector<int>(order.size(), 0);
    for (int i = 1; i < order.size(); ++i) {
        VertexID u = order[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = order[j];
            if (query.getEdgeID(u, u2) != -1)
                vertexParents[i].push_back(u2);
        }
        if (vertexParents[i].empty()) {
            size_t minDist = std::numeric_limits<size_t>::max();
            for (VertexID j = 0; j < i; ++j) {
                VertexID u2 = order[j];
                if (cs.dist[u][u2] < minDist) {
                    minDist = cs.dist[u][u2];
                    cartesianParent[i] = u2;
                }
            }
        }
    }
    VertexID maxVertex = 0;
    for (VertexID u: order) {
        if (u > maxVertex)
            maxVertex = u;
    }
    for (const auto &rule: symmetryRules) {
        VertexID u1 = rule[0];
        for (int i = 1; i < rule.size(); ++i) {
            VertexID u2 = rule[i];
            if (std::find(order.begin(), order.end(), u1) != order.end() &&
                std::find(order.begin(), order.end(), u2) != order.end()) {
                std::vector<int> attrToPos(maxVertex + 1);
                for (int j = 0; j < order.size(); ++j) {
                    attrToPos[order[j]] = j;
                }
                if (attrToPos[u1] < attrToPos[u2])
                    largerAttrs[attrToPos[u2]].push_back(u1);
                else
                    smallerAttrs[attrToPos[u1]].push_back(u2);
            }
        }
    }
    if (skipReuse) return;
    std::vector<std::vector<VertexID>> oldBefore = vertexParents;
    int maxCoverSize = 1, maxPos = -1;
    for (int j = 2; j < order.size(); ++j) {
        for (int k = 2; k < j; ++k) {
            if (orderedSubset(oldBefore[k], vertexParents[j]) &&
                unorderedSubset(largerAttrs[k], largerAttrs[j]) &&
                unorderedSubset(smallerAttrs[k], smallerAttrs[j])) {
                if (oldBefore[k].size() > maxCoverSize) {
                    maxCoverSize = oldBefore[k].size();
                    maxPos = k;
                }
            }
        }
        if (maxPos != -1) {
            std::unordered_set<VertexID> exists(oldBefore[maxPos].begin(), oldBefore[maxPos].end());
            vertexParents[j].clear();
            for (VertexID u: oldBefore[j]) {
                if (exists.find(u) == exists.end())
                    vertexParents[j].push_back(u);
            }
            candidatesBefore[j] = maxPos;
        }
    }
}

void intersectionInfoLastBag(const std::vector<VertexID> &order, const HyperTree &t, std::vector<std::vector<VertexID>> &vertexParents,
                             std::vector<VertexID> &cartesianParent, std::vector<std::vector<VertexID>> &nIDs, const Graph &query, const CandidateSpace &cs,
                             std::vector<std::vector<VertexID>> &symmetryRules, std::vector<int> &candidatesBefore,
                             std::vector<std::vector<VertexID>> &largerAttrs, std::vector<std::vector<VertexID>> &smallerAttrs) {
    const HyperNode &bag = t.nodes[t.numNodes - 1];
    vertexParents = std::vector<std::vector<VertexID>>(order.size());
    cartesianParent = std::vector<VertexID>(order.size(), 99);
    largerAttrs = std::vector<std::vector<VertexID>>(order.size());
    smallerAttrs = std::vector<std::vector<VertexID>>(order.size());
    bool cartesian = false;
    for (int i = 0; i < order.size(); ++i) {
        VertexID u = order[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = order[j];
            if (query.getEdgeID(u, u2) != -1)
                vertexParents[i].push_back(u2);
        }
        if (i > 0 && vertexParents[i].empty()) cartesian = true;
    }
    VertexID maxVertex = 0;
    for (VertexID u: order) {
        if (u > maxVertex)
            maxVertex = u;
    }
    for (const auto &rule: symmetryRules) {
        VertexID u1 = rule[0];
        for (int i = 1; i < rule.size(); ++i) {
            VertexID u2 = rule[i];
            if (std::find(order.begin(), order.end(), u1) != order.end() &&
                std::find(order.begin(), order.end(), u2) != order.end()) {
                std::vector<int> attrToPos(maxVertex + 1);
                for (int j = 0; j < order.size(); ++j) {
                    attrToPos[order[j]] = j;
                }
                if (attrToPos[u1] < attrToPos[u2])
                    largerAttrs[attrToPos[u2]].push_back(u1);
                else
                    smallerAttrs[attrToPos[u1]].push_back(u2);
            }
        }
    }
    nIDs = std::vector<std::vector<VertexID>>(order.size());
    // remove attributes that are covered by materialized bags
    for (int i = 0; i < order.size(); ++i) {
        std::set<VertexID> coveredAttrs;
        VertexID u = order[i];
        for (VertexID nID = 0; nID < t.numNodes - 1; ++nID) {
            if (std::find(t.nodes[nID].attributes, t.nodes[nID].attributes + t.nodes[nID].numAttributes, u) !=
                    t.nodes[nID].attributes + t.nodes[nID].numAttributes) {
                nIDs[i].push_back(nID);
                for (int j = 0; j < t.nodes[nID].numAttributes; ++j) {
                    VertexID u2 = t.nodes[nID].attributes[j];
                    coveredAttrs.insert(u2);
                }
            }
        }
        std::vector<ui> attrsToJoin;
        for (VertexID u2: vertexParents[i]) {
            if (coveredAttrs.find(u2) == coveredAttrs.end())
                attrsToJoin.push_back(u2);
        }
        vertexParents[i] = attrsToJoin;
    }
    if (cartesian) {
        for (ui i = 1; i < order.size(); ++i) {
            VertexID u = order[i];
            if (vertexParents[i].empty()) {
                size_t minDist = std::numeric_limits<size_t>::max();
                for (VertexID j = 0; j < i; ++j) {
                    VertexID u2 = order[j];
                    if (cs.dist[u][u2] < minDist) {
                        minDist = cs.dist[u][u2];
                        cartesianParent[i] = u2;
                    }
                }
            }
        }
    }
}

void
cardEstimateWrapper(const PatternGraph &p, const Graph &g, CandidateSpace &cs, std::vector<VertexID> &vertices,
                    bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount,
                    VertexID **totalCandidates, ui *totalCandCount, std::vector<ui> &poses,
                    std::vector<std::vector<VertexID>> &symmetryRules, std::vector<VertexID> &tmpCand) {
    std::vector<VertexID> order = RIOrder(p, vertices);
    if (order.empty()) {
        return;
    }
    CanonType id = getSubsetID(order);
    if (subsetToCard.find(id) != subsetToCard.end()) return;
    if (order.size() == 1) {
        subsetToCard[1 << order[0]] = cs.candSize(order[0]);
        return;
    }
    std::vector<std::vector<VertexID>> vertexParents, largerAttrs, smallerAttrs;
    std::vector<VertexID> cartesianParent;
    std::vector<int> candidatesBefore;
    intersectionInfo(order, vertexParents, cartesianParent, p, cs, symmetryRules, candidatesBefore,
                    largerAttrs, smallerAttrs);
    if (order.size() == 2) {
        if (largerAttrs[1].empty() && smallerAttrs[1].empty())
            subsetToCard[id] = (double)g.getNumEdges();
        else
            subsetToCard[id] = (double)g.getNumEdges() / 2;
        return;
    }
    cardEstimateUnlabeled(order, vertexParents, cartesianParent, candidatesBefore, largerAttrs, smallerAttrs,
                          cs, visited, partMatch, candidates, candCount, totalCandidates, totalCandCount, poses, tmpCand);
}

void
cardEstimateUnlabeled(const std::vector<VertexID> &order, const std::vector<std::vector<VertexID>> &vertexParents,
                      const std::vector<VertexID> &cartesianParent, std::vector<int> candidatesBefore,
                      std::vector<std::vector<VertexID>> &largerAttrs, std::vector<std::vector<VertexID>> &smallerAttrs,
                      CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount,
                      VertexID **totalCandidates, ui *totalCandCount, std::vector<ui> &poses, std::vector<VertexID> &tmpCand) {
    int depth = 0;
    VertexID u = order[0];
    const VertexID **neighbors = new const VertexID *[order.size()];
    ui maxSize = cs.getMaxSize();
    ui *neighborCount = new ui[order.size()];
    std::vector<ui> numSamples(order.size());
    std::vector<ui> numUsed(order.size(), 0);
    std::vector<std::vector<double>> estimate(order.size());
    for (int i = 0; i < order.size(); ++i) estimate[i] = std::vector<double>(order.size() - i, 0.0);
    numSamples[0] = SAMPLE_SIZE;
    totalCandCount[0] = cs.candSize(u);
    memcpy(candidates[0], cs.candSet(u), cs.candSize(u) * sizeof(VertexID));
    // randomly sample a subset of candidates
    if (totalCandCount[0] < SAMPLE_SIZE) candCount[0] = totalCandCount[0];
    else candCount[0] = SAMPLE_SIZE;
    memcpy(candidates[0], cs.candSet(u), cs.candSize(u) * sizeof(VertexID));
    if (candCount[0] < cs.candSize(u))
        sampleKElements(candidates[0], cs.candSize(u), candCount[0]);
    numSamples[1] = SAMPLE_SIZE / candCount[0];
    poses[0] = 0;
    while (depth >= 0) {
        while (poses[depth] < candCount[depth]) {
            VertexID v = candidates[depth][poses[depth]];
            ++poses[depth];
            visited[v] = true;
            partMatch[order[depth]] = v;
            ++depth;
            numUsed[depth] = 0;
            for (int i = 0; i < estimate[depth].size(); ++i) {
                estimate[depth][i] = 0.0;
            }
            poses[depth] = 0;
            u = order[depth];
            const std::vector<VertexID> &parents = vertexParents[depth];
            int reuse = 0;
            if (candidatesBefore[depth]) {
                reuse = 1;
                neighbors[0] = totalCandidates[candidatesBefore[depth]];
                neighborCount[0] = totalCandCount[candidatesBefore[depth]];
            }
            for (int i = 0; i < parents.size(); ++i) {
                VertexID pU = parents[i];
                neighbors[i + reuse] = cs.getNeighbors(pU, partMatch[pU], u, neighborCount[i + reuse]);
            }
            generateCandidates(neighbors, neighborCount, parents.size() + reuse, tmpCand, visited,
                               totalCandidates[depth], totalCandCount[depth], partMatch, cs, u, largerAttrs[depth],
                               smallerAttrs[depth], cartesianParent[depth]);
            if (depth == order.size() - 1) {
                --depth;
                estimate[depth][0] += 1.0;
                estimate[depth][1] += static_cast<double>(totalCandCount[depth + 1]);
                ++numUsed[depth];
                visited[partMatch[order[depth]]] = false;
            }
            else {
                // randomly sample a subset of candidates. set the number of samples
                ui sampleSize;
                sampleSize = totalCandCount[depth] / SAMPLE_PORTION;
                sampleSize = std::min(numSamples[depth], sampleSize);
                if (totalCandCount[depth] != 0 && sampleSize == 0) sampleSize = 1;
                if (sampleSize < totalCandCount[depth])
                    sampleKElements(totalCandidates[depth], totalCandCount[depth], candidates[depth], sampleSize);
                else memcpy(candidates[depth], totalCandidates[depth], totalCandCount[depth] * sizeof(VertexID));
                candCount[depth] = sampleSize;
                if (sampleSize != 0) numSamples[depth + 1] = numSamples[depth] / sampleSize;
                else numSamples[depth + 1] = 0;
            }
        }
        --depth;
        if (depth >= 0) {
            // update estimations at the current depth
            estimate[depth][0] += 1.0;
            for (int i = 1; i < estimate[depth].size(); ++i) {
                estimate[depth][i] += static_cast<double>(estimate[depth + 1][i - 1]) * totalCandCount[depth] / candCount[depth];
            }
            numUsed[depth] += numUsed[depth + 1];
            if (poses[depth] < candCount[depth])
                numSamples[depth + 1] = (numSamples[depth] - numUsed[depth]) / (candCount[depth] - poses[depth]);
            visited[partMatch[order[depth]]] = false;
        }
    }

    double estimation = 0.0;
    for (double est: estimate[0])
        if (est > estimation)
            estimation = est;
    if (estimation < 1.0) estimation = 1.0;
    CanonType id = getSubsetID(order);
    subsetToCard[id] = estimation;
    delete[] neighbors;
    delete[] neighborCount;
    if (gTriCnt == 0.0 && order.size() == 3 && vertexParents[2].size() == 2 && vertexParents[1].size() == 1) {
        gTriCnt = estimation;
    }
}

double estimateCartesian(VertexID u1, VertexID u2, CandidateSpace &cs) {
    double card = 0.0;
    const std::vector<VertexID> &path = cs.reconstructPath(u1, u2);
    VertexID *candidates = new VertexID [cs.candSize(u2)];
    ui candCount = 0;
    // sample vertices for cartesian product
    VertexID *vSamples = new VertexID [cs.candSize(u1)];
    memcpy(vSamples, cs.candSet(u1), cs.candSize(u1) * sizeof(VertexID));
    ui numSamples = cs.candSize(u1) / SAMPLE_PORTION;
    if (numSamples >= SAMPLE_SIZE) numSamples = SAMPLE_SIZE;
    if (numSamples == 0) numSamples = 1;
    sampleKElements(vSamples, cs.candSize(u1), numSamples);
    for (int i = 0; i < numSamples; ++i) {
        optimizedCartesianProduct(cs, vSamples[i], path, candidates, candCount, 0, 0, nullptr);
        card += double(candCount) * cs.candSize(u1) / numSamples;
    }
    delete[] candidates;
    delete[] vSamples;
    if (card < cs.candSize(u1)) card = static_cast<double>(cs.candSize(u1));
    if (card < cs.candSize(u2)) card = static_cast<double>(cs.candSize(u2));
    return card;
}

std::vector<VertexID> RIOrder(const Graph &query) {
    std::vector<VertexID> order;
    VertexID maxDegreeU = 0;
    ui maxDegree = 0;
    for (VertexID u = 1; u < query.getNumVertices(); ++u) {
        if (query.getDegree(u) > maxDegree) {
            maxDegree = query.getDegree(u);
            maxDegreeU = u;
        }
    }
    order.push_back(maxDegreeU);
    while (order.size() != query.getNumVertices()) {
        std::vector<int> numBackWard(query.getNumVertices(), -1);
        std::vector<int> breakTieNum1(query.getNumVertices(), 0);
        std::vector<int> breakTieNum2(query.getNumVertices(), 0);
        std::vector<VertexID> remaining;
        for (VertexID u = 0; u < query.getNumVertices(); ++u) {
            if (std::find(order.begin(), order.end(), u) != order.end())
                continue;
            else remaining.push_back(u);
            numBackWard[u] = 0;
            for (VertexID u2 : order) {
                if (query.getEdgeID(u, u2) != -1) ++numBackWard[u];
            }
        }
        for (VertexID u: remaining) {
            for (VertexID u2: order) {
                for (VertexID u3: remaining) {
                    if (query.getEdgeID(u, u3) != -1 && query.getEdgeID(u2, u3) != -1) {
                        ++breakTieNum1[u];
                        break;
                    }
                }
            }
        }
        for (VertexID u: remaining) {
            for (VertexID u2: remaining) {
                if (query.getEdgeID(u, u2) != -1) {
                    bool notConnectedToPrev = true;
                    for (VertexID u3: order) {
                        if (query.getEdgeID(u2, u3) != -1)
                            notConnectedToPrev = false;
                    }
                    if (notConnectedToPrev) ++breakTieNum2[u2];
                }
            }
        }
        int maxNumBack = -1, maxTie1 = 0, maxTie2 = 0;
        VertexID nextU = remaining[0];
        for (int i = 0; i < numBackWard.size(); ++i) {
            if (numBackWard[i] > maxNumBack) {
                maxNumBack = numBackWard[i];
                nextU = i;
            }
        }
        for (VertexID u: remaining) {
            if (numBackWard[u] == maxNumBack && (breakTieNum1[u] > maxTie1 || breakTieNum1[u] == maxTie1 &&
            breakTieNum2[u] >= maxTie2)) {
                nextU = u;
                maxTie1 = breakTieNum1[u];
                maxTie2 = breakTieNum2[u];
            }
        }
        order.push_back(nextU);
    }

    return order;
}

std::vector<VertexID> RIOrder(const Graph &query, const std::vector<VertexID> &vertices) {
    std::vector<VertexID> order;
    VertexID maxDegreeU = vertices[0];
    ui maxDegree = 0;
    for (VertexID u: vertices) {
        if (query.getDegree(u) > maxDegree) {
            maxDegree = query.getDegree(u);
            maxDegreeU = u;
        }
    }
    order.push_back(maxDegreeU);
    while (order.size() != vertices.size()) {
        std::vector<int> numBackWard(query.getNumVertices(), -1);
        std::vector<int> breakTieNum1(query.getNumVertices(), 0);
        std::vector<int> breakTieNum2(query.getNumVertices(), 0);
        std::vector<VertexID> remaining;
        for (VertexID u: vertices) {
            if (std::find(order.begin(), order.end(), u) != order.end())
                continue;
            else remaining.push_back(u);
            numBackWard[u] = 0;
            for (VertexID u2 : order) {
                if (query.getEdgeID(u, u2) != -1) ++numBackWard[u];
            }
        }
        for (VertexID u: remaining) {
            for (VertexID u2: order) {
                for (VertexID u3: remaining) {
                    if (query.getEdgeID(u, u3) != -1 && query.getEdgeID(u2, u3) != -1) {
                        ++breakTieNum1[u];
                        break;
                    }
                }
            }
        }
        for (VertexID u: remaining) {
            for (VertexID u2: remaining) {
                if (query.getEdgeID(u, u2) != -1) {
                    bool notConnectedToPrev = true;
                    for (VertexID u3: order) {
                        if (query.getEdgeID(u2, u3) != -1)
                            notConnectedToPrev = false;
                    }
                    if (notConnectedToPrev) ++breakTieNum2[u2];
                }
            }
        }
        int maxNumBack = -1, maxTie1 = 0, maxTie2 = 0;
        VertexID nextU = remaining[0];
        for (int i = 0; i < numBackWard.size(); ++i) {
            if (numBackWard[i] > maxNumBack) {
                maxNumBack = numBackWard[i];
                nextU = i;
            }
        }
        for (VertexID u: remaining) {
            if (numBackWard[u] == maxNumBack && (breakTieNum1[u] > maxTie1 || breakTieNum1[u] == maxTie1 &&
                                                                              breakTieNum2[u] > maxTie2)) {
                nextU = u;
                maxTie1 = breakTieNum1[u];
                maxTie2 = breakTieNum2[u];
            }
        }
        order.push_back(nextU);
    }

    return order;
}

std::vector<VertexID> RIOrder(const Graph &query, const std::vector<VertexID> &vertices, const std::vector<int> &repetitions) {
    std::vector<VertexID> order;
    VertexID maxDegreeU = vertices[0];
    ui maxDegree = 0;
    int maxRepetition = 0;

    // Find the vertex with the highest repetition count
    for (VertexID u : vertices) {
        if (repetitions[u] > maxRepetition) {
            maxRepetition = repetitions[u];
            maxDegreeU = u;
            maxDegree = query.getDegree(u);
        } else if (repetitions[u] == maxRepetition && query.getDegree(u) > maxDegree) {
            maxDegree = query.getDegree(u);
            maxDegreeU = u;
        }
    }

    order.push_back(maxDegreeU);

    while (order.size() != vertices.size()) {
        std::vector<int> numBackWard(query.getNumVertices(), -1);
        std::vector<int> breakTieNum1(query.getNumVertices(), 0);
        std::vector<int> breakTieNum2(query.getNumVertices(), 0);
        std::vector<VertexID> remaining;

        for (VertexID u : vertices) {
            if (std::find(order.begin(), order.end(), u) != order.end())
                continue;
            else remaining.push_back(u);
            numBackWard[u] = 0;
            for (VertexID u2 : order) {
                if (query.getEdgeID(u, u2) != -1) ++numBackWard[u];
            }
        }

        for (VertexID u : remaining) {
            for (VertexID u2 : order) {
                for (VertexID u3 : remaining) {
                    if (query.getEdgeID(u, u3) != -1 && query.getEdgeID(u2, u3) != -1) {
                        ++breakTieNum1[u];
                        break;
                    }
                }
            }
        }

        for (VertexID u : remaining) {
            for (VertexID u2 : remaining) {
                if (query.getEdgeID(u, u2) != -1) {
                    bool notConnectedToPrev = true;
                    for (VertexID u3 : order) {
                        if (query.getEdgeID(u2, u3) != -1)
                            notConnectedToPrev = false;
                    }
                    if (notConnectedToPrev) ++breakTieNum2[u2];
                }
            }
        }

        int maxNumBack = -1;
        int maxRepetitionNext = -1;
        int maxTie1 = 0;
        int maxTie2 = 0;
        VertexID nextU = remaining[0];

        for (VertexID u : remaining) {
//            if (numBackWard[u] > 0) {
//                if (repetitions[u] > maxRepetitionNext ||
//                    (repetitions[u] == maxRepetitionNext && numBackWard[u] > maxNumBack) ||
//                    (numBackWard[u] == maxNumBack && repetitions[u] == maxRepetitionNext && breakTieNum1[u] > maxTie1) ||
//                    (numBackWard[u] == maxNumBack && repetitions[u] == maxRepetitionNext && breakTieNum1[u] == maxTie1 && breakTieNum2[u] > maxTie2)) {
//                    nextU = u;
//                    maxNumBack = numBackWard[u];
//                    maxRepetitionNext = repetitions[u];
//                    maxTie1 = breakTieNum1[u];
//                    maxTie2 = breakTieNum2[u];
//                }
//            }
//            else{
                if (numBackWard[u] > maxNumBack ||
                    (numBackWard[u] == maxNumBack && repetitions[u] > maxRepetitionNext) ||
                    (numBackWard[u] == maxNumBack && repetitions[u] == maxRepetitionNext && breakTieNum1[u] > maxTie1) ||
                    (numBackWard[u] == maxNumBack && repetitions[u] == maxRepetitionNext && breakTieNum1[u] == maxTie1 && breakTieNum2[u] > maxTie2)) {
                    nextU = u;
                    maxNumBack = numBackWard[u];
                    maxRepetitionNext = repetitions[u];
                    maxTie1 = breakTieNum1[u];
                    maxTie2 = breakTieNum2[u];
                }
//            }
        }

        order.push_back(nextU);
    }

    return order;
}

std::vector<VertexID> RIOrder(const Graph &query, const std::vector<VertexID> &prefix, const std::vector<VertexID> &other) {
    if (prefix.empty()) return RIOrder(query, other);
    std::vector<VertexID> order = prefix;
    size_t total = prefix.size() + other.size();
    while (order.size() != total) {
        std::vector<int> numBackWard(query.getNumVertices(), -1);
        std::vector<int> breakTieNum1(query.getNumVertices(), 0);
        std::vector<int> breakTieNum2(query.getNumVertices(), 0);
        std::vector<VertexID> remaining;
        for (VertexID u: other) {
            if (std::find(order.begin(), order.end(), u) != order.end())
                continue;
            else remaining.push_back(u);
            numBackWard[u] = 0;
            for (VertexID u2 : order) {
                if (query.getEdgeID(u, u2) != -1) ++numBackWard[u];
            }
        }
        for (VertexID u: remaining) {
            for (VertexID u2: order) {
                for (VertexID u3: remaining) {
                    if (query.getEdgeID(u, u3) != -1 && query.getEdgeID(u2, u3) != -1) {
                        ++breakTieNum1[u];
                        break;
                    }
                }
            }
        }
        for (VertexID u: remaining) {
            for (VertexID u2: remaining) {
                if (query.getEdgeID(u, u2) != -1) {
                    bool notConnectedToPrev = true;
                    for (VertexID u3: order) {
                        if (query.getEdgeID(u2, u3) != -1)
                            notConnectedToPrev = false;
                    }
                    if (notConnectedToPrev) ++breakTieNum2[u2];
                }
            }
        }
        int maxNumBack = -1, maxTie1 = 0, maxTie2 = 0;
        VertexID nextU = remaining[0];
        for (int i = 0; i < numBackWard.size(); ++i) {
            if (numBackWard[i] > maxNumBack) {
                maxNumBack = numBackWard[i];
                nextU = i;
            }
        }
        for (VertexID u: remaining) {
            if (numBackWard[u] == maxNumBack && (breakTieNum1[u] > maxTie1 || breakTieNum1[u] == maxTie1 &&
                                                                              breakTieNum2[u] > maxTie2)) {
                nextU = u;
                maxTie1 = breakTieNum1[u];
                maxTie2 = breakTieNum2[u];
            }
        }
        order.push_back(nextU);
    }

    return order;
}

std::vector<VertexID> GQLOrder(const Graph &query, const CandidateSpace &cs) {
    std::vector<VertexID> order;
    ui minSize = cs.candSize(0);
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        if (cs.candSize(u) < minSize) minSize = cs.candSize(u);
    }
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        if (cs.candSize(u) == minSize) {
            order.push_back(u);
            break;
        }
    }
    while (order.size() != query.getNumVertices()) {
        std::vector<VertexID> candidates;
        for (VertexID u = 0; u < query.getNumVertices(); ++u) {
            if (std::find(order.begin(), order.end(), u) != order.end()) continue;
            bool connected = false;
            for (VertexID u2: order) {
                if (query.getEdgeID(u, u2) != -1) {
                    connected = true;
                    break;
                }
            }
            if (connected) candidates.push_back(u);
        }
        minSize = cs.candSize(candidates[0]);
        VertexID nextU = candidates[0];
        for (VertexID u: candidates) {
            if (cs.candSize(u) < minSize) {
                minSize = cs.candSize(u);
                nextU = u;
            }
        }
        order.push_back(nextU);
    }

    return order;
}

std::vector<VertexID> simpleOrder(const Graph &query, const CandidateSpace &cs, const std::vector<VertexID> &vertices,
                                  const std::vector<int> &repetitions) {
    ui n = query.getNumVertices();
    std::map<VertexID, VertexPriority> vset;
    std::vector<bool> visited(query.getNumVertices(), true);
    for (VertexID u: vertices) visited[u] = false;
    std::vector<std::vector<VertexID>> components;
    query.computeConnectedComponents(vertices, components);
    std::vector<size_t> componentWeight(components.size());
    for (int i = 0; i < components.size(); ++i) {
        componentWeight[i] = 1;
        for (int j = 0; j < components[i].size(); ++j) {
            componentWeight[i] *= cs.candSize(components[i][j]);
        }
    }
    for (int i = 0; i < components.size(); ++i) {
        for (int j = i + 1; j < components.size(); ++j) {
            if (componentWeight[j] > componentWeight[i]) {
                std::swap(components[i], components[j]);
            }
        }
    }
    std::vector<ui> prefixSum(components.size());
    prefixSum[0] = 0;
    for (int i = 1; i < prefixSum.size(); ++i) {
        prefixSum[i] = prefixSum[i - 1] + components[i - 1].size();
    }
    std::vector<VertexID> order;
    std::vector<ui> degrees(query.getNumVertices(), 0);
    for (int i = 0; i < vertices.size(); ++i) {
        VertexID u = vertices[i];
        for (int j = 0; j < vertices.size(); ++j) {
            VertexID u2 = vertices[j];
            if (query.getEdgeID(u, u2) != -1) {
                ++degrees[u];
                ++degrees[u2];
            }
        }
    }
    for (int i = 0; i < components.size(); ++i) {
        const std::vector<VertexID> &component = components[i];
        for (VertexID u: component) {
            VertexPriority vp;
            vp.num = repetitions[u];
            vp.backwardNbr = 0;
            vp.degree = degrees[u];
            vp.candSize = cs.candSize(u);
            vset[u] = vp;
        }
        VertexID next;
        VertexPriority min;

        while (!vset.empty()) {
            for (auto &it : vset) {
                VertexID u = it.first;
                VertexPriority priority = it.second;
                bool connected = false;
                if (order.size() == prefixSum[i]) connected = true;
                else {
                    for (VertexID u2: order) {
                        if (query.getEdgeID(u, u2) != -1) {
                            connected = true;
                            break;
                        }
                    }
                }
                if (!connected) continue;
                if (priority < min) {
                    min = priority;
                    next = u;
                }
            }
            order.push_back(next);
            vset.erase(next);
            visited[next] = true;
            ui numNbr;
            const VertexID *neighbors = query.getNeighbors(next, numNbr);
            for (int j = 0; j < numNbr; ++j) {
                VertexID u2 = neighbors[j];
                if (!visited[u2]) ++vset[u2].backwardNbr;
            }

            min.num = min.backwardNbr = min.degree = 0;
            min.candSize = std::numeric_limits<ui>::max();
        }
    }

    return order;
}

std::vector<VertexID>
simpleOrder(const Graph &query, const CandidateSpace &cs, const std::vector<VertexID> &prefix,
            std::vector<ui> &numBackWard, const std::vector<VertexID> &prevOrder, int noChangePos) {
    std::map<VertexID, VertexPriority> vset;
    std::vector<VertexID> vertices;
    for (int i = 0; i < noChangePos; ++i) {
        vertices.push_back(prevOrder[i]);
    }
    std::vector<VertexID> order;
    if (vertices.empty()) return order;
    std::vector<bool> visited(query.getNumVertices(), true);
    for (VertexID u: vertices) visited[u] = false;
    std::vector<std::vector<VertexID>> components;
    query.computeConnectedComponents(vertices, components);
    std::vector<size_t> componentWeight(components.size());
    for (int i = 0; i < components.size(); ++i) {
        componentWeight[i] = 1;
        for (int j = 0; j < components[i].size(); ++j) {
            componentWeight[i] *= cs.candSize(components[i][j]);
        }
    }
    for (int i = 0; i < components.size(); ++i) {
        for (int j = i + 1; j < components.size(); ++j) {
            if (componentWeight[j] > componentWeight[i]) {
                std::swap(components[i], components[j]);
            }
        }
    }
    std::vector<ui> prefixSum(components.size());
    prefixSum[0] = 0;
    for (int i = 1; i < prefixSum.size(); ++i) {
        prefixSum[i] = prefixSum[i - 1] + components[i - 1].size();
    }
    std::vector<ui> degrees(query.getNumVertices(), 0);
    for (int i = 0; i < vertices.size(); ++i) {
        VertexID u = vertices[i];
        for (int j = 0; j < vertices.size(); ++j) {
            VertexID u2 = vertices[j];
            if (query.getEdgeID(u, u2) != -1) {
                ++degrees[u];
                ++degrees[u2];
            }
        }
    }
    for (int i = 0; i < components.size(); ++i) {
        const std::vector<VertexID> &component = components[i];
        for (VertexID u: component) {
            VertexPriority vp;
            vp.num = 0;
            vp.backwardNbr = 0;
            for (VertexID u2: prefix) {
                if (query.getEdgeID(u, u2) != -1)
                    ++vp.backwardNbr;
            }
            vp.degree = degrees[u];
            vp.candSize = cs.candSize(u);
            vset[u] = vp;
        }
        VertexID next;
        VertexPriority min;

        while (!vset.empty()) {
            for (auto &it : vset) {
                VertexID u = it.first;
                VertexPriority priority = it.second;
                bool connected = false;
                if (order.size() == prefixSum[i]) connected = true;
                else {
                    for (VertexID u2: order) {
                        if (query.getEdgeID(u, u2) != -1) {
                            connected = true;
                            break;
                        }
                    }
                }
                if (!connected) continue;
                if (priority < min) {
                    min = priority;
                    next = u;
                }
            }
            order.push_back(next);
            numBackWard.push_back(vset[next].backwardNbr);
            vset.erase(next);
            visited[next] = true;
            ui numNbr;
            const VertexID *neighbors = query.getNeighbors(next, numNbr);
            for (int j = 0; j < numNbr; ++j) {
                VertexID u2 = neighbors[j];
                if (!visited[u2]) ++vset[u2].backwardNbr;
            }

            min.num = min.backwardNbr = min.degree = 0;
            min.candSize = std::numeric_limits<ui>::max();
        }
    }
    for (int i = noChangePos + 1; i < prevOrder.size(); ++i) {
        order.push_back(prevOrder[i]);
    }
    return order;
}

std::vector<ui>
computeNumBackWard(const Graph &query, const std::vector<VertexID> &prefix, const std::vector<VertexID> &localOrder) {
    std::vector<ui> numBackWard(localOrder.size(), 0);
    for (int i = 0; i < localOrder.size(); ++i) {
        VertexID u = localOrder[i];
        for (int j = 0; j < prefix.size(); ++j) {
            if (query.getEdgeID(prefix[j], u) != -1)
                ++numBackWard[i];
        }
        for (int j = 0; j < i; ++j) {
            if (query.getEdgeID(localOrder[j], u) != -1)
                ++numBackWard[i];
        }
    }
    return numBackWard;
}

std::vector<double> computeCost(const std::vector<VertexID> &order, const Graph &query, CandidateSpace &cs,
                                bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount) {
    std::vector<std::vector<VertexID>> vertexParents;
    std::vector<VertexID> cartesianParent;
    initPoses(order, query, vertexParents, cartesianParent, cs.dist);
    std::vector<double> cards;
    uint64_t id = 0;
    int maxPos = -1;
    for (int i = 0; i < order.size(); ++i) {
        id += 1 << order[i];
        if (subsetToCard.find(id) == subsetToCard.end()) maxPos = i;
    }
    std::vector<VertexID> subsequence(order.begin(), order.begin() + maxPos + 1);
    std::vector<ui> poses(query.getNumVertices(), 0);
    if (visited != nullptr)
        cardEstimateLabeled(subsequence, vertexParents, cartesianParent, cs, visited, partMatch, candidates, candCount, poses, cards);
    id = 0;
    for (int i = 0; i <= maxPos; ++i) {
        id += 1 << subsequence[i];
        subsetToCard[id] = cards[i];
    }
    for (int i = maxPos + 1; i < order.size(); ++i) {
        id += 1 << order[i];
        cards.push_back(subsetToCard[id]);
    }
    std::vector<double> costs(order.size());
    for (int i = 0; i < order.size(); ++i) {
        VertexID u = order[i];
        if (i == 0) {
            costs[i] = subsetToCard[1 << u];
            continue;
        }
        double listSize = 0;
        if (vertexParents[i].empty()) {
            VertexID u2 = cartesianParent[i];
            uint64_t id2 = (1 << u) + (1 << u2);
            if (subsetToCard.find(id2) == subsetToCard.end()) {
                subsetToCard[id2] = estimateCartesian(u2, u, cs);
            }
            listSize += subsetToCard[id2] / cs.candSize(u2);
        }
        else {
            for (ui u2: vertexParents[i]) {
                uint64_t edgeID = (1 << u) + (1 << u2);
                listSize += subsetToCard[edgeID] / subsetToCard[1 << u2];
            }
        }
        costs[i] += listSize * cards[i - 1];
    }

    return costs;
}

double computeCost(const std::vector<VertexID> &order, const std::vector<std::vector<VertexID>> &vertexParents,
                   const std::vector<VertexID> &cartesianParent, const std::vector<double> &cards, CandidateSpace &cs) {
    double cost = 0.0;
    for (int i = 0; i < order.size(); ++i) {
        VertexID u = order[i];
        if (i == 0) {
            cost += cs.candSize(u);
            continue;
        }
        double listSize = 0.0;
        if (vertexParents[i].empty()) {
            VertexID u2 = cartesianParent[i];
            uint64_t id2 = (1 << u) + (1 << u2);
            if (subsetToCard.find(id2) == subsetToCard.end()) {
                subsetToCard[id2] = estimateCartesian(u2, u, cs);
            }
            listSize = subsetToCard[id2] / cs.candSize(u2);
        }
        else {
            for (ui u2: vertexParents[i]) {
                uint64_t edgeID = (1 << u) + (1 << u2);
                listSize += subsetToCard[edgeID] / subsetToCard[1 << u2];
            }
        }
        cost += listSize * cards[i - 1];
    }

    return cost;
}

// order should extend prefix. The cost does not include the prefix cost
double computeCost(const std::vector<VertexID> &prefix, const std::vector<VertexID> &order, const Graph &query,
                   CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount) {
    double cost = 0.0;
    uint64_t prefixID = 0;
    int maxCardPos = 0;
    for (VertexID u : prefix) prefixID += 1 << u;
    uint64_t id = prefixID;
    for (int i = 0; i < order.size(); ++i) {
        id += 1 << order[i];
        if (subsetToCard.find(id) == subsetToCard.end())
            maxCardPos = i;
    }
    std::vector<std::vector<VertexID>> vertexParents;
    std::vector<VertexID> cartesianParent;
    initPoses(order, query, vertexParents, cartesianParent, cs.dist);
    std::vector<VertexID> subsequence(maxCardPos);
    for (int i = 0; i < maxCardPos; ++i) subsequence[i] = order[i];
    std::vector<double> estimation;
    std::vector<ui> poses(query.getNumVertices(), 0);
    if (!subsequence.empty())
        cardEstimateLabeled(subsequence, vertexParents, cartesianParent, cs, visited, partMatch, candidates, candCount, poses, estimation);
    id = prefixID;
    for (int i = 0; i < subsequence.size(); ++i) {
        id += 1 << order[i];
        if (subsetToCard.find(id) == subsetToCard.end()) subsetToCard[id] = estimation[i];
    }
    id = prefixID;
    std::vector<VertexID> nodeOrder = prefix;
    for (VertexID u: order) nodeOrder.push_back(u);
    initPoses(nodeOrder, query, vertexParents, cartesianParent, cs.dist);
    for (int i = 0; i < order.size(); ++i) {
        double card = 1.0;
        if (id != 0) card = subsetToCard[id];
        VertexID u = order[i];
        if (prefix.empty() && i == 0) {
            cost += cs.candSize(u);
            id += 1 << order[i];
            continue;
        }
        double listSize = std::numeric_limits<double>::max();
        if (vertexParents[i + prefix.size()].empty()) {
            VertexID u2 = cartesianParent[i + prefix.size()];
            uint64_t id2 = (1 << u) + (1 << u2);
            if (subsetToCard.find(id2) == subsetToCard.end()) {
                subsetToCard[id2] = estimateCartesian(u2, u, cs);
            }
            listSize = subsetToCard[id2] / cs.candSize(u2);
        }
        else {
            for (VertexID u2: vertexParents[i + prefix.size()]) {
                uint64_t edgeID = (1 << u) + (1 << u2);
                listSize += subsetToCard[edgeID] / subsetToCard[1 << u2];
            }
        }
        cost += listSize * card;
        id += 1 << order[i];
    }

    return cost;
}

double computeCost(const PrefixNode *pt, const std::vector<std::vector<VertexID>> &localOrders, const Graph &query,
                   CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount,
                   const std::vector<std::vector<VertexID>> &matchedAttrs, const std::vector<VertexID> &prevAttrs,
                   const std::vector<double> &factors) {
    // pt is the virtual root
    std::vector<VertexID> attrsInPath;
    int depth = 0;
    const PrefixNode *pn = pt;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes(height, nullptr);
    // generate root-to-leaf paths
    std::vector<std::vector<VertexID>> paths;
    std::vector<bool> computed(query.getNumVertices(), false);
    double cost = 0.0;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            attrsInPath.push_back(u);
            nodes[depth] = current;
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else {
                std::vector<VertexID> path = prevAttrs;
                for (VertexID u2 : attrsInPath) path.push_back(u2);
                paths.push_back(path);
                if (childPoses[depth] < pn -> children.size())
                    attrsInPath.pop_back();
            }
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = pt;
            else pn = nodes[depth - 1];
            attrsInPath.pop_back();
            if (childPoses[depth] < pn -> children.size()) attrsInPath.pop_back();
        }
    }
    for (const auto &path : paths) {
        std::vector<double> costs = computeCost(path, query, cs, visited, partMatch, candidates, candCount);
        for (int i = 0; i < costs.size(); ++i) {
            VertexID u = path[i];
            if (!computed[u]) {
                computed[u] = true;
                cost += costs[i];
            }
        }
    }
    pn = pt;
    depth = 0;
    childPoses = std::vector<ui>(height, 0);
    attrsInPath.clear();
    for (VertexID nID : pt -> nIDsToCall)
        cost += computeCost(matchedAttrs[nID], localOrders[nID], query, cs, visited, partMatch, candidates, candCount) * factors[nID];
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            attrsInPath.push_back(u);
            nodes[depth] = current;
            for (VertexID nID: current -> nIDsToCall) {
                std::vector<VertexID> prefix = matchedAttrs[nID];
                for (VertexID u2 : attrsInPath) prefix.push_back(u2);
                cost += computeCost(prefix, localOrders[nID], query, cs, visited, partMatch, candidates, candCount) * factors[nID];
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else if (childPoses[depth] < pn -> children.size())
                attrsInPath.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = pt;
            else pn = nodes[depth - 1];
            attrsInPath.pop_back();
            if (childPoses[depth] < pn -> children.size()) attrsInPath.pop_back();
        }
    }

    return cost;
}

ui maxNumBackWard(const std::vector<VertexID> &order, const Graph &query) {
    ui maxNum = 1;
    for (int i = 1; i < order.size(); ++i) {
        ui num = 0;
        VertexID u = order[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = order[j];
            if (query.getEdgeID(u, u2) != -1) ++num;
        }

        if (maxNum < num) maxNum = num;
    }

    return maxNum;
}

bool saveSubsetToCard(std::ofstream& ofs) {
    if (!ofs) return false;
    uint64_t map_size = subsetToCard.size();
    ofs.write(reinterpret_cast<const char*>(&map_size), sizeof(map_size));
    for (const auto& pair : subsetToCard) {
        ofs.write(reinterpret_cast<const char*>(&pair.first), sizeof(pair.first));
        double scaled = static_cast<double>(static_cast<int64_t>(pair.second * 100.0)) / 100.0;
        ofs.write(reinterpret_cast<const char*>(&scaled), sizeof(scaled));
    }

    return ofs.good();
}

bool loadSubsetToCard(std::ifstream& ifs) {
    if (!ifs) return false;
    subsetToCard.clear();
    subsetToCard[0] = 0;
    uint64_t map_size;
    ifs.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
    for (uint64_t i = 0; i < map_size; ++i) {
        uint64_t key;
        double value;
        ifs.read(reinterpret_cast<char*>(&key), sizeof(key));
        ifs.read(reinterpret_cast<char*>(&value), sizeof(value));
        subsetToCard[key] = value;
    }

    return ifs.good();
}