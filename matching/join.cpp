//
// Created by anonymous authors on 2024/3/8.
//

#include "join.h"

size_t gNumResult = 0;

void twoIntersection(const VertexID *lArray, ui lBegin, ui lEnd, const VertexID *rArray, ui rBegin, ui rEnd, ui **iters, VertexID *cn, ui &cn_count) {
#ifdef COLLECT_STATISTICS
    ++gNumInterSection;
#endif
    ui li = lBegin;
    ui ri = rBegin;
    cn_count = 0;
    ui lCount = lEnd - lBegin, rCount = rEnd - rBegin;
    if (lCount / rCount > 50) {
        while (true) {
            if (rArray[ri] < lArray[li]) {
                ++ri;
                if (ri >= rEnd) return;
            }
            li = ComputeSetIntersection::GallopingSearch(lArray, li, lEnd, rArray[ri]);
            if (li >= lEnd) return;
            if (lArray[li] == rArray[ri]) {
                iters[cn_count][0] = li;
                iters[cn_count][1] = ri;
                cn[cn_count++] = lArray[li];
                ++li;
                ++ri;
                if (li >= lEnd || ri >= rEnd) {
                    return;
                }
            }
        }
    }
    else if (rCount / lCount > 50) {
        while (true) {
            if (lArray[li] < rArray[ri]) {
                ++li;
                if (li >= lEnd) {
                    return;
                }
            }
            ri = ComputeSetIntersection::GallopingSearch(rArray, ri, rEnd, lArray[li]);
            if (ri >= rEnd) return;
            if (lArray[li] == rArray[ri]) {
                iters[cn_count][0] = li;
                iters[cn_count][1] = ri;
                cn[cn_count++] = lArray[li];
                ++li;
                ++ri;
                if (li >= lEnd || ri >= rEnd) {
                    return;
                }
            }
        }
    }
    else {
        while (true) {
            if (lArray[li] < rArray[ri]) {
                ++li;
                if (li >= lEnd) {
                    return;
                }
            }
            else if (lArray[li] > rArray[ri]) {
                ++ri;
                if (ri >= rEnd) {
                    return;
                }
            }
            else {
                iters[cn_count][0] = li;
                iters[cn_count][1] = ri;
                cn[cn_count++] = lArray[li];
                ++li;
                ++ri;
                if (li >= lEnd || ri >= rEnd) {
                    return;
                }
            }
        }
    }
}

void multiIntersection(const VertexID **values, ui num, ui **iters, ui &iterSize,
                       std::vector<ui> &beginPoses, std::vector<ui> &endPoses, VertexID *cn, ui &cn_count, ui numBag) {
    if (num == 1) {
        cn_count = endPoses[0] - beginPoses[0];
        for (int i = 0; i < cn_count; ++i) {
            cn[i] = values[0][beginPoses[0] + i];
        }
        if (numBag == 1) {
            for (int i = 0; i < cn_count; ++i)
                iters[i][0] = beginPoses[0] + i;
        }
        return;
    }
#ifdef COLLECT_STATISTICS
    ++gNumInterSection;
#endif
    std::vector<ui> indices(num);
    std::vector<ui> reverseIndices(num);
    iterSize = 0;
    for (int i = 0; i < num; ++i) {
        if (endPoses[i] - beginPoses[i] == 0)
            return;
    }
    for (ui i = 0; i < num; ++i) {
        indices[i] = i;
    }
    for (int i = 0; i < num; ++i) {
        for (int j = i + 1; j < num; ++j) {
            if (endPoses[j] - beginPoses[j] < endPoses[i] - beginPoses[i]) {
                std::swap(values[i], values[j]);
                std::swap(indices[i], indices[j]);
                std::swap(beginPoses[i], beginPoses[j]);
                std::swap(endPoses[i], endPoses[j]);
            }
        }
    }
    for (ui i = 0; i < num; ++i) {
        reverseIndices[indices[i]] = i;
    }
    twoIntersection(values[0], beginPoses[0], endPoses[0], values[1], beginPoses[1], endPoses[1], iters, cn, cn_count);
    for (int i = 2; i < num; ++i) {
        if (cn_count == 0) return;
        ComputeSetIntersection::ComputeCandidates(cn, cn_count, values[i] + beginPoses[i], endPoses[i] - beginPoses[i], cn, cn_count);
    }
    if (cn_count > 0) {
        for (int i = 0; i < numBag; ++i)
            setIters(values[reverseIndices[i]], beginPoses[reverseIndices[i]], endPoses[reverseIndices[i]], iters, i, cn, cn_count);
    }
}

void
generateCandidates(const HyperTree &t, VertexID nID, int depth, const VertexID **arrays, ui *counts, ui num,
                   VertexID *cn, ui &cn_count, VertexID *partMatch, CandidateSpace &cs, VertexID current,
                   VertexID cartesianParent, std::vector<ui> &beginPoses, std::vector<ui> &endPoses) {
    const std::vector<VertexID> &largerAttrs = t.nodes[nID].largerAttrs[depth];
    const std::vector<VertexID> &smallerAttrs = t.nodes[nID].smallerAttrs[depth];
    VertexID maxCompared = 0, minCompared = (VertexID) - 1;
    for (int i = 0; i < largerAttrs.size(); ++i) {
        if (partMatch[largerAttrs[i]] > maxCompared)
            maxCompared = partMatch[largerAttrs[i]];
    }
    for (int i = 0; i < smallerAttrs.size(); ++i) {
        if (partMatch[smallerAttrs[i]] < minCompared)
            minCompared = partMatch[smallerAttrs[i]];
    }
    if (num != 0) {
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
}

void
generateCandidates(const HyperTree &t, VertexID nID, int depth, const VertexID **arrays, ui *counts, ui num,
                   ui **iters, ui &iterSize, VertexID *cn, ui &cn_count, ui numBags,
                   VertexID *partMatch, CandidateSpace &cs, VertexID current, VertexID cartesianParent,
                   std::vector<ui> &beginPoses, std::vector<ui> &endPoses) {
    const std::vector<VertexID> &largerAttrs = t.nodes[nID].largerAttrs[depth];
    const std::vector<VertexID> &smallerAttrs = t.nodes[nID].smallerAttrs[depth];
    VertexID maxCompared = 0, minCompared = (VertexID) - 1;
    for (int i = 0; i < largerAttrs.size(); ++i) {
        if (partMatch[largerAttrs[i]] > maxCompared)
            maxCompared = partMatch[largerAttrs[i]];
    }
    for (int i = 0; i < smallerAttrs.size(); ++i) {
        if (partMatch[smallerAttrs[i]] < minCompared)
            minCompared = partMatch[smallerAttrs[i]];
    }
    if (num != 0) {
        for (int i = 0; i < num; ++i) {
            beginPoses[i] = 0;
            endPoses[i] = counts[i];
        }
        if (!largerAttrs.empty()) {
            for (int i = 0; i < num; ++i) {
                beginPoses[i] = setBeginPos(arrays[i], counts[i], maxCompared);
            }
        }
        if (!smallerAttrs.empty()) {
            for (int i = 0; i < num; ++i) {
                endPoses[i] = setEndPos(arrays[i], counts[i], minCompared);
            }
        }
        for (int i = 0; i < num; ++i) {
            if (beginPoses[i] == endPoses[i]) {
                cn_count = iterSize = 0;
                return;
            }
        }
        if (num == 2) {
            twoIntersection(arrays[0], beginPoses[0], endPoses[0], arrays[1], beginPoses[1], endPoses[1], iters, cn, cn_count);
        }
        else {
            multiIntersection(arrays, num, iters, iterSize, beginPoses, endPoses, cn, cn_count, numBags);
        }
    }
    else if (cartesianParent != 99) {
        const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, current);
        optimizedCartesianProduct(cs, partMatch[cartesianParent], path, cn, cn_count, maxCompared, minCompared, partMatch);
        for (int i = 0; i < counts[0]; ++i) {
            iters[i][0] = i;
        }
    }
    else {
        // iterNum should be 0
        const VertexID *candSet = cs.candSet(current);
        counts[0] = (cs.candSize(current));
        for (int i = 0; i < cs.candSize(current); ++i) {
            cn[i] = candSet[i];
            iters[i][0] = i;
        }
        cn_count = iterSize;
    }
    iterSize = cn_count;
}

void setIters(const VertexID *values, ui begin, ui end, ui **iters, int id, VertexID *cn, ui cn_count) {
    ui pos = 0, pos2 = begin, iterSize = 0;
    if ((end - begin) / cn_count < BINARY_SEARCH_THRESHOLD) {
        while (pos < cn_count) {
            if (cn[pos] == values[pos2]) {
                iters[iterSize][id] = pos2;
                ++iterSize;
                ++pos;
            }
            else ++pos2;
        }
    }
    else {
        while (pos < cn_count) {
            pos2 = ComputeSetIntersection::BinarySearch(values, pos2, end, cn[pos]);
            iters[iterSize][id] = pos2;
            ++iterSize;
            ++pos;
        }
    }
}

bool generateCandidates(const DataGraph &g, VertexID v, VertexID **allCandidates, ui *allCandCount, const std::vector<int> &verticesAfter,
                        const std::vector<std::vector<VertexID>> &verticesLarger, const std::vector<std::vector<VertexID>> &verticesSmaller, const std::vector<int> &copyAfter,
                        const std::vector<int> &copyTypes, VertexID *partMatch) {
    for (int i = 0; i < verticesAfter.size(); ++i) {
        int pos = verticesAfter[i];
        ui neighborCount;
        VertexID *neighbors;
        neighbors = g.getNeighbors(v, neighborCount);
        ui beginPos = 0, endPos = neighborCount;
        VertexID maxCompared = 0, minCompared = (VertexID) - 1;
        for (VertexID u2: verticesLarger[i]) {
            if (partMatch[u2] > maxCompared)
                maxCompared = partMatch[u2];
        }
        for (VertexID u2: verticesSmaller[i]) {
            if (partMatch[u2] < minCompared)
                minCompared = partMatch[u2];
        }
        if (!verticesLarger[i].empty()) {
            beginPos = setBeginPos(neighbors, neighborCount, maxCompared);
        }
        if (!verticesSmaller[i].empty()) {
            endPos = setEndPos(neighbors, neighborCount, minCompared);
        }
        ComputeSetIntersection::ComputeCandidates(allCandidates[pos - 1], allCandCount[pos - 1], neighbors + beginPos,
                                                  endPos - beginPos, allCandidates[pos], allCandCount[pos]);
        if (allCandCount[pos] == 0) return false;
    }
    for (int i = 0; i < copyAfter.size(); ++i) {
        int pos = copyAfter[i];
        int type = copyTypes[i];
        ui neighborCount;
        VertexID *neighbors;
        if (type == 0) neighbors = g.getNeighbors(v, neighborCount);
        else if (type == 1) neighbors = g.getNeighborsLargerID(v, neighborCount);
        else neighbors = g.getNeighborsSmallerID(v, neighborCount);
        memcpy(allCandidates[pos], neighbors, neighborCount * sizeof(VertexID));
        allCandCount[pos] = neighborCount;
    }

    return true;
}

void nodeJoin(const HyperTree &t, VertexID nID, CandidateSpace &cs, bool *visited, VertexID *partMatch,
              int mappingSize, VertexID **allCandidates, ui *allCandCount, ui prefixSum, const VertexID **neighbors, ui *neighborCount, ui *allPoses,
              std::vector<std::vector<VertexID>> &tuples, TrieLevel &level, ui &length,
              std::vector<ui> &beginPoses, std::vector<ui> &endPoses) {
    VertexID **candidates = allCandidates + prefixSum;
    ui *candCount = allCandCount + prefixSum;
    ui *poses = allPoses + prefixSum;
    const HyperNode &tau = t.nodes[nID];
    const std::vector<VertexID> &trieOrder = t.trieOrder[nID];
    const VertexID *order = tau.attributes;
    if (mappingSize == tau.numAttributes) {
        if (t.trieOrder[nID].empty()) return;
        if (length < tuples.size()) {
            std::vector<VertexID> &tuple = tuples[length];
            for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                tuple[i] = partMatch[tau.attributes[trieOrder[i] + tau.prefixSize]];
            }
        }
        else {
            tuples.emplace_back(tau.numAttributes - tau.prefixSize);
            std::vector<VertexID> &tuple = tuples.back();
            for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                tuple[i] = partMatch[tau.attributes[trieOrder[i] + tau.prefixSize]];
            }
        }
        ++length;
        return;
    }
    const std::vector<std::vector<VertexID>> &vertexParents = tau.attributesBefore;
    const std::vector<int> &reuseCandidates = tau.candidatesBefore;
    // handle the first level
    VertexID u = order[mappingSize];
    if (mappingSize == 0) {
        candCount[0] = cs.candSize(u);
        memcpy(candidates[0], cs.candSet(u), cs.candSize(u) * sizeof(VertexID));
    }
    else {
        const std::vector<VertexID> &parents = vertexParents[mappingSize];
        int reuse = 0;
        if (reuseCandidates[mappingSize] != 0) {
            reuse = 1;
            neighbors[0] = allCandidates[reuseCandidates[mappingSize]];
            neighborCount[0] = allCandCount[reuseCandidates[mappingSize]];
        }
        for (int i = 0; i < parents.size(); ++i) {
            VertexID pU = parents[i];
            neighbors[i + reuse] = cs.getNeighbors(partMatch[pU], neighborCount[i + reuse]);
        }
        generateCandidates(t, nID, mappingSize, neighbors, neighborCount, parents.size() + reuse,
                           candidates[mappingSize], candCount[mappingSize], partMatch, cs, u, tau.cartesianParent[mappingSize],
                           beginPoses, endPoses);
    }
    if (trieOrder.size() == 1 && mappingSize == tau.numAttributes - 1) {
        level.values = candidates[mappingSize];
        level.length = candCount[mappingSize];
        return;
    }
    int depth = mappingSize;
    poses[depth] = 0;
    while (depth >= mappingSize) {
        while (poses[depth] < candCount[depth]) {
            VertexID v = candidates[depth][poses[depth]];
            ++poses[depth];
            if (visited[v]) continue;
            visited[v] = true;
            partMatch[order[depth]] = v;
            if (depth + 1 == tau.numAttributes) {
                if (length < tuples.size()) {
                    std::vector<VertexID> &tuple = tuples[length];
                    for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                        tuple[i] = partMatch[tau.attributes[trieOrder[i] + tau.prefixSize]];
                    }
                }
                else {
                    tuples.emplace_back(tau.numAttributes - tau.prefixSize);
                    std::vector<VertexID> &tuple = tuples.back();
                    for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                        tuple[i] = partMatch[tau.attributes[trieOrder[i] + tau.prefixSize]];
                    }
                }
                ++length;
                visited[partMatch[tau.attributes[depth]]] = false;
            }
            else {
                ++depth;
                poses[depth] = 0;
                u = order[depth];
                const std::vector<VertexID> &parents = vertexParents[depth];
                int reuse = 0;
                if (reuseCandidates[depth] != 0) {
                    reuse = 1;
                    neighbors[0] = allCandidates[reuseCandidates[depth]];
                    neighborCount[0] = allCandCount[reuseCandidates[depth]];
                }
                for (int i = 0; i < parents.size(); ++i) {
                    VertexID pU = parents[i];
                    neighbors[i + reuse] = cs.getNeighbors(partMatch[pU], neighborCount[i + reuse]);
                }
                generateCandidates(t, nID, depth, neighbors, neighborCount, parents.size() + reuse,
                                   candidates[depth], candCount[depth], partMatch, cs, u, tau.cartesianParent[depth],
                                   beginPoses, endPoses);
                if (trieOrder.size() == 1 && depth + 1 == tau.numAttributes) {
                    level.values = candidates[depth];
                    level.length = candCount[depth];
                    --depth;
                    continue;
                }
            }
        }
        --depth;
        if (depth >= mappingSize) visited[partMatch[order[depth]]] = false;
    }
}

void nodeJoin(const HyperTree &t, VertexID nID, const DataGraph &g, bool *visited, VertexID *partMatch,
              int mappingSize, VertexID **nodeCandidates, ui *nodeCandCount,
              VertexID **allCandidates, ui *allCandCount, ui prefixSum,
              const VertexID **neighbors, ui *neighborCount, ui *nodePoses,
              std::vector<std::vector<VertexID>> &tuples, TrieLevel &level, ui &length) {
    const auto &candIndex = t.candIndex[nID];
    VertexID **candidates = nodeCandidates + prefixSum;
    ui *candCount = nodeCandCount + prefixSum;
    ui *poses = nodePoses + prefixSum;
    const HyperNode &tau = t.nodes[nID];
    const std::vector<VertexID> &trieOrder = t.trieOrder[nID];
    const VertexID *order = tau.attributes;
    if (mappingSize == tau.numAttributes) {
        if (t.trieOrder[nID].empty()) return;
        if (length < tuples.size()) {
            std::vector<VertexID> &tuple = tuples[length];
            for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                tuple[i] = partMatch[tau.attributes[trieOrder[i] + tau.prefixSize]];
            }
        }
        else {
            tuples.emplace_back(tau.numAttributes - tau.prefixSize);
            std::vector<VertexID> &tuple = tuples.back();
            for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                tuple[i] = partMatch[tau.attributes[trieOrder[i] + tau.prefixSize]];
            }
        }
        ++length;
        return;
    }
    const std::vector<std::vector<int>> &verticesAfter = tau.attributesAfter;
    const std::vector<std::vector<std::vector<VertexID>>> &verticesLarger = tau.attrAfterLarger;
    const std::vector<std::vector<std::vector<VertexID>>> &verticesSmaller = tau.attrAfterSmaller;
    const std::vector<std::vector<int>> &copyAfter = tau.copyAfter;
    const std::vector<std::vector<int>> &copyTypes = tau.copyAfterTypes;
    const std::vector<std::vector<int>> &candidatesAfter = tau.candidatesAfter;
    // handle the first level
    VertexID u = order[mappingSize];
    if (mappingSize == 0) {
        candCount[0] = g.getNumVertices();
        for (VertexID v = 0; v < g.getNumVertices(); ++v)
            candidates[0][v] = v;
    }
    int depth = mappingSize;
    ui beginPos = 0, endPos = allCandCount[candIndex[depth]];
    if (!t.nodes[nID].largerAttrs[depth].empty()) {
        VertexID maxCompared = 0;
        for (VertexID u2: t.nodes[nID].largerAttrs[depth]) {
            if (partMatch[u2] > maxCompared)
                maxCompared = partMatch[u2];
        }
        beginPos = setBeginPos(candidates[depth], endPos, maxCompared);
    }
    if (!t.nodes[nID].smallerAttrs[depth].empty()) {
        VertexID minCompared = g.getNumVertices();
        for (VertexID u2: t.nodes[nID].smallerAttrs[depth]) {
            if (partMatch[u2] < minCompared)
                minCompared = partMatch[u2];
        }
        endPos = setEndPos(candidates[depth], endPos, minCompared);
    }
    poses[depth] = beginPos;
    candCount[depth] = endPos;
    for (int i = 0; i < candidatesAfter[depth].size(); ++i) {
        int pos = candidatesAfter[depth][i];
        allCandCount[pos] = endPos - beginPos;
        memcpy(allCandidates[pos], candidates[depth] + beginPos, allCandCount[pos] * sizeof(VertexID));
    }
    if (trieOrder.size() == 1 && mappingSize == tau.numAttributes - 1) {
        level.values = candidates[mappingSize];
        level.length = candCount[mappingSize];


        return;
    }
    poses[depth] = 0;
    while (depth >= mappingSize) {
        while (poses[depth] < candCount[depth]) {
            VertexID v = candidates[depth][poses[depth]];
            ++poses[depth];
            if (visited[v]) continue;
            partMatch[order[depth]] = v;
            if (!generateCandidates(g, v, allCandidates, allCandCount, verticesAfter[depth], verticesLarger[depth],
                                    verticesSmaller[depth], copyAfter[depth], copyTypes[depth], partMatch))
                continue;
            visited[v] = true;
            if (depth + 1 == tau.numAttributes) {
                if (length < tuples.size()) {
                    std::vector<VertexID> &tuple = tuples[length];
                    for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                        tuple[i] = partMatch[tau.attributes[trieOrder[i] + tau.prefixSize]];
                    }
                }
                else {
                    tuples.emplace_back(tau.numAttributes - tau.prefixSize);
                    std::vector<VertexID> &tuple = tuples.back();
                    for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                        tuple[i] = partMatch[tau.attributes[trieOrder[i] + tau.prefixSize]];
                    }
                }
                ++length;
                visited[partMatch[tau.attributes[depth]]] = false;
            }
            else {
                ++depth;
                beginPos = 0, endPos = allCandCount[candIndex[depth]];
                if (!t.nodes[nID].largerAttrs[depth].empty()) {
                    VertexID maxCompared = 0;
                    for (VertexID u2: t.nodes[nID].largerAttrs[depth]) {
                        if (partMatch[u2] > maxCompared)
                            maxCompared = partMatch[u2];
                    }
                    beginPos = setBeginPos(candidates[depth], endPos, maxCompared);
                }
                if (!t.nodes[nID].smallerAttrs[depth].empty()) {
                    VertexID minCompared = g.getNumVertices();
                    for (VertexID u2: t.nodes[nID].smallerAttrs[depth]) {
                        if (partMatch[u2] < minCompared)
                            minCompared = partMatch[u2];
                    }
                    endPos = setEndPos(candidates[depth], endPos, minCompared);
                }
                poses[depth] = beginPos;
                candCount[depth] = endPos;
                if (trieOrder.size() == 1 && depth + 1 == tau.numAttributes) {
                    level.values = candidates[depth] + beginPos;
                    level.length = candCount[depth];
                    --depth;
                    continue;
                }
                for (int i = 0; i < candidatesAfter[depth].size(); ++i) {
                    int pos = candidatesAfter[depth][i];
                    allCandCount[pos] = endPos - beginPos;
                    memcpy(allCandidates[pos], candidates[depth] + beginPos, allCandCount[pos] * sizeof(VertexID));
                }
            }
        }
        --depth;
        if (depth >= mappingSize) visited[partMatch[order[depth]]] = false;
    }
}

void globalJoin(std::vector<std::vector<VertexID>> &result, size_t &count, const Graph &query, const HyperTree &t,
                const DataGraph &g, CandidateSpace &cs, VertexID *partMatch, int mappingSize, int pathLength,
                VertexID **nodeCandidates, ui *nodeCandCount, VertexID **allCandidates, ui *allCandCount,
                ui prefixSum, bool *visited, const std::vector<std::vector<VertexID>> &nIDs,
                const std::vector<std::vector<VertexID>> &vertexParents, const std::vector<VertexID> &cartesianParent,
                ui ***iters, ui *iterSize, bool traversal, const VertexID **neighbors, ui *neighborCount,
                std::vector<std::vector<TrieLevel *>> &traversedLevels, std::vector<TrieLevel *>& lastLevels,
                ui *nodePoses, std::vector<ui> &traversePoses, std::vector<ui> &numBranches,
                TrieLevel &level, std::vector<std::vector<VertexID>> &tuples, ui &length,
                std::vector<ui> &beginPoses, std::vector<ui> &endPoses) {
    const HyperNode &globalNode = t.nodes[t.numNodes - 1];
    if (globalNode.numAttributes == mappingSize) {
#ifdef ALL_LEVEL
        if (traversal) enumerate(result, t, t.globalOrder, partMatch, mappingSize, traversedLevels, visited);
#else
        traverse(count, query, t, t.globalOrder, partMatch, mappingSize, lastLevels, traversedLevels, traversePoses, numBranches, visited, t.extendLevel, cs.labeled);
#endif
        return;
    }
    const auto &candIndex = t.candIndex[t.numNodes - 1];
    VertexID **candidates = nodeCandidates + prefixSum;
    ui *candCount = nodeCandCount + prefixSum;
    ui *poses = nodePoses + prefixSum;
    const std::vector<std::vector<int>> &candidatesAfter = globalNode.candidatesAfter;
    const std::vector<VertexID> &trieOrder = t.trieOrder.back();
    const VertexID *order = globalNode.attributes;
    for (int i = mappingSize; i < globalNode.numAttributes; ++i) poses[i] = 0;
    const std::vector<int> &reuseCandidates = t.nodes[t.numNodes - 1].candidatesBefore;
    // for each level and each nID, the trie node is traversedNodes[]
    // initialize the first level
    if (mappingSize == 0 && nIDs[0].empty()) {
        VertexID u = order[0];
        const VertexID *candSet = cs.candSet(u);
        for (int i = 0; i < cs.candSize(u); ++i) {
            candidates[0][i] = candSet[i];
        }
        iterSize[0] = cs.candSize(u);
    }
    else if (nIDs[mappingSize].empty()) {
        ui beginPos = 0, endPos = allCandCount[candIndex[pathLength]];
        if (!globalNode.largerAttrs[pathLength].empty()) {
            VertexID maxCompared = 0;
            for (VertexID u2: globalNode.largerAttrs[pathLength]) {
                if (partMatch[u2] > maxCompared)
                    maxCompared = partMatch[u2];
            }
            beginPos = setBeginPos(candidates[pathLength], endPos, maxCompared);
        }
        if (!globalNode.smallerAttrs[pathLength].empty()) {
            VertexID minCompared = g.getNumVertices();
            for (VertexID u2: globalNode.smallerAttrs[pathLength]) {
                if (partMatch[u2] < minCompared)
                    minCompared = partMatch[u2];
            }
            endPos = setEndPos(candidates[pathLength], endPos, minCompared);
        }
        poses[pathLength] = beginPos;
        candCount[pathLength] = iterSize[pathLength] = endPos;
        for (int i = 0; i < candidatesAfter[pathLength].size(); ++i) {
            int pos = candidatesAfter[pathLength][i];
            allCandCount[pos] = endPos - beginPos;
            memcpy(allCandidates[pos], candidates[pathLength] + beginPos, allCandCount[pos] * sizeof(VertexID));
        }
    }
    else {
        for (int i = 0; i < nIDs[mappingSize].size(); ++i) {
            VertexID id = nIDs[mappingSize][i];
            neighbors[i] = traversedLevels[id].back()->values;
            neighborCount[i] = traversedLevels[id].back()->length;
        }
        for (int i = 0; i < vertexParents[mappingSize].size(); ++i) {
            VertexID pU = vertexParents[mappingSize][i];
            VertexID pV = partMatch[pU];
            VertexID u = order[mappingSize];
            ui numNeighbors;
            const VertexID *pVNeighbors = cs.getNeighbors(pV, numNeighbors);
            neighbors[nIDs[mappingSize].size() + i] = pVNeighbors;
            neighborCount[nIDs[mappingSize].size() + i] = numNeighbors;
        }
        if (nIDs[mappingSize].size() + vertexParents[mappingSize].size() == 0 && cartesianParent[mappingSize] == 99) {
            neighbors[0] = cs.candSet(order[0]);
            neighborCount[0] = cs.candSize(order[0]);
        }
        generateCandidates(t, t.numNodes - 1, mappingSize, neighbors, neighborCount, nIDs[mappingSize].size() +
                                                                                     vertexParents[mappingSize].size(),
                           iters[pathLength], iterSize[pathLength], candidates[pathLength], candCount[pathLength], nIDs[mappingSize].size(),
                           partMatch, cs, order[mappingSize], cartesianParent[mappingSize], beginPoses, endPoses);
    }
    if (mappingSize == globalNode.numAttributes - 1 && mappingSize == t.extendLevel) {
        level.values = candidates[pathLength];
        level.length = candCount[pathLength];
        traverse(count, query, t, t.globalOrder, partMatch, mappingSize, lastLevels, traversedLevels, traversePoses, numBranches, visited,
                 t.extendLevel, cs.labeled);
        return;
    }
    const std::vector<std::vector<int>> &verticesAfter = globalNode.attributesAfter;
    const std::vector<std::vector<std::vector<VertexID>>> &verticesLarger = globalNode.attrAfterLarger;
    const std::vector<std::vector<std::vector<VertexID>>> &verticesSmaller = globalNode.attrAfterSmaller;
    const std::vector<std::vector<int>> &copyAfter = globalNode.copyAfter;
    const std::vector<std::vector<int>> &copyTypes = globalNode.copyAfterTypes;
    int depth = 0;
    while (depth >= 0) {
        while (poses[mappingSize + depth] < iterSize[pathLength + depth]) {
            VertexID v;
            if (nIDs[mappingSize + depth].empty()) v = candidates[pathLength + depth][poses[mappingSize + depth]];
            else {
                ui childPos = iters[pathLength + depth][poses[mappingSize + depth]][0];
                VertexID nID = nIDs[mappingSize + depth][0];
                v = traversedLevels[nID].back()->values[childPos];
            }
            ++poses[mappingSize + depth];
            if (visited[v]) continue;
            partMatch[order[mappingSize + depth]] = v;
            if (!generateCandidates(g, v, allCandidates, allCandCount, verticesAfter[pathLength + depth],
                                    verticesLarger[pathLength + depth],
                                    verticesSmaller[pathLength + depth], copyAfter[pathLength + depth],
                                    copyTypes[pathLength + depth], partMatch))
                continue;
            visited[v] = true;
            for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                VertexID nID = nIDs[mappingSize + depth][i];
                ui childPos = iters[pathLength + depth][poses[mappingSize + depth]-1][i];
                if (traversedLevels[nID].back()->children) {
                    traversedLevels[nID].push_back(&traversedLevels[nID].back()->children[childPos]);
                    lastLevels[nID] = traversedLevels[nID].back();
                }
                else {
                    traversedLevels[nID].push_back(nullptr);
                    lastLevels[nID] = nullptr;
                }
            }
            if (depth + 1 + mappingSize == globalNode.numAttributes) {
#ifdef ALL_LEVEL
                if (traversal) enumerate(result, t, t.globalOrder, partMatch, mappingSize + depth + 1, traversedLevels, visited);
#else
                if (t.extendLevel != globalNode.numAttributes) {
                    if (length < tuples.size()) {
                        std::vector<VertexID> &tuple = tuples[length];
                        for (int i = 0; i < globalNode.numAttributes - t.extendLevel; ++i) {
                            tuple[i] = partMatch[globalNode.attributes[trieOrder[i] + t.extendLevel]];
                        }
                    }
                    else {
                        tuples.emplace_back(globalNode.numAttributes - t.extendLevel);
                        std::vector<VertexID> &tuple = tuples.back();
                        for (int i = 0; i < globalNode.numAttributes - t.extendLevel; ++i) {
                            tuple[i] = partMatch[globalNode.attributes[trieOrder[i] + t.extendLevel]];
                        }
                    }
                    ++length;
                }
                else {
                    // when a match of shared attributes is found, traverse the remaining attributes
                    traverse(count, query, t, t.globalOrder, partMatch, mappingSize + depth + 1, lastLevels, traversedLevels,
                             traversePoses, numBranches, visited,
                             t.extendLevel, cs.labeled);
                }
#endif
                for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                    VertexID nID = nIDs[mappingSize + depth][i];
                    traversedLevels[nID].pop_back();
                    lastLevels[nID] = traversedLevels[nID].back();
                }
                visited[v] = false;
            }
            else {
                ++depth;
                poses[mappingSize + depth] = 0;
                ui beginPos = 0;
                if (!nIDs[mappingSize + depth].empty()) {
                    for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                        VertexID nID = nIDs[mappingSize + depth][i];
                        neighbors[i] = lastLevels[nID]->values;
                        neighborCount[i] = lastLevels[nID]->length;
                    }
                    for (int i = 0; i < vertexParents[mappingSize + depth].size(); ++i) {
                        VertexID pU = vertexParents[mappingSize + depth][i];
                        VertexID pV = partMatch[pU];
                        VertexID u = order[mappingSize + depth];
                        ui numNeighbors;
                        const VertexID *pVNeighbors = cs.getNeighbors(pV, numNeighbors);
                        neighbors[nIDs[mappingSize + depth].size() + i] = pVNeighbors;
                        neighborCount[nIDs[mappingSize + depth].size() + i] = numNeighbors;
                    }

                    generateCandidates(t, t.numNodes - 1, mappingSize + depth, neighbors, neighborCount,
                                       nIDs[mappingSize + depth].size() + vertexParents[mappingSize + depth].size(),
                                       iters[pathLength + depth], iterSize[pathLength + depth], candidates[pathLength + depth], candCount[pathLength + depth],
                                       nIDs[mappingSize + depth].size(), partMatch,
                                       cs, order[mappingSize + depth], cartesianParent[mappingSize + depth],
                                       beginPoses, endPoses);
                }
                else {
                    ui endPos = allCandCount[candIndex[pathLength + depth]];
                    if (!globalNode.largerAttrs[pathLength + depth].empty()) {
                        VertexID maxCompared = 0;
                        for (VertexID u2: globalNode.largerAttrs[pathLength + depth]) {
                            if (partMatch[u2] > maxCompared)
                                maxCompared = partMatch[u2];
                        }
                        beginPos = setBeginPos(candidates[pathLength + depth], endPos, maxCompared);
                    }
                    if (!globalNode.smallerAttrs[pathLength + depth].empty()) {
                        VertexID minCompared = g.getNumVertices();
                        for (VertexID u2: globalNode.smallerAttrs[pathLength + depth]) {
                            if (partMatch[u2] < minCompared)
                                minCompared = partMatch[u2];
                        }
                        endPos = setEndPos(candidates[pathLength + depth], endPos, minCompared);
                    }
                    poses[pathLength + depth] = beginPos;
                    candCount[pathLength + depth] = iterSize[pathLength + depth] = endPos;
                    for (int i = 0; i < candidatesAfter[pathLength + depth].size(); ++i) {
                        int pos = candidatesAfter[pathLength + depth][i];
                        allCandCount[pos] = endPos - beginPos;
                        memcpy(allCandidates[pos], candidates[pathLength + depth] + beginPos, allCandCount[pos] * sizeof(VertexID));
                    }
                }
                if (level.oneLevel && depth + 1 + mappingSize == globalNode.numAttributes && t.extendLevel == globalNode.numAttributes - 1) {
                    level.values = candidates[pathLength + depth] + beginPos;
                    level.length = iterSize[pathLength + depth];
                    traverse(count, query, t, t.globalOrder, partMatch, mappingSize + depth, lastLevels, traversedLevels, traversePoses, numBranches, visited,
                             t.extendLevel, cs.labeled);
                    break;
                }
            }
        }
        --depth;
#ifndef ALL_LEVEL
        if (depth + 1 + mappingSize == t.extendLevel && !tuples.empty()) {
            if (!level.oneLevel) buildTrie(tuples, level, length);
            traverse(count, query, t, t.globalOrder, partMatch, mappingSize + depth + 1, lastLevels, traversedLevels, traversePoses, numBranches, visited,
                     t.extendLevel, cs.labeled);
        }
#endif
        if (depth >= 0) {
            VertexID v = partMatch[order[mappingSize + depth]];
            for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                VertexID nID = nIDs[mappingSize + depth][i];
                traversedLevels[nID].pop_back();
                lastLevels[nID] = traversedLevels[nID].back();
            }
            visited[v] = false;
        }
    }
}

void sharedJoin(const HyperTree &t, const PrefixNode *pt, const Graph &query, CandidateSpace &cs,
                std::vector<TrieLevel> &trieLevels, bool *visited, std::vector<std::vector<VertexID>> &result,
                size_t &count, std::vector<ui> &beginPoses, std::vector<ui> &endPoses, bool traverse) {
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
    const HyperNode &globalNode = t.nodes[t.numNodes - 1];
    const PrefixNode *pn = pt;
    ui maxSize = cs.getMaxSize();
    ui ***iters = new ui **[t.numAttributes];
    ui *iterSizes = new ui [t.numAttributes];
    for (int i = 0; i < t.numAttributes; ++i) {
        iters[i] = new ui *[maxSize];
        for (int j = 0; j < maxSize; ++j) {
            iters[i][j] = new ui[t.numAttributes];
        }
        iterSizes[i] = 0;
    }
    for (int i = 0; i < prefixSum.back(); ++i) {
        nodeCandidates[i] = new VertexID[cs.getMaxSize()];
        nodeCandCount[i] = 0;
        nodePoses[i] = 0;
    }
    std::vector<ui> traversePoses(t.numAttributes);
    std::vector<ui> numBranches(t.numAttributes);
    ui height = pt->getHeight();
    VertexID **pCandidates = new VertexID *[height];
    ui *pCandCount = new ui[height];
    pCandidates[0] = new VertexID[cs.getMaxSize()];
    for (int i = 0; i < height; ++i) {
        pCandCount[i] = 0;
    }
    std::vector<ui> pPoses(height, 0);
    std::vector<ui> childPoses(height, 0);
    const VertexID **neighbors = new const VertexID *[t.numAttributes];
    ui *neighborCount = new ui[t.numAttributes];
    ui maxEdgeSize = maxNumBackWard(t.globalOrder, query);
    ui maxNum = height + globalNode.numAttributes - globalNode.prefixSize;
    std::vector<std::vector<TrieLevel *>> traversedLevels(trieLevels.size());
    std::vector<TrieLevel *> lastLevels(trieLevels.size());
    for (VertexID nID = 0; nID < trieLevels.size(); ++nID) {
        traversedLevels[nID].push_back(&trieLevels[nID]);
        lastLevels[nID] = &trieLevels[nID];
    }
    pn = pt;
    std::vector<const PrefixNode *> nodes(height, nullptr);
    std::vector<int> mappingSizes = getMappingSizes(t, pt);
    std::map<const PrefixNode *, std::vector<VertexID>> bagsBelow;
    pt->addBagsBelow(bagsBelow, t.numNodes - 1);
    int depth = 0;
    for (VertexID nID: pn -> nIDsToCall) {
        if (nID != t.numNodes - 1) {
            nodeJoin(t, nID, cs, visited, partMatch, 0, nodeCandidates, nodeCandCount, prefixSum[nID], neighbors, neighborCount,
                     nodePoses, tuples[nID], trieLevels[nID], lengths[nID], beginPoses, endPoses);
        }
    }

    // create candidates for the first prefixNode
    if (!pn -> children.empty()) {
        if (pn -> children[0] -> pathToGlobal) {
            for (VertexID nID : pn -> nIDsToBuild) {
                if (lengths[nID] == 0) return;
                if (!trieLevels[nID].oneLevel) buildTrie(tuples[nID], trieLevels[nID], lengths[nID]);
            }
            if (pn-> children[0] -> nIDsToJoin.empty()) {
                VertexID u = pn -> children[0] -> u;
                pCandCount[0] = iterSizes[0] = cs.candSize(u);
                memcpy(pCandidates[0], cs.candSet(u), pCandCount[0] * sizeof(VertexID));
            }
            else {
                for (int i = 0; i < pn -> children[0] -> nIDsToJoin.size(); ++i) {
                    VertexID id = pn -> children[0] -> nIDsToJoin[i];
                    neighbors[i] = trieLevels[id].values;
                    neighborCount[i] = trieLevels[id].length;
                }
                generateCandidates(t, bagsBelow[pn][0], 0, neighbors, neighborCount, pn->children[0]->nIDsToJoin.size(), iters[0],
                                   iterSizes[0], pCandidates[0], pCandCount[0], pn-> children[0] -> nIDsToJoin.size(), partMatch, cs, pn->children[0]->u, 0,
                                   beginPoses, endPoses);
            }
        }
        else {
            pCandCount[0] = cs.candSize(pn->children[0]->u);
            memcpy(pCandidates[0], cs.candSet(pn->children[0]->u), pCandCount[0] * sizeof(VertexID));
        }
    }
    while (depth >= 0) {
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
                                 nodeCandCount, prefixSum[nID], neighbors, neighborCount, nodePoses, tuples[nID],
                                 trieLevels[nID], lengths[nID], beginPoses, endPoses);
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
                        // all joined relations extend one level
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
                                     nodePoses, tuples[nID], trieLevels[nID], lengths[nID], beginPoses, endPoses);
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
                                    if (lengths[nID2] == 0) {
                                        flag = false;
                                        break;
                                    }
                                }
                                if (flag) {
                                    for (VertexID nID2 : current -> nIDsToBuild)
                                        if (!trieLevels[nID2].oneLevel) buildTrie(tuples[nID2], trieLevels[nID2], lengths[nID2]);
                                    globalJoin(result, count, query, t, cs, partMatch, mappingSizes[nID], depth + 1,
                                               nodeCandidates, nodeCandCount, prefixSum[nID],
                                               visited, globalNode.nIDs, globalNode.attributesBefore,
                                               globalNode.cartesianParent, iters, iterSizes, traverse, neighbors, neighborCount,
                                               traversedLevels, lastLevels, nodePoses, traversePoses, numBranches,
                                               trieLevels[nID], tuples[nID], lengths[nID], beginPoses, endPoses);
                                }
                                else {
                                    for (VertexID nID2 : current -> nIDsToBuild) {
//                                        tuples[nID2].clear();
                                        lengths[nID2] = 0;
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
                    if (lengths[nID] == 0) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    for (VertexID nID : pn -> nIDsToBuild)
                        if (!trieLevels[nID].oneLevel) buildTrie(tuples[nID], trieLevels[nID], lengths[nID]);
                }
                else {
                    for (VertexID nID2 : pn -> nIDsToBuild) {
//                        tuples[nID2].clear();
                        lengths[nID2] = 0;
                    }
                    visited[partMatch[pn->u]] = false;
                    for (VertexID nID : pn -> nIDsToJoin) {
                        traversedLevels[nID].pop_back();
                        lastLevels[nID] = traversedLevels[nID].back();
                    }
                    --depth;
                    if (depth == 0) pn = pt;
                    else pn = nodes[depth - 1];
                    continue;
                }
                VertexID firstBag = bagsBelow[current][0];
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
                        neighbors[i + reuse] = cs.getNeighbors(partMatch[pU], neighborCount[i + reuse]);
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
                        const VertexID *pVNeighbors = cs.getNeighbors(pV, numNeighbors);
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
                    VertexID firstBag = bagsBelow[current][0];
                    const std::vector<VertexID> &parents = t.nodes[firstBag].attributesBefore[depth];
                    int reuse = 0;
                    if (t.nodes[firstBag].candidatesBefore[depth] != 0) {
                        reuse = 1;
                        neighbors[0] = nodeCandidates[t.nodes[firstBag].candidatesBefore[depth]];
                        neighborCount[0] = nodeCandCount[t.nodes[firstBag].candidatesBefore[depth]];
                    }
                    for (int i = 0; i < parents.size(); ++i) {
                        VertexID pU = parents[i];
                        neighbors[i + reuse] = cs.getNeighbors(partMatch[pU], neighborCount[i + reuse]);
                    }
                    generateCandidates(t, firstBag, depth, neighbors, neighborCount, parents.size() + reuse,
                                       nodeCandidates[prefixSum[firstBag] + depth], nodeCandCount[prefixSum[firstBag] + depth],
                                       partMatch, cs, u, current->cartesianParent, beginPoses, endPoses);
                    pCandidates[depth] = nodeCandidates[prefixSum[firstBag] + depth];
                    pCandCount[depth] = nodeCandCount[prefixSum[firstBag] + depth];
                }
                else {
                    pCandCount[0] = cs.candSize(u);
                    memcpy(pCandidates[0], cs.candSet(u), pCandCount[0] * sizeof(VertexID));
                }
            }
        }
        --depth;
        if (depth >= 0) {
            const PrefixNode *current = pn;
            if (depth == 0) pn = pt;
            else pn = nodes[depth - 1];
            VertexID u = current -> u;
            VertexID v = partMatch[u];
            for (VertexID nID: current -> nIDsToCall) {
                if (nID == t.numNodes - 1) {
                    bool flag = true;
                    for (VertexID nID2 : current -> nIDsToBuild) {
                        if (lengths[nID2] == 0) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        for (VertexID nID2 : current -> nIDsToBuild)
                            if (!trieLevels[nID2].oneLevel) buildTrie(tuples[nID2], trieLevels[nID2], lengths[nID2]);
                        globalJoin(result, count, query, t, cs, partMatch, mappingSizes[nID], depth + 1, nodeCandidates, nodeCandCount,
                                   prefixSum[nID], visited, globalNode.nIDs, globalNode.attributesBefore,
                                   globalNode.cartesianParent, iters, iterSizes, traverse, neighbors, neighborCount,
                                   traversedLevels, lastLevels, nodePoses, traversePoses, numBranches,
                                   trieLevels[nID], tuples[nID], lengths[nID], beginPoses, endPoses);
                    }
                    else {
                        for (VertexID nID2 : current -> nIDsToBuild) {
                            lengths[nID2] = 0;
                        }
                    }
                }
            }
            visited[v] = false;
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
                if (lengths[nID2] == 0) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                for (VertexID nID2 : pn -> nIDsToBuild)
                    if (!trieLevels[nID2].oneLevel) buildTrie(tuples[nID2], trieLevels[nID2], lengths[nID2]);
                globalJoin(result, count, query, t, cs, partMatch, mappingSizes[nID], depth + 1, nodeCandidates, nodeCandCount,
                           prefixSum[nID], visited, globalNode.nIDs, globalNode.attributesBefore,
                           globalNode.cartesianParent, iters, iterSizes, traverse, neighbors, neighborCount,
                           traversedLevels, lastLevels, nodePoses, traversePoses, numBranches,
                           trieLevels[nID], tuples[nID], lengths[nID], beginPoses, endPoses);
            }
            else {
                for (VertexID nID2 : pn -> nIDsToBuild) {
//                    tuples[nID2].clear();
                    lengths[nID2] = 0;
                }
            }
        }
    }
    delete[] partMatch;
    delete[] neighbors;
    for (int i = 0; i < trieLevels.size(); ++i) {
        delete[] nodeCandidates[i];
    }
    delete[] nodeCandidates;
    delete[] nodeCandCount;
    delete[] nodePoses;
    for (int i = 0; i < t.numAttributes; ++i) {
        for (int j = 0; j < maxSize; ++j) {
            delete[] iters[i][j];
        }
        delete[] iters[i];
    }
    delete[] iters;
    delete[] iterSizes;
}

void globalJoin(std::vector<std::vector<VertexID>> &result, size_t &count, const Graph &query, const HyperTree &t,
                CandidateSpace &cs, VertexID *partMatch, int mappingSize, int pathLength, VertexID **allCandidates, ui *allCandCount,
                ui prefixSum, bool *visited, const std::vector<std::vector<VertexID>> &nIDs,
                const std::vector<std::vector<VertexID>> &vertexParents, const std::vector<VertexID> &cartesianParent,
                ui ***iters, ui *iterSize, bool traversal, const VertexID **neighbors, ui *neighborCount,
                std::vector<std::vector<TrieLevel *>> &traversedLevels, std::vector<TrieLevel *>& lastLevels,
                ui *allPoses, std::vector<ui> &traversePoses, std::vector<ui> &numBranches,
                TrieLevel &level, std::vector<std::vector<VertexID>> &tuples, ui &length,
                std::vector<ui> &beginPoses, std::vector<ui> &endPoses) {
    const HyperNode &globalNode = t.nodes[t.numNodes - 1];
    if (globalNode.numAttributes == mappingSize) {
#ifdef ALL_LEVEL
        if (traversal) enumerate(result, t, t.globalOrder, partMatch, mappingSize, traversedLevels, visited);
#else
        traverse(count, query, t, t.globalOrder, partMatch, mappingSize, lastLevels, traversedLevels, traversePoses, numBranches, visited, t.extendLevel, cs.labeled);
#endif
        return;
    }
    VertexID **candidates = allCandidates + prefixSum;
    ui *candCount = allCandCount + prefixSum;
    ui *poses = allPoses + prefixSum;
    const std::vector<VertexID> &trieOrder = t.trieOrder.back();
    const VertexID *order = globalNode.attributes;
    for (int i = mappingSize; i < globalNode.numAttributes; ++i) poses[i] = 0;
    const std::vector<int> &reuseCandidates = t.nodes[t.numNodes - 1].candidatesBefore;
    // for each level and each nID, the trie node is traversedNodes[]
    // initialize the first level
    if (mappingSize == 0 && nIDs[0].empty()) {
        VertexID u = order[0];
        const VertexID *candSet = cs.candSet(u);
        for (int i = 0; i < cs.candSize(u); ++i) {
            candidates[0][i] = candSet[i];
        }
        iterSize[0] = cs.candSize(u);
    }
    else if (nIDs[mappingSize].empty()) {
        VertexID u = order[mappingSize];
        const std::vector<VertexID> &parents = vertexParents[mappingSize];
        int reuse = 0;
        if (reuseCandidates[mappingSize] != 0) {
            reuse = 1;
            neighbors[0] = allCandidates[reuseCandidates[mappingSize]];
            neighborCount[0] = allCandCount[reuseCandidates[mappingSize]];
        }
        for (int i = 0; i < parents.size(); ++i) {
            VertexID pU = parents[i];
            neighbors[i + reuse] = cs.getNeighbors(partMatch[pU], neighborCount[i + reuse]);
        }
        generateCandidates(t, t.numNodes - 1, mappingSize, neighbors, neighborCount, parents.size() + reuse,
                           candidates[pathLength], candCount[pathLength], partMatch, cs, u, t.nodes[t.numNodes - 1].cartesianParent[mappingSize],
                           beginPoses, endPoses);
        iterSize[pathLength] = candCount[pathLength];
    }
    else {
        for (int i = 0; i < nIDs[mappingSize].size(); ++i) {
            VertexID id = nIDs[mappingSize][i];
            neighbors[i] = traversedLevels[id].back()->values;
            neighborCount[i] = traversedLevels[id].back()->length;
        }
        for (int i = 0; i < vertexParents[mappingSize].size(); ++i) {
            VertexID pU = vertexParents[mappingSize][i];
            VertexID pV = partMatch[pU];
            VertexID u = order[mappingSize];
            ui numNeighbors;
            const VertexID *pVNeighbors = cs.getNeighbors(pV, numNeighbors);
            neighbors[nIDs[mappingSize].size() + i] = pVNeighbors;
            neighborCount[nIDs[mappingSize].size() + i] = numNeighbors;
        }
        if (nIDs[mappingSize].size() + vertexParents[mappingSize].size() == 0 && cartesianParent[mappingSize] == 99) {
            neighbors[0] = cs.candSet(order[0]);
            neighborCount[0] = cs.candSize(order[0]);
        }
        generateCandidates(t, t.numNodes - 1, mappingSize, neighbors, neighborCount, nIDs[mappingSize].size() +
                                                                                 vertexParents[mappingSize].size(),
                           iters[pathLength], iterSize[pathLength], candidates[pathLength], candCount[pathLength], nIDs[mappingSize].size(),
                           partMatch, cs, order[mappingSize], cartesianParent[mappingSize], beginPoses, endPoses);
    }
    if (mappingSize == globalNode.numAttributes - 1 && mappingSize == t.extendLevel) {
        level.values = candidates[pathLength];
        level.length = candCount[pathLength];
        traverse(count, query, t, t.globalOrder, partMatch, mappingSize, lastLevels, traversedLevels, traversePoses, numBranches, visited,
                 t.extendLevel, cs.labeled);
        return;
    }
    int depth = 0;
    while (depth >= 0) {
        while (poses[mappingSize + depth] < iterSize[pathLength + depth]) {
            VertexID v;
            if (nIDs[mappingSize + depth].empty()) v = candidates[pathLength + depth][poses[mappingSize + depth]];
            else {
                ui childPos = iters[pathLength + depth][poses[mappingSize + depth]][0];
                VertexID nID = nIDs[mappingSize + depth][0];
                v = traversedLevels[nID].back()->values[childPos];
                // all joined relations extend one level
                if (!visited[v]) {
                    for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                        nID = nIDs[mappingSize + depth][i];
                        childPos = iters[pathLength + depth][poses[mappingSize + depth]][i];
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
            ++poses[mappingSize + depth];
            if (visited[v]) continue;
            visited[v] = true;
            partMatch[order[mappingSize + depth]] = v;
            if (depth + 1 + mappingSize == globalNode.numAttributes) {
#ifdef ALL_LEVEL
                if (traversal) enumerate(result, t, t.globalOrder, partMatch, mappingSize + depth + 1, traversedLevels, visited);
#else
                if (t.extendLevel != globalNode.numAttributes) {
                    if (length < tuples.size()) {
                        std::vector<VertexID> &tuple = tuples[length];
                        for (int i = 0; i < globalNode.numAttributes - t.extendLevel; ++i) {
                            tuple[i] = partMatch[globalNode.attributes[trieOrder[i] + t.extendLevel]];
                        }
                    }
                    else {
                        tuples.emplace_back(globalNode.numAttributes - t.extendLevel);
                        std::vector<VertexID> &tuple = tuples.back();
                        for (int i = 0; i < globalNode.numAttributes - t.extendLevel; ++i) {
                            tuple[i] = partMatch[globalNode.attributes[trieOrder[i] + t.extendLevel]];
                        }
                    }
                    ++length;
                }
                else {
                    // when a match of shared attributes is found, traverse the remaining attributes
                    traverse(count, query, t, t.globalOrder, partMatch, mappingSize + depth + 1, lastLevels, traversedLevels,
                             traversePoses, numBranches, visited,
                             t.extendLevel, cs.labeled);
                }
#endif
                for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                    VertexID nID = nIDs[mappingSize + depth][i];
                    traversedLevels[nID].pop_back();
                    lastLevels[nID] = traversedLevels[nID].back();
                }
                visited[v] = false;
            }
            else {
                ++depth;
                poses[mappingSize + depth] = 0;
                if (nIDs[mappingSize + depth].empty()) {
                    VertexID u = order[mappingSize + depth];
                    const std::vector<VertexID> &parents = vertexParents[mappingSize + depth];
                    int reuse = 0;
                    if (reuseCandidates[mappingSize + depth] != 0) {
                        reuse = 1;
                        neighbors[0] = allCandidates[reuseCandidates[mappingSize + depth]];
                        neighborCount[0] = allCandCount[reuseCandidates[mappingSize + depth]];
                    }
                    for (int i = 0; i < parents.size(); ++i) {
                        VertexID pU = parents[i];
                        neighbors[i + reuse] = cs.getNeighbors(partMatch[pU], neighborCount[i + reuse]);
                    }
                    generateCandidates(t, t.numNodes - 1, mappingSize + depth, neighbors, neighborCount, parents.size() + reuse,
                                       candidates[pathLength + depth], candCount[pathLength + depth], partMatch, cs, u,
                                       t.nodes[t.numNodes - 1].cartesianParent[mappingSize + depth], beginPoses, endPoses);
                    iterSize[pathLength + depth] = candCount[pathLength + depth];
                }
                else {
                    for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                        VertexID nID = nIDs[mappingSize + depth][i];
                        neighbors[i] = lastLevels[nID]->values;
                        neighborCount[i] = lastLevels[nID]->length;
                    }
                    for (int i = 0; i < vertexParents[mappingSize + depth].size(); ++i) {
                        VertexID pU = vertexParents[mappingSize + depth][i];
                        VertexID pV = partMatch[pU];
                        VertexID u = order[mappingSize + depth];
                        ui numNeighbors;
                        const VertexID *pVNeighbors = cs.getNeighbors(pV, numNeighbors);
                        neighbors[nIDs[mappingSize + depth].size() + i] = pVNeighbors;
                        neighborCount[nIDs[mappingSize + depth].size() + i] = numNeighbors;
                    }

                    generateCandidates(t, t.numNodes - 1, mappingSize + depth, neighbors, neighborCount,
                                       nIDs[mappingSize + depth].size() + vertexParents[mappingSize + depth].size(),
                                       iters[pathLength + depth], iterSize[pathLength + depth], candidates[pathLength + depth], candCount[pathLength + depth],
                                       nIDs[mappingSize + depth].size(), partMatch,
                                       cs, order[mappingSize + depth], cartesianParent[mappingSize + depth],
                                       beginPoses, endPoses);
                }
                if (level.oneLevel && depth + 1 + mappingSize == globalNode.numAttributes && t.extendLevel == globalNode.numAttributes - 1) {
                    level.values = candidates[pathLength + depth];
                    level.length = iterSize[pathLength + depth];
                    traverse(count, query, t, t.globalOrder, partMatch, mappingSize + depth, lastLevels, traversedLevels, traversePoses, numBranches, visited,
                             t.extendLevel, cs.labeled);
                    break;
                }
            }
        }
        --depth;
#ifndef ALL_LEVEL
        if (depth + 1 + mappingSize == t.extendLevel && !tuples.empty()) {
            if (!level.oneLevel) buildTrie(tuples, level, length);
            traverse(count, query, t, t.globalOrder, partMatch, mappingSize + depth + 1, lastLevels, traversedLevels, traversePoses, numBranches, visited,
                     t.extendLevel, cs.labeled);
        }
#endif
        if (depth >= 0) {
            VertexID v = partMatch[order[mappingSize + depth]];
            for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                VertexID nID = nIDs[mappingSize + depth][i];
                traversedLevels[nID].pop_back();
                lastLevels[nID] = traversedLevels[nID].back();
            }
            visited[v] = false;
        }
    }
}

// for enumerating all results
void enumerate(std::vector<std::vector<VertexID>> &result, const HyperTree &t, const std::vector<VertexID> &order,
               VertexID *partMatch, int mappingSize, std::vector<std::vector<TrieLevel *>> &traversedLevels, bool *visited) {
    if (t.defaultPartition.size() > mappingSize) mappingSize = t.defaultPartition.size();
    ui numNodes = traversedLevels.size();
    ui totalLevel = 0;
    std::vector<ui> numLevels(numNodes, 0);
    for (int i = 0; i < numNodes; ++i) {
        ui level = 0;
        if (traversedLevels[i].empty()) continue;
        const TrieLevel* l = traversedLevels[i].back();
        while (l && l->length) {
            l = &(l->children[0]);
            ++level;
        }
        totalLevel += level;
        numLevels[i] = level;
    }
    if (totalLevel == 0) {
        produceResult(result, partMatch, t.numAttributes);
        return;
    }
    ui matched = t.numAttributes - totalLevel;
    std::vector<ui> poses(totalLevel, 0);
    std::vector<ui> counts(totalLevel, 0);
    std::vector<VertexID> maxCompared(totalLevel, 0), minCompared(totalLevel, (VertexID) - 1);
    std::vector<ui> beginPoses(totalLevel), endPoses(totalLevel);
    // initialize the first level
    VertexID firstLevelNID = t.nIDs[mappingSize][0];
    for (int i = 0; i < t.largerAttrs[matched].size(); ++i) {
        if (partMatch[t.largerAttrs[matched][i]] > maxCompared[matched])
            maxCompared[0] = partMatch[t.largerAttrs[matched][i]];
    }
    for (int i = 0; i < t.smallerAttrs[matched].size(); ++i) {
        if (partMatch[t.smallerAttrs[matched][i]] < minCompared[0])
            minCompared[0] = partMatch[t.smallerAttrs[matched][i]];
    }
    endPoses[0] = traversedLevels[firstLevelNID].back()->length;
    if (maxCompared[0] != 0) {
        beginPoses[0] = setBeginPos(traversedLevels[firstLevelNID].back()->values, traversedLevels[firstLevelNID].back()->length, maxCompared[0]);
    }
    if (minCompared[0] != (VertexID) -  1) {
        endPoses[0] = setEndPos(traversedLevels[firstLevelNID].back()->values, traversedLevels[firstLevelNID].back()->length, minCompared[0]);
    }
    counts[0] = endPoses[0] - beginPoses[0];
    poses[0] = beginPoses[0];
    int depth = 0;
    while (depth >= 0) {
        while (poses[depth] < counts[depth]) {
            VertexID nID = t.nIDs[mappingSize + depth][0];
            VertexID v = traversedLevels[nID].back()->values[poses[depth]];
            ui childPos = poses[depth];
            ++poses[depth];
            if (visited[v]) continue;
            visited[v] = true;
            partMatch[order[mappingSize + depth]] = v;
            if (traversedLevels[nID].back()->children)
                traversedLevels[nID].push_back(&(traversedLevels[nID].back()->children[childPos]));
            else traversedLevels[nID].push_back(nullptr);
            if (depth == totalLevel - 1) {
                produceResult(result, partMatch, t.numAttributes);
                traversedLevels[nID].pop_back();
                visited[v] = false;
            }
            else {
                ++depth;
                beginPoses[depth] = 0;
                nID = t.nIDs[mappingSize + depth][0];
                endPoses[depth] = traversedLevels[nID].back()->length;
                maxCompared[depth] = 0;
                minCompared[depth] = (VertexID) - 1;
                for (int i = 0; i < t.largerAttrs[depth + matched].size(); ++i) {
                    if (partMatch[t.largerAttrs[depth + matched][i]] > maxCompared[depth])
                        maxCompared[depth] = partMatch[t.largerAttrs[depth + matched][i]];
                }
                for (int i = 0; i < t.smallerAttrs[depth + matched].size(); ++i) {
                    if (partMatch[t.smallerAttrs[depth + matched][i]] < minCompared[depth])
                        minCompared[depth] = partMatch[t.smallerAttrs[depth + matched][i]];
                }
                if (maxCompared[depth] != 0) {
                    beginPoses[depth] = setBeginPos(traversedLevels[nID].back()->values, endPoses[nID], maxCompared[depth]);
                }
                if (minCompared[depth] != (VertexID) - 1) {
                    endPoses[depth] = setEndPos(traversedLevels[nID].back()->values, endPoses[depth], minCompared[depth]);
                }
                counts[depth] = endPoses[depth] - beginPoses[depth];
                poses[depth] = beginPoses[depth];
            }
        }
        --depth;
        if (depth >= 0) {
            traversedLevels[t.nIDs[mappingSize + depth][0]].pop_back();
            visited[partMatch[t.globalOrder[mappingSize + depth]]] = false;
        }
    }
}


void
traverse(size_t &count, const Graph &query, const HyperTree &t, const std::vector<VertexID> &order, VertexID *partMatch,
         int mappingSize, const std::vector<TrieLevel *> &levels, std::vector<std::vector<TrieLevel *>> &traversedLevels,
         std::vector<ui> &poses, std::vector<ui> &numBranches, bool *visited, int extendLevel, bool labeled) {
    const std::vector<VertexID> &depthToNID = t.depthToNID;
    const std::vector<ui> &numLevels = t.levels;
    const std::vector<std::vector<VertexID>> &groups = t.groups;
    if (t.defaultPartition.size() > mappingSize) mappingSize = t.defaultPartition.size();
    ui totalLevel = 0;
    for (VertexID nID = 0; nID < levels.size(); ++nID) {
        totalLevel += numLevels[nID];
    }
    if (totalLevel == 0) {
        size_t num = 1;
        if (!t.symmLastLevel.empty()) {
            VertexID nID2 = t.nIDs.back()[0];
            num = levels[nID2] -> numResults(partMatch, t.attributesToCheck);
            num = choosec(num, t.symmLastLevel.size());
        }
        else if (!t.subsetLastLevel.empty()) {
            for (int i = 0; i < t.subsetLastLevel.size(); ++i) {
                VertexID nID2 = t.subsetLastLevel[i];
                size_t num2 = levels[nID2] -> numResults(partMatch, t.subsetToCheck[i], t.largerAttrs[mappingSize + i], t.smallerAttrs[mappingSize + i]);
                if (num2 <= i) {
                    num = 0;
                    break;
                }
                num *= (num2 - i);
            }
        }
        for (auto &group : groups) {
            if (group.size() == 1) {
                VertexID nID2 = group[0];
                if (!labeled) num *= levels[nID2] -> numResults(partMatch, t.attributesToCheck, t.largerAttrs.back(), t.smallerAttrs.back());
                else num *= levels[nID2] -> numTuples(visited);
            }
            else {
                std::vector<std::vector<VertexID>> vsets(group.size());
                for (int i = 0; i < group.size(); ++i) {
                    VertexID nID2 = group[i];
                    vsets[i].reserve(levels[nID2]->length);
                    for (int j = 0; j < levels[nID2]->length; ++j) {
                        VertexID v2 = levels[nID2]->values[j];
                        if (!visited[v2]) vsets[i].push_back(v2);
                    }
                }

                num *= numUniqueCombine(vsets);
            }
        }

        count += num;
    }
    else {
        for (int i = 0; i < totalLevel; ++i) {
            poses[i] = 0;
            numBranches[i] = 0;
        }
        VertexID firstLevelNID = depthToNID[0];
        numBranches[0] = levels[firstLevelNID]->length;
        if (t.hasLargerAttrs[extendLevel]) {
            ui maxCompared = 0;
            for (int i = 0; i < t.largerAttrs[extendLevel].size(); ++i) {
                if (partMatch[t.largerAttrs[extendLevel][i]] > maxCompared)
                    maxCompared = partMatch[t.largerAttrs[extendLevel][i]];
            }
            poses[0] = setBeginPos(levels[firstLevelNID]->values, levels[firstLevelNID]->length, maxCompared);
        }
        if (t.hasSmallerAttrs[extendLevel]) {
            ui minCompared = (VertexID) - 1;
            for (int i = 0; i < t.smallerAttrs[extendLevel].size(); ++i) {
                if (partMatch[t.smallerAttrs[extendLevel][i]] < minCompared)
                    minCompared = partMatch[t.smallerAttrs[extendLevel][i]];
            }
            numBranches[0] = setEndPos(levels[firstLevelNID]->values, levels[firstLevelNID]->length, minCompared);
        }
        int depth = 0;
        while (depth >= 0) {
            const TrieLevel *level = traversedLevels[depthToNID[depth]].back();
            VertexID *values = level->values;
            while (poses[depth] < numBranches[depth]) {
                VertexID nID = depthToNID[depth];
                VertexID v = values[poses[depth]];
                ui childPos = poses[depth];
                ++poses[depth];
                if (!visited[v]) {
                    visited[v] = true;
                    partMatch[order[mappingSize + depth]] = v;
                    if (level->children) {
                        traversedLevels[nID].push_back(&(level->children[childPos]));
                    }
                    else
                        traversedLevels[nID].push_back(nullptr);
                    if (depth + 1 < totalLevel) {
                        ++depth;
                        poses[depth] = 0;
                        nID = depthToNID[depth];
                        level = traversedLevels[nID].back();
                        values = level->values;
                        numBranches[depth] = level->length;
                        if (t.hasLargerAttrs[extendLevel + depth]) {
                            ui maxCompared = 0;
                            for (int i = 0; i < t.largerAttrs[extendLevel + depth].size(); ++i) {
                                if (partMatch[t.largerAttrs[extendLevel + depth][i]] > maxCompared)
                                    maxCompared = partMatch[t.largerAttrs[extendLevel + depth][i]];
                            }
                            poses[depth] = setBeginPos(values, level->length, maxCompared);
                        }
                        if (t.hasSmallerAttrs[extendLevel + depth]) {
                            ui minCompared = (VertexID) - 1;
                            for (int i = 0; i < t.smallerAttrs[extendLevel + depth].size(); ++i) {
                                if (partMatch[t.smallerAttrs[extendLevel + depth][i]] < minCompared)
                                    minCompared = partMatch[t.smallerAttrs[extendLevel + depth][i]];
                            }
                            numBranches[depth] = setEndPos(values, level->length, minCompared);
                        }
                    }
                    else {
                        size_t num = 1;
                        if (!t.symmLastLevel.empty()) {
                            VertexID nID2 = t.nIDs.back()[0];
                            num = traversedLevels[nID2].back() -> numResults(partMatch, t.attributesToCheck);
                            num = choosec(num, t.symmLastLevel.size());
                        }
                        else if (!t.subsetLastLevel.empty()) {
                            for (int i = 0; i < t.subsetLastLevel.size(); ++i) {
                                VertexID nID2 = t.subsetLastLevel[i];
                                size_t num2 = traversedLevels[nID2].back() -> numResults(partMatch, t.subsetToCheck[i], t.largerAttrs[depth + mappingSize + i], t.smallerAttrs[depth + mappingSize + i]);
                                if (num2 <= i) {
                                    num = 0;
                                    break;
                                }
                                num *= (num2 - i);
                            }
                        }
                        for (auto &group : groups) {
                            if (group.size() == 1) {
                                VertexID nID2 = group[0];
                                if (!labeled) num *= traversedLevels[nID2].back() -> numResults(partMatch, t.attributesToCheck, t.largerAttrs.back(), t.smallerAttrs.back());
                                else num *= traversedLevels[nID2].back() -> numTuples(visited);
                            }
                            else {
                                std::vector<std::vector<VertexID>> vsets(group.size());
                                for (int i = 0; i < group.size(); ++i) {
                                    VertexID nID2 = group[i];
                                    vsets[i].reserve(traversedLevels[nID2].back()->length);
                                    for (int j = 0; j < traversedLevels[nID2].back()->length; ++j) {
                                        VertexID v2 = traversedLevels[nID2].back()->values[j];
                                        if (!visited[v2]) vsets[i].push_back(v2);
                                    }
                                }

                                num *= numUniqueCombine(vsets);
                            }
                        }

                        count += num;
                        traversedLevels[nID].pop_back();
                        visited[v] = false;
                    }
                }
            }
            --depth;
            if (depth >= 0) {
                VertexID nID = depthToNID[depth];
                VertexID v = partMatch[order[mappingSize + depth]];
                traversedLevels[nID].pop_back();
                visited[v] = false;
            }
        }
    }
}

void
produceResult(std::vector<std::vector<VertexID>> &result, VertexID *partMatch, ui size) {
#ifdef COLLECT_RESULT
    result.emplace_back(partMatch, partMatch + size);
#endif
#ifdef COLLECT_STATISTICS
    ++gNumResult;
#endif
}

void storeMatches(const std::vector<std::vector<VertexID>> &result, std::ofstream &outStream) {
    outStream << result.size() << std::endl;
    for (int i = 0; i < result.size(); ++i) {
        for (int j = 0; j < result[i].size(); ++j) {
            outStream << result[i][j] << " ";
        }
        outStream << std::endl;
    }
}

void buildTrie(std::vector<std::vector<VertexID>> &tuples, TrieLevel &root, ui &length) {
    std::sort(tuples.begin(), tuples.begin() + length);
    root.buildTrieFromSortedMatchesBatch(tuples, 0, length, 0);
    length = 0;
}
