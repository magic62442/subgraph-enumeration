//
// Created by Qiyan LI on 2024/2/27.
//

#include "decomposition.h"

void HyperNode::initPoses(const Graph &query, const std::vector<std::vector<size_t>> &dist) {
    attributesBefore = std::vector<std::vector<VertexID>>(numAttributes);
    cartesianParent = std::vector<VertexID>(numAttributes);
    bool cartesian = false;
    for (int i = 0; i < numAttributes; ++i) {
        VertexID u = attributes[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = attributes[j];
            if (query.getEdgeID(u, u2) != -1)
                attributesBefore[i].push_back(u2);
        }
        if (i > 0 && attributesBefore[i].empty()) cartesian = true;
    }
    if (cartesian) {
        for (ui i = 1; i < numAttributes; ++i) {
            VertexID u = attributes[i];
            if (attributesBefore[i].empty()) {
                size_t minDist = std::numeric_limits<size_t>::max();
                for (VertexID j = 0; j < i; ++j) {
                    VertexID u2 = attributes[j];
                    if (dist[u][u2] < minDist) {
                        minDist = dist[u][u2];
                        cartesianParent[i] = u2;
                    }
                }
            }
        }
    }
}

void HyperNode::initPoses(const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query,
                          const std::vector<std::vector<size_t>> &dist, bool global) {
    attributesBefore = std::vector<std::vector<VertexID>>(numAttributes);
    cartesianParent = std::vector<VertexID>(numAttributes, 99);
    bool cartesian = false;
    for (int i = 0; i < numAttributes; ++i) {
        VertexID u = attributes[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = attributes[j];
            if (query.getEdgeID(u, u2) != -1)
                attributesBefore[i].push_back(u2);
        }
        if (i > 0 && attributesBefore[i].empty()) cartesian = true;
    }
    if (global) {
        nIDs = std::vector<std::vector<VertexID>>(numAttributes);
        // remove attributes that are covered by materialized bags
        for (int i = 0; i < numAttributes; ++i) {
            std::set<VertexID> coveredAttrs;
            VertexID u = attributes[i];
            for (VertexID nID = 0; nID < sharedAttrs.size() - 1; ++nID) {
                if (std::find(sharedAttrs[nID].begin(), sharedAttrs[nID].end(), u) != sharedAttrs[nID].end()) {
                    nIDs[i].push_back(nID);
                    for (VertexID u2: sharedAttrs[nID])
                        coveredAttrs.insert(u2);
                }
            }
            std::vector<ui> attrsToJoin;
            for (VertexID u2: attributesBefore[i]) {
                if (coveredAttrs.find(u2) == coveredAttrs.end())
                    attrsToJoin.push_back(u2);
            }
            attributesBefore[i] = attrsToJoin;
        }
    }
    if (cartesian) {
        for (ui i = 1; i < numAttributes; ++i) {
            VertexID u = attributes[i];
            if (attributesBefore[i].empty()) {
                size_t minDist = std::numeric_limits<size_t>::max();
                for (VertexID j = 0; j < i; ++j) {
                    VertexID u2 = attributes[j];
                    if (dist[u][u2] < minDist) {
                        minDist = dist[u][u2];
                        cartesianParent[i] = u2;
                    }
                }
            }
        }
    }
}

void HyperNode::copyTo(HyperNode &other) const {
    delete[] other.attributes;
    delete[] other.prefix;
    other.numAttributes = numAttributes;
    other.prefixSize = prefixSize;
    other.fw = fw;
    if (numAttributes != 0) {
        other.attributes = new VertexID[numAttributes];
        std::copy(attributes, attributes + numAttributes, other.attributes);
    } else {
        other.attributes = nullptr;
    }
    if (prefixSize != 0) {
        other.prefix = new VertexID[numAttributes];
        std::copy(prefix, prefix + prefixSize, other.prefix);
    } else {
        other.prefix = nullptr;
    }
    other.attributesBefore = attributesBefore;
    other.cartesianParent = cartesianParent;
    other.nIDs = nIDs;
    other.id = id;
    other.canonValue = canonValue;
}

bool HyperNode::hasVertex(VertexID u) const {
    for (ui i = 0; i < numAttributes; i++)
        if (attributes[i] == u) return true;
    return false;
}

bool HyperNode::isTriangle(const Graph &query) const {
    if (numAttributes != 3) return false;
    VertexID u0 = attributes[0], u1 = attributes[1], u2 = attributes[2];
    if (query.getEdgeID(u0, u1) != -1 && query.getEdgeID(u1, u2) != -1 && query.getEdgeID(u0, u2) != -1) return true;
    else return false;
}

HyperNode::HyperNode(CanonType id, VertexID *attributes, ui numAttributes) : id(id), attributes(attributes),
                                                                             numAttributes(numAttributes) {
    canonValue = 0;
    prefix = nullptr;
    prefixSize = 0;
    fw = 0.0;
}

HyperNode::HyperNode(const HyperNode &other) {
    id = other.id;
    canonValue = other.canonValue;
    numAttributes = other.numAttributes;
    attributes = new VertexID[numAttributes];
    memcpy(attributes, other.attributes, numAttributes * sizeof(VertexID));
    attributesBefore = other.attributesBefore;
    cartesianParent = other.cartesianParent;
    nIDs = other.nIDs;
    prefixSize = other.prefixSize;
    if (prefixSize != 0) {
        prefix = new VertexID[prefixSize];
        memcpy(prefix, other.prefix, prefixSize * sizeof(VertexID));
    }
    else prefix = nullptr;
    fw = other.fw;
    largerAttrs = other.largerAttrs;
    smallerAttrs = other.smallerAttrs;
    hasLargerAttrs = other.hasLargerAttrs;
    hasSmallerAttrs = other.hasSmallerAttrs;
    candidatesBefore = other.candidatesBefore;
    attributesAfter = other.attributesAfter;
    attrAfterLarger = other.attrAfterLarger;
    attrAfterSmaller = other.attrAfterSmaller;
    copyAfter = other.copyAfter;
    copyAfterTypes = other.copyAfterTypes;
    candidatesAfter = other.candidatesAfter;
}

HyperNode &HyperNode::operator=(const HyperNode &other) {
    if (this == &other) return *this;
    id = other.id;
    canonValue = other.canonValue;
    delete[] attributes;
    numAttributes = other.numAttributes;
    if (numAttributes > 0 && other.attributes != nullptr) {
        attributes = new VertexID[numAttributes];
        std::copy(other.attributes, other.attributes + numAttributes, attributes);
    } else {
        attributes = nullptr;
    }
    delete[] prefix;
    prefixSize = other.prefixSize;
    if (prefixSize > 0 && other.prefix != nullptr) {
        prefix = new VertexID[prefixSize];
        std::copy(other.prefix, other.prefix + prefixSize, prefix);
    } else {
        prefix = nullptr;
    }
    attributesBefore = other.attributesBefore;
    cartesianParent = other.cartesianParent;
    nIDs = other.nIDs;
    fw = other.fw;
    largerAttrs = other.largerAttrs;
    smallerAttrs = other.smallerAttrs;
    hasLargerAttrs = other.hasLargerAttrs;
    hasSmallerAttrs = other.hasSmallerAttrs;
    candidatesBefore = other.candidatesBefore;
    attributesAfter = other.attributesAfter;
    attrAfterLarger = other.attrAfterLarger;
    attrAfterSmaller = other.attrAfterSmaller;
    copyAfter = other.copyAfter;
    copyAfterTypes = other.copyAfterTypes;
    candidatesAfter = other.candidatesAfter;
    return *this;
}

HyperTree::HyperTree(const HyperTree &other) {
    numAttributes = other.numAttributes;
    numNodes = other.numNodes;
    nodes = new HyperNode[numNodes];
    for (ui i = 0; i < numNodes; ++i)
        nodes[i] = other.nodes[i];
    v2n = new std::vector<VertexID>[numAttributes];
    for (ui i = 0; i < numAttributes; ++i)
        v2n[i] = other.v2n[i];
    edges = other.edges;
    compressionSizes = other.compressionSizes;
    nIDs = other.nIDs;
    extendLevel = other.extendLevel;
    globalOrder = other.globalOrder;
    trieOrder = other.trieOrder;
    nodesAtStep = other.nodesAtStep;
    defaultPartition = other.defaultPartition;
    attributesBefore = other.attributesBefore;
    newGlobalNode = other.newGlobalNode;
    cartesianParent = other.cartesianParent;
    depthToNID = other.depthToNID;
    trieParents = other.trieParents;
    levels = other.levels;
    groups = other.groups;
    adaptiveDepthToNID = other.adaptiveDepthToNID;
    adaptiveLevels = other.adaptiveLevels;
    adaptiveGroups = other.adaptiveGroups;
    largerAttrs = other.largerAttrs;
    smallerAttrs = other.smallerAttrs;
    hasLargerAttrs = other.hasLargerAttrs;
    hasSmallerAttrs = other.hasSmallerAttrs;
    symmetryRules = other.symmetryRules;
    symmLastLevel = other.symmLastLevel;
    subsetLastLevel = other.subsetLastLevel;
    subsetToCheck = other.subsetToCheck;
    candIndex = other.candIndex;
    divideFactor = other.divideFactor;
    attributesToCheck = other.attributesToCheck;
}

HyperTree& HyperTree::operator=(const HyperTree& other) {
    if (this == &other) {
        return *this;
    }

    delete[] nodes;
    delete[] v2n;
    numAttributes = other.numAttributes;
    numNodes = other.numNodes;
    if (numNodes > 0) {
        nodes = new HyperNode[numNodes];
        for (ui i = 0; i < numNodes; ++i) {
            nodes[i] = other.nodes[i];
        }
    } else {
        nodes = nullptr;
    }
    if (numAttributes > 0) {
        v2n = new std::vector<VertexID>[numAttributes];
        for (ui i = 0; i < numAttributes; ++i) {
            v2n[i] = other.v2n[i];
        }
    } else {
        v2n = nullptr;
    }
    edges = other.edges;
    compressionSizes = other.compressionSizes;
    nIDs = other.nIDs;
    extendLevel = other.extendLevel;
    globalOrder = other.globalOrder;
    trieOrder = other.trieOrder;
    nodesAtStep = other.nodesAtStep;
    defaultPartition = other.defaultPartition;
    attributesBefore = other.attributesBefore;
    newGlobalNode = other.newGlobalNode;
    cartesianParent = other.cartesianParent;
    depthToNID = other.depthToNID;
    trieParents = other.trieParents;
    levels = other.levels;
    groups = other.groups;
    adaptiveDepthToNID = other.adaptiveDepthToNID;
    adaptiveLevels = other.adaptiveLevels;
    adaptiveGroups = other.adaptiveGroups;
    largerAttrs = other.largerAttrs;
    smallerAttrs = other.smallerAttrs;
    hasLargerAttrs = other.hasLargerAttrs;
    hasSmallerAttrs = other.hasSmallerAttrs;
    symmetryRules = other.symmetryRules;
    symmLastLevel = other.symmLastLevel;
    subsetLastLevel = other.subsetLastLevel;
    subsetToCheck = other.subsetToCheck;
    candIndex = other.candIndex;
    divideFactor = other.divideFactor;
    attributesToCheck = other.attributesToCheck;

    return *this;
}


void HyperTree::initPoses(const Graph &query, const CandidateSpace &cs, bool handleNode) {
    numAttributes = query.getNumVertices();
    delete[] v2n;
    v2n = new std::vector<VertexID> [numAttributes];
    for (ui i = 0; i < numNodes; ++i) {
        for (ui j = 0; j < nodes[i].numAttributes; ++j) {
            v2n[nodes[i].attributes[j]].push_back(i);
        }
    }
    if (handleNode) {
        for (int i = 0; i < numNodes; ++i) {
            nodes[i].initPoses(query, cs.dist);
        }
    }

    // set the compression sizes
    compressionSizes = std::vector<ui>(numNodes, 0);
//    extendLevel = 0;
//    for (int i = 0; i < globalOrder.size(); ++i) {
//        if (v2n[globalOrder[i]].size() > 1)
//            extendLevel = i + 1;
//    }
//    if (!newGlobalNode && extendLevel < nodes[numNodes - 1].prefixSize) extendLevel = nodes[numNodes - 1].prefixSize;
    buildTrieOrder();
    nodesAtStep = std::vector<std::vector<VertexID>>(defaultPartition.size());
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < defaultPartition.size(); ++i) {
            if (std::includes(defaultPartition.begin(), defaultPartition.begin() + i + 1,
                              nodes[nID].prefix, nodes[nID].prefix + nodes[nID].prefixSize)) {
                nodesAtStep[i].push_back(nID);
                break;
            }
        }
    }
    nIDs = std::vector<std::vector<VertexID>>(globalOrder.size());
    for (int i = 0; i < globalOrder.size(); ++i) {
        VertexID u = globalOrder[i];
        nIDs[i] = v2n[u];
    }
    attributesBefore = std::vector<std::vector<VertexID>>(globalOrder.size());
    cartesianParent = std::vector<VertexID>(globalOrder.size());
    for (int i = 0; i < globalOrder.size(); ++i) {
        VertexID u = globalOrder[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = globalOrder[j];
            if (query.getEdgeID(u, u2) != -1)
                attributesBefore[i].push_back(u2);
        }
    }
    for (ui i = 1; i < globalOrder.size(); ++i) {
        VertexID u = globalOrder[i];
        if (attributesBefore[i].empty()) {
            size_t minDist = std::numeric_limits<size_t>::max();
            for (VertexID j = 0; j < i; ++j) {
                VertexID u2 = globalOrder[j];
                if (cs.dist[u][u2] < minDist) {
                    minDist = cs.dist[u][u2];
                    cartesianParent[i] = u2;
                }
            }
        }
    }
    adaptiveDepthToNID.clear();
    adaptiveGroups.clear();
    adaptiveLevels.clear();
//    buildTraverseStruct(query);
}

void HyperTree::buildFromTD(FHD &fhd) {
    numNodes = fhd.X.size();
    nodes = new HyperNode[numNodes + 1];
    for (VertexID nID = 0; nID < fhd.X.size(); ++nID) {
        std::vector<size_t> tmpV;
        fhd.X[nID].getelement(tmpV);
        nodes[nID].numAttributes = tmpV.size();
        nodes[nID].attributes = new VertexID[nodes[nID].numAttributes];
        for (int i = 0; i < tmpV.size(); ++i) {
            nodes[nID].attributes[i] = tmpV[i];
        }
    }
    edges = std::vector<std::vector<VertexID>>(numNodes);
    for (auto ei: fhd.eg)
        edges[ei.first].push_back(ei.second);
}

void HyperTree::print(const Graph &query) const {
    for (int i = 0; i < numNodes; ++i) {
        std::cout << "node " << i << " : ";
        for (int j = 0; j < nodes[i].numAttributes; ++j) {
            std::cout << nodes[i].attributes[j] << " ";
        }
        std::cout << std::endl;
        for (int j = 0; j < nodes[i].numAttributes; ++j) {
            for (int k = j + 1; k < nodes[i].numAttributes; ++k) {
                if (query.getEdgeID(nodes[i].attributes[j], nodes[i].attributes[k]) != -1) {
                    std::cout << nodes[i].attributes[j] << " " << nodes[i].attributes[k] << std::endl;
                }
            }
        }
        std::cout << std::endl;
    }
}

void HyperTree::addGlobalNode(const std::vector<VertexID> &globalAttrs) {
    HyperNode *old = nodes;
    nodes = new HyperNode[numNodes + 1];
    for (int i = 0; i < numNodes; ++i) {
        old[i].copyTo(nodes[i]);
    }
    delete[] old;
    nodes[numNodes] = HyperNode();
    nodes[numNodes].attributes = new VertexID[globalAttrs.size()];
    nodes[numNodes].numAttributes = globalAttrs.size();
    for (int i = 0; i < globalAttrs.size(); ++i) {
        VertexID u = globalAttrs[i];
        nodes[numNodes].attributes[i] = u;
        if (v2n[u].back() != numNodes) v2n[u].push_back(numNodes);
    }

    ++numNodes;
    newGlobalNode = true;
}

void HyperTree::buildTrieOrder() {
    trieOrder = std::vector<std::vector<VertexID>>(numNodes);
    for (ui i = 0; i < numNodes - 1; ++i) {
        for (ui j = 0; j < globalOrder.size(); ++j) {
            VertexID u = globalOrder[j];
            ui k = nodes[i].prefixSize;
            for (; k < nodes[i].numAttributes; ++k) {
                if (nodes[i].attributes[k] == u) {
                    break;
                }
            }
            if (k != nodes[i].numAttributes) trieOrder[i].push_back(k - nodes[i].prefixSize);
        }
    }
    for (ui j = 0; j < globalOrder.size(); ++j) {
        VertexID u = globalOrder[j];
        ui k = extendLevel;
        for (; k < nodes[numNodes - 1].numAttributes; ++k) {
            if (nodes[numNodes - 1].attributes[k] == u) {
                break;
            }
        }
        if (k != nodes[numNodes - 1].numAttributes) trieOrder[numNodes - 1].push_back(k - extendLevel);
    }
}

void HyperTree::buildTraverseStruct(const Graph &query) {
    buildTraverseAdaptive(query, extendLevel);
    depthToNID = adaptiveDepthToNID[extendLevel];
    levels = adaptiveLevels[extendLevel];
    groups = adaptiveGroups[extendLevel];
    if (!newGlobalNode && extendLevel > nodes[numNodes - 1].prefixSize) {
        trieOrder[numNodes - 1].clear();
        for (int i = 0; i < globalOrder.size(); ++i) {
            VertexID u = globalOrder[i];
            ui k = extendLevel;
            for (; k < nodes[numNodes - 1].numAttributes; ++k) {
                if (nodes[numNodes - 1].attributes[k] == u)
                    break;
            }
            if (k != nodes[numNodes - 1].numAttributes) trieOrder[numNodes - 1].push_back(k - extendLevel);
        }
    }
}

void HyperTree::buildTraverseAdaptive(const Graph &query, int currentExtendLevel) {
    if (adaptiveDepthToNID.find(currentExtendLevel) != adaptiveDepthToNID.end()) return;
    std::vector<VertexID> toExtend;
    int pos = currentExtendLevel;
    if (defaultPartition.size() > nodes[numNodes - 1].numAttributes) pos = defaultPartition.size();
    std::vector<ui> currentLevels = std::vector<ui>(numNodes, 0);
    std::vector<std::vector<VertexID>> localAttr(numNodes);
    for (int i = pos; i < numAttributes; ++i) {
        VertexID u = globalOrder[i];
        VertexID nID = v2n[u][0];
        ++currentLevels[nID];
        localAttr[nID].push_back(u);
    }
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        if (currentLevels[nID] != 0) toExtend.push_back(nID);
    }
    std::sort(toExtend.begin(), toExtend.end(), [currentLevels](VertexID a, VertexID b) {
        return currentLevels[a] > currentLevels[b];
    });
    // divide bags into groups based on labels
    // for each bag, if the last levels have distinct labels, directly multiply the num
    // else, extend to the last level and do intersection for each label
    std::map<LabelID, std::vector<VertexID>> labelToNodes;
    std::set<VertexID> independents;
    // all bags start from one level
    for (VertexID nID : toExtend) {
        VertexID u = localAttr[nID].back();
        LabelID l = query.getVertexLabel(u);
        labelToNodes[l].push_back(nID);
    }
    std::vector<std::vector<VertexID>> currentGroups;
    for (auto &item : labelToNodes) {
        if (item.second.size() > 1) currentGroups.push_back(item.second);
        else independents.insert(item.second[0]);
    }
    for (VertexID nID : toExtend) {
        --currentLevels[nID];
        if (independents.find(nID) == independents.end()) continue;
        ui numAttrs = nodes[nID].numAttributes - nodes[nID].prefixSize;
        VertexID u;
        LabelID l;
        if (currentLevels[nID] > 0) {
            u = localAttr[nID][currentLevels[nID] - 1];
            l = query.getVertexLabel(u);
        }
        std::vector<VertexID> group = {nID};
        currentGroups.push_back(group);
        // expend the level of nID
        while (currentLevels[nID] > 0 && (labelToNodes.find(l) == labelToNodes.end() || labelToNodes[l][0] == nID)) {
            if (labelToNodes.find(l) == labelToNodes.end()) labelToNodes[l].push_back(nID);
            --currentLevels[nID];
            if (currentLevels[nID] > 0) {
                u = localAttr[nID][currentLevels[nID] - 1];
                l = query.getVertexLabel(u);
            }
        }
    }
    std::vector<VertexID> currentdepthToNID;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < currentLevels[nID]; ++i) {
            currentdepthToNID.push_back(nID);
        }
    }
    adaptiveDepthToNID[currentExtendLevel] = currentdepthToNID;
    adaptiveLevels[currentExtendLevel] = currentLevels;
    adaptiveGroups[currentExtendLevel] = currentGroups;
}

void HyperTree::writeToStream(std::ostream &outStream) {
    outStream << "num bags: " << numNodes << std::endl;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < nodes[nID].numAttributes; ++i) {
            outStream << nodes[nID].attributes[i];
            if (i != nodes[nID].numAttributes - 1) outStream << ",";
        }
        outStream << std::endl;
    }
    for (auto rule: symmetryRules) {
        outStream << rule[0] << " < ";
        for (int i = 1; i < rule.size(); ++i)
            outStream << rule[i] << ",";
        outStream << std::endl;
    }
}

//void HyperTree::selectSymmetry(const PatternGraph &p) {
//    std::vector<std::vector<VertexID>> sortedAttributes;
//    for (VertexID nID = 0; nID < numNodes; ++nID) {
//        sortedAttributes.emplace_back(nodes[nID].attributes, nodes[nID].attributes + nodes[nID].numAttributes);
//        std::sort(sortedAttributes.back().begin(), sortedAttributes.back().end());
//    }
//    const std::vector<std::vector<std::vector<VertexID>>> &candidateRules = p.getCandidateRules();
//    // priorities: 1. the size of rules not in bag: smaller better
//    // 2. the last level to have a rule: earlier is better
//    int minNotContained = p.getNumVertices(), minRuleSize = p.getNumVertices(), minPos = p.getNumVertices(), minIdx = 0;
//    for (int i = 0; i < candidateRules.size(); ++i) {
//        const auto &ruleSet = candidateRules[i];
//        int numNotContained = ruleSet.size();
//        int ruleSize = 0, lastPos = 0;
//        std::vector<std::pair<VertexID, VertexID>> pairs;
//        for (int j = 0; j < ruleSet.size(); ++j) {
//            std::vector<VertexID> rule = ruleSet[j];
//            bool flag = true;
//            VertexID u = rule[0];
//            for (int k = 1; k < rule.size(); ++k) {
//                VertexID u2 = rule[k];
//                bool contain = false;
//                for (VertexID nID = 0; nID < numNodes; ++nID) {
//                    if (std::find(sortedAttributes[nID].begin(), sortedAttributes[nID].end(), u) != sortedAttributes[nID].end() &&
//                        std::find(sortedAttributes[nID].begin(), sortedAttributes[nID].end(), u2) != sortedAttributes[nID].end())
//                        contain = true;
//                }
//                if (!contain) {
//                    flag = false;
//                    break;
//                }
//            }
//            if (flag)  --numNotContained;
//            else {
//                if (ruleSize < rule.size())
//                    ruleSize = rule.size();
//                for (int k = 0; k < globalOrder.size(); ++k) {
//                    u = globalOrder[k];
//                    if (std::find(rule.begin(), rule.end(), u) != rule.end()) {
//                        if (lastPos < k) lastPos = k;
//                    }
//                }
//            }
//        }
//        if (numNotContained < minNotContained) {
//            minNotContained = numNotContained;
//            minRuleSize = ruleSize;
//            minPos = lastPos;
//            minIdx = i;
//        }
//        else if (numNotContained == minNotContained) {
//            if (ruleSize < minRuleSize) {
//                minRuleSize = ruleSize;
//                minPos = lastPos;
//                minIdx = i;
//            }
//            else if (ruleSize == minRuleSize) {
//                if (lastPos < minPos) {
//                    minPos = lastPos;
//                    minIdx = i;
//                }
//            }
//        }
//    }
//    symmetryRules = candidateRules[minIdx];
//}

void HyperTree::selectSymmetry(const PatternGraph &p) {
    // Sort all HyperNode ids based on specified criteria
    std::vector<VertexID> sortedNodeIds;

    // Initialize with all node IDs from 0 to numNodes-1
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        sortedNodeIds.push_back(nID);
    }

    // If newGlobalNode is true, put the last node (numNodes-1) at the front
    if (newGlobalNode) {
        // Move the last node to the front
        VertexID lastNodeId = numNodes - 1;
        auto it = std::find(sortedNodeIds.begin(), sortedNodeIds.end(), lastNodeId);
        if (it != sortedNodeIds.end()) {
            sortedNodeIds.erase(it);
            sortedNodeIds.insert(sortedNodeIds.begin(), lastNodeId);
        }
    }

    // Sort the remaining nodes (excluding the first one if newGlobalNode is true)
    size_t startIdx = newGlobalNode ? 1 : 0;
    std::sort(sortedNodeIds.begin() + startIdx, sortedNodeIds.end(),
              [this, &p](VertexID id1, VertexID id2) {
        const HyperNode &node1 = nodes[id1];
        const HyperNode &node2 = nodes[id2];

        // First priority: fw larger first
        if (node1.fw != node2.fw) {
            return node1.fw > node2.fw;
        }

        // If fw is the same, count edges in each node
        auto countEdges = [&p](const HyperNode &node) -> int {
            int edgeCount = 0;
            for (int i = 0; i < node.numAttributes; ++i) {
                for (int j = i + 1; j < node.numAttributes; ++j) {
                    if (p.isEdge(node.attributes[i], node.attributes[j])) {
                        edgeCount++;
                    }
                }
            }
            return edgeCount;
        };

        int edges1 = countEdges(node1);
        int edges2 = countEdges(node2);

        // If fw is the same and not 0, fewer edges come first
        if (node1.fw != 0.0 && node1.fw == node2.fw) {
            if (edges1 != edges2) {
                return edges1 < edges2;
            }
        }
        // Otherwise, more numAttributes come first
        else {
            if (node1.numAttributes != node2.numAttributes) {
                return node1.numAttributes > node2.numAttributes;
            }
        }

        // If all criteria are equal, maintain original order
        return id1 < id2;
    });
    std::vector<std::vector<VertexID>> sortedAttributes;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        sortedAttributes.emplace_back(nodes[nID].attributes, nodes[nID].attributes + nodes[nID].numAttributes);
        std::sort(sortedAttributes.back().begin(), sortedAttributes.back().end());
    }
    const std::vector<std::vector<std::vector<VertexID>>> &candidateRules = p.getCandidateRules();

    // Create group structure for sortedNodeIds based on canonValue
    std::vector<int> startPoses;
    startPoses.push_back(0);

    for (int i = 1; i < sortedNodeIds.size(); ++i) {
        VertexID prevNodeId = sortedNodeIds[i-1];
        VertexID currNodeId = sortedNodeIds[i];

        if (nodes[prevNodeId].canonValue != nodes[currNodeId].canonValue) {
            startPoses.push_back(i);
        }
    }
    startPoses.push_back(numNodes); // Add end position

    // Calculate score for each candidateRule and collect unique scores
    std::map<std::vector<int>, int> scoreToIdx; // map from score to idx

    for (int i = 0; i < candidateRules.size(); ++i) {
        const auto &ruleSet = candidateRules[i];

        // Count how many <u, u2> pairs each HyperNode contains
        std::vector<int> nodePairCounts(numNodes, 0);

        for (int j = 0; j < ruleSet.size(); ++j) {
            const std::vector<VertexID> &rule = ruleSet[j];
            VertexID u = rule[0];
            for (int k = 1; k < rule.size(); ++k) {
                VertexID u2 = rule[k];

                // Check which HyperNodes contain this <u, u2> pair
                for (VertexID nID = 0; nID < numNodes; ++nID) {
                    if (std::find(sortedAttributes[nID].begin(), sortedAttributes[nID].end(), u) != sortedAttributes[nID].end() &&
                        std::find(sortedAttributes[nID].begin(), sortedAttributes[nID].end(), u2) != sortedAttributes[nID].end()) {
                        nodePairCounts[nID]++;
                    }
                }
            }
        }

        // Create score vector based on sortedNodeIds order
        std::vector<int> currentScore(numNodes);
        for (int j = 0; j < numNodes; ++j) {
            VertexID nodeId = sortedNodeIds[j];
            currentScore[j] = nodePairCounts[nodeId];
        }

        // Keep only one idx for each unique score
        if (scoreToIdx.find(currentScore) == scoreToIdx.end()) {
            scoreToIdx[currentScore] = i;
        }
    }

    // Filter by canonValue grouping and find the best score
    std::map<std::vector<std::vector<int>>, int> canonGroupedScoreToIdx;

    for (const auto& entry : scoreToIdx) {
        const std::vector<int>& score = entry.first;
        int idx = entry.second;

        // Group scores by canonValue
        std::map<CanonType, std::vector<int>> canonToPairCounts;

        for (int j = 0; j < score.size(); ++j) {
            VertexID nodeId = sortedNodeIds[j];
            CanonType canonVal = nodes[nodeId].canonValue;
            canonToPairCounts[canonVal].push_back(score[j]);
        }

        // Create 2D vector: for each canonValue, sort pair counts in ascending order
        std::vector<std::vector<int>> canonGroupedScore;
        for (auto& pair : canonToPairCounts) {
            std::vector<int>& counts = pair.second;
            std::sort(counts.begin(), counts.end()); // ascending order
            canonGroupedScore.push_back(counts);
        }

        // Keep only one idx for each unique canonGroupedScore
        if (canonGroupedScoreToIdx.find(canonGroupedScore) == canonGroupedScoreToIdx.end()) {
            canonGroupedScoreToIdx[canonGroupedScore] = idx;
        }
    }

    // Convert to pairs for sorting
    std::vector<std::pair<std::vector<int>, int>> scoreIdxPairs;
    for (const auto& entry : canonGroupedScoreToIdx) {
        // Find the original score for this idx
        int idx = entry.second;
        for (const auto& originalEntry : scoreToIdx) {
            if (originalEntry.second == idx) {
                scoreIdxPairs.push_back({originalEntry.first, idx});
                break;
            }
        }
    }

    // Sort by score using group-wise multiplication comparison
    std::sort(scoreIdxPairs.begin(), scoreIdxPairs.end(),
              [&startPoses](const std::pair<std::vector<int>, int>& a, const std::pair<std::vector<int>, int>& b) {
        const std::vector<int>& scoreA = a.first;
        const std::vector<int>& scoreB = b.first;

        // Compare each group by multiplying pair counts within the group
        for (int groupIdx = 0; groupIdx < startPoses.size() - 1; ++groupIdx) {
            int start = startPoses[groupIdx];
            int end = startPoses[groupIdx + 1];

            long long productA = 1;
            long long productB = 1;

            for (int i = start; i < end; ++i) {
                productA *= scoreA[i];
                productB *= scoreB[i];
            }

            if (productA != productB) {
                return productA > productB; // larger product first
            }
        }
        return false; // equal scores
    });

    // Select the best candidateRule
    int maxIdx = scoreIdxPairs.empty() ? 0 : scoreIdxPairs[0].second;
    symmetryRules = candidateRules[maxIdx];
}

void HyperTree::setSymmetry() {
    std::vector<int> globalAttrToPos;
    for (VertexID u = 0; u < globalOrder.size(); ++u) {
        globalAttrToPos.push_back(std::find(globalOrder.begin(), globalOrder.end(), u) - globalOrder.begin());
        largerAttrs[u].clear();
        smallerAttrs[u].clear();
    }
    std::vector<std::vector<int>> attrToPos;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        attrToPos.emplace_back(nodes[nID].numAttributes);
        for (int i = 0; i < nodes[nID].numAttributes; ++i) {
            attrToPos[nID][i] = globalAttrToPos[nodes[nID].attributes[i]];
        }
    }
    for (const auto &rule: symmetryRules) {
        bool flag = false;
        VertexID u1 = rule[0];
        for (int i = 1; i < rule.size(); ++i) {
            VertexID u2 = rule[i];
            bool flag = false;
            for (VertexID nID = 0; nID < numNodes; ++nID) {
                if (std::find(nodes[nID].attributes, nodes[nID].attributes + nodes[nID].numAttributes, u1) != nodes[nID].attributes + nodes[nID].numAttributes &&
                    std::find(nodes[nID].attributes, nodes[nID].attributes + nodes[nID].numAttributes, u2) != nodes[nID].attributes + nodes[nID].numAttributes) {
                    flag = true;
                    std::vector<int> attrToPos(numAttributes);
                    for (int j = 0; j < nodes[nID].numAttributes; ++j) {
                        attrToPos[nodes[nID].attributes[j]] = j;
                    }
                    if (attrToPos[u1] < attrToPos[u2]) {
                        nodes[nID].largerAttrs[attrToPos[u2]].push_back(u1);
                        nodes[nID].hasLargerAttrs[attrToPos[u2]] = true;
                    }
                    else {
                        nodes[nID].smallerAttrs[attrToPos[u1]].push_back(u2);
                        nodes[nID].hasSmallerAttrs[attrToPos[u1]] = true;
                    }
                }
            }
            if (!flag) {
                if (globalAttrToPos[u1] < globalAttrToPos[u2]) {
                    largerAttrs[globalAttrToPos[u2]].push_back(u1);
                    hasLargerAttrs[globalAttrToPos[u2]] = true;
                }
                else {
                    smallerAttrs[globalAttrToPos[u1]].push_back(u2);
                    hasSmallerAttrs[globalAttrToPos[u1]] = true;
                }
            }
        }
    }
}

void HyperTree::refineSymmetryAttrs() {
    // Build partial order from symmetry rules
    std::map<VertexID, std::set<VertexID>> partialOrder;

    // Initialize empty sets for all vertices that appear in rules
    for (const auto &rule : symmetryRules) {
        VertexID u = rule[0];
        if (partialOrder.find(u) == partialOrder.end()) {
            partialOrder[u] = std::set<VertexID>();
        }
    }

    // Build transitive closure of partial order: u < v means u is smaller than v
    for (const auto &rule : symmetryRules) {
        for (int j = 1; j < rule.size(); ++j) {
            partialOrder[rule[0]].insert(rule[j]);
        }
    }

    // Compute transitive closure
    bool changed = true;
    while (changed) {
        changed = false;
        for (auto &[u, larger] : partialOrder) {
            std::set<VertexID> newLarger = larger;
            for (VertexID v : larger) {
                for (VertexID w : partialOrder[v]) {
                    if (newLarger.find(w) == newLarger.end()) {
                        newLarger.insert(w);
                        changed = true;
                    }
                }
            }
            larger = newLarger;
        }
    }

    // Helper function to get maximal elements
    auto getMaximalElements = [&](const std::vector<VertexID> &vertices) -> std::vector<VertexID> {
        std::vector<VertexID> maximal;
        for (VertexID u : vertices) {
            bool isMaximal = true;
            for (VertexID v : vertices) {
                if (u != v && partialOrder[u].count(v) > 0) {
                    isMaximal = false;
                    break;
                }
            }
            if (isMaximal) {
                maximal.push_back(u);
            }
        }
        return maximal;
    };

    // Helper function to get minimal elements
    auto getMinimalElements = [&](const std::vector<VertexID> &vertices) -> std::vector<VertexID> {
        std::vector<VertexID> minimal;
        for (VertexID u : vertices) {
            bool isMinimal = true;
            for (VertexID v : vertices) {
                if (u != v && partialOrder[v].count(u) > 0) {
                    isMinimal = false;
                    break;
                }
            }
            if (isMinimal) {
                minimal.push_back(u);
            }
        }
        return minimal;
    };

    // Refine global largerAttrs and smallerAttrs
    for (int i = 0; i < largerAttrs.size(); ++i) {
        if (largerAttrs[i].size() > 1) {
            largerAttrs[i] = getMaximalElements(largerAttrs[i]);
        }
    }

    for (int i = 0; i < smallerAttrs.size(); ++i) {
        if (smallerAttrs[i].size() > 1) {
            smallerAttrs[i] = getMinimalElements(smallerAttrs[i]);
        }
    }

    // Refine node-level largerAttrs and smallerAttrs
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < nodes[nID].largerAttrs.size(); ++i) {
            if (nodes[nID].largerAttrs[i].size() > 1) {
                nodes[nID].largerAttrs[i] = getMaximalElements(nodes[nID].largerAttrs[i]);
            }
        }

        for (int i = 0; i < nodes[nID].smallerAttrs.size(); ++i) {
            if (nodes[nID].smallerAttrs[i].size() > 1) {
                nodes[nID].smallerAttrs[i] = getMinimalElements(nodes[nID].smallerAttrs[i]);
            }
        }
        for (int i = 0; i < nodes[nID].attrAfterLarger.size(); ++i) {
            for (int j = 0; j < nodes[nID].attrAfterLarger[i].size(); ++j) {
                if (nodes[nID].attrAfterLarger[i][j].size() > 1)
                    nodes[nID].attrAfterLarger[i][j] = getMaximalElements(nodes[nID].attrAfterLarger[i][j]);
            }
            for (int j = 0; j < nodes[nID].attrAfterSmaller[i].size(); ++j) {
                if (nodes[nID].attrAfterSmaller[i][j].size() > 1)
                    nodes[nID].attrAfterSmaller[i][j] = getMinimalElements(nodes[nID].attrAfterSmaller[i][j]);
            }
        }
    }
}

void HyperTree::buildTraverseUnlabeled(const Graph &query, CandidateSpace &cs, bool iep) {
    if (numNodes < 2) {
        iep = false;
        return;
    }
    int pos = extendLevel;
    if (defaultPartition.size() > nodes[numNodes - 1].numAttributes) pos = defaultPartition.size();
    std::vector<ui> currentLevels = std::vector<ui>(numNodes, 0);
    std::vector<std::vector<VertexID>> currentGroups;
    std::vector<VertexID> currentDepthToNID;
    std::vector<std::vector<VertexID>> globalSymmetry;
    int pos1 = numAttributes - 1;
    while (v2n[globalOrder[pos1]].size() == 1) --pos1;
    std::vector<VertexID> newGlobalOrder(globalOrder.begin(), globalOrder.begin() + pos1 + 1);
    std::vector<VertexID> notTraversed;
    std::vector<int> numIncluded(numNodes, 0);
    for (VertexID u: newGlobalOrder) {
        for (VertexID nID: v2n[u]) ++numIncluded[nID];
    }
    for (int i = pos1 + 1; i < numAttributes; ++i) {
        VertexID u = globalOrder[i];
        VertexID nID = v2n[u][0];
        if (newGlobalNode && nID == numNodes - 1) continue;
        if (numIncluded[nID] < nodes[nID].numAttributes - 1) {
            ++numIncluded[nID];
            newGlobalOrder.push_back(u);
        }
        else notTraversed.push_back(u);
    }
    for (const auto &rule: symmetryRules) {
        if (unorderedSubset(rule, notTraversed))
            globalSymmetry.push_back(rule);
    }
    std::vector<std::vector<VertexID>> vertexParents(notTraversed.size());
    for (int i = 0; i < notTraversed.size(); ++i) {
        VertexID u = notTraversed[i];
        VertexID nID = v2n[u][0];
        for (int j = 0; j < nodes[nID].numAttributes; ++j) {
            if (query.getEdgeID(u, nodes[nID].attributes[j]) != -1)
                vertexParents[i].push_back(nodes[nID].attributes[j]);
        }
    }
    if (iep && !globalSymmetry.empty()) {
        int symPos = 0;
        ui maxSize = 0;
        for (int i = 0; i < globalSymmetry.size(); ++i) {
            ui size = globalSymmetry[i].size();
            if (size > maxSize) {
                maxSize = size;
                symPos = i;
            }
        }
        const std::vector<VertexID> &rule = globalSymmetry[symPos];
        bool flag = true;
        for (VertexID u: rule) {
            VertexID nID = v2n[u][0];
            if (u != nodes[nID].attributes[nodes[nID].numAttributes - 1])
                flag = false;
        }
        if (flag) {
            for (int i = 0; i < rule.size(); ++i) {
                VertexID nID = v2n[rule[i]][0];
                symmLastLevel.push_back(nID);
            }
            for (VertexID u: notTraversed) {
                if (std::find(rule.begin(), rule.end(), u) == rule.end())
                    newGlobalOrder.push_back(u);
            }
            for (VertexID u: rule) newGlobalOrder.push_back(u);
        }
    }
    // the largest subset
    else if (iep){
        std::vector<std::pair<int, int>> subsetRelation;
        for (int i = 0; i < notTraversed.size(); ++i) {
            VertexID u1 = notTraversed[i];
            VertexID nID1 = v2n[u1][0];
            for (int j = i + 1; j < notTraversed.size(); ++j) {
                VertexID u2 = notTraversed[j];
                VertexID nID2 = v2n[u2][0];
                if (unorderedSubset(vertexParents[i], vertexParents[j]) && unorderedSubset(vertexParents[j], vertexParents[i])) continue;
                if (unorderedSubset(vertexParents[i], vertexParents[j]) && unorderedSubset(nodes[nID1].largerAttrs.back(),
                   nodes[nID2].largerAttrs.back()) && unorderedSubset(nodes[nID1].smallerAttrs.back(), nodes[nID2].smallerAttrs.back())) {
                    subsetRelation.emplace_back(j, i);
                }
                if (unorderedSubset(vertexParents[j], vertexParents[i]) && unorderedSubset(nodes[nID2].largerAttrs.back(),
                   nodes[nID1].largerAttrs.back()) && unorderedSubset(nodes[nID2].smallerAttrs.back(), nodes[nID1].smallerAttrs.back()))
                    subsetRelation.emplace_back(i, j);
            }
        }
        if (!subsetRelation.empty()) {
            std::vector<VertexID> lastVertices;
            std::vector<int> chain = buildMaxSubsetChain(subsetRelation, numAttributes);
            for (int i: chain) {
                lastVertices.push_back(notTraversed[i]);
                subsetLastLevel.push_back(v2n[notTraversed[i]][0]);
            }
            for (VertexID u: notTraversed) {
                if (std::find(lastVertices.begin(), lastVertices.end(), u) == lastVertices.end())
                    newGlobalOrder.push_back(u);
            }
            for (VertexID u: lastVertices) newGlobalOrder.push_back(u);
        }
    }
    if (symmLastLevel.empty() && subsetLastLevel.empty()) {
        // the vertex with the largest card
        VertexID lastVertex = notTraversed[0];
        ui minNumBackWard = numAttributes, minNumRule = numAttributes;
        for (int i = 0; i < notTraversed.size(); ++i) {
            VertexID u = notTraversed[i];
            VertexID nID = v2n[u][0];
            ui numRule = nodes[nID].largerAttrs.back().size() + nodes[nID].smallerAttrs.back().size();
            if (vertexParents[i].size() < minNumBackWard || vertexParents[i].size() == minNumBackWard &&
                numRule < minNumRule) {
                lastVertex = u;
                minNumBackWard = vertexParents[i].size();
                minNumRule = numRule;
            }
        }
        for (VertexID u: notTraversed)
            if (u != lastVertex) newGlobalOrder.push_back(u);
        newGlobalOrder.push_back(lastVertex);
    }
    globalOrder = newGlobalOrder;
    initPoses(query, cs, false);
    setSymmetry();
    if (!iep) {
        depthToNID.clear();
        // a bag count one level, other levels are traversed
        for (int i = pos; i < numAttributes - 1; ++i) {
            VertexID u = globalOrder[i];
            VertexID nID = v2n[u][0];
            depthToNID.push_back(nID);
            ++currentLevels[nID];
        }
        VertexID u = globalOrder.back();
        VertexID nID = v2n[u][0];
        levels = currentLevels;
        groups.clear();
        groups.emplace_back();
        groups[0].push_back(nID);
        adaptiveDepthToNID[extendLevel] = depthToNID;
        adaptiveGroups[extendLevel] = groups;
        adaptiveLevels[extendLevel] = levels;
        for (VertexID u2 = 0; u2 < numAttributes; ++u2) {
            if (u2 == u || query.getEdgeID(u2, u) != -1) continue;
            bool uLastHasRule = false;
            for (auto rule: symmetryRules) {
                if ((rule[0] == u2 && std::find(rule.begin(), rule.end(), u) != rule.end()) || (rule[0] == u && std::find(rule.begin(), rule.end(), u2) != rule.end()))
                    uLastHasRule = true;
            }
            if (uLastHasRule) continue;
            attributesToCheck.push_back(u2);
        }
        return;
    }
    VertexID last = globalOrder.back();
    for (VertexID u = 0; u < numAttributes; ++u) {
        if (u == last || query.getEdgeID(u, last) != -1) continue;
        bool uLastHasRule = false;
        for (auto rule: symmetryRules) {
            if ((rule[0] == u && std::find(rule.begin(), rule.end(), last) != rule.end()) || (rule[0] == last && std::find(rule.begin(), rule.end(), u) != rule.end()))
                uLastHasRule = true;
        }
        if (uLastHasRule) continue;
        attributesToCheck.push_back(u);
    }
    if (!symmLastLevel.empty()) {
        for (int i = extendLevel; i < numAttributes - symmLastLevel.size(); ++i) {
            VertexID u = globalOrder[i];
            VertexID nID = v2n[u][0];
            ++currentLevels[nID];
        }
    }
    else if (!subsetLastLevel.empty()) {
        for (int i = extendLevel; i < numAttributes - subsetLastLevel.size(); ++i) {
            VertexID u = globalOrder[i];
            VertexID nID = v2n[u][0];
            ++currentLevels[nID];
        }
        std::vector<VertexID> lastAttrs(subsetLastLevel.size());
        for (int i = 0; i < subsetLastLevel.size(); ++i) {
            VertexID nID = subsetLastLevel[i];
            VertexID u = nodes[nID].attributes[nodes[nID].numAttributes - 1];
            lastAttrs[i] = u;
        }
        subsetToCheck = std::vector<std::vector<VertexID>>(subsetLastLevel.size());
        for (int i = 0; i < subsetLastLevel.size(); ++i) {
            VertexID u = lastAttrs[i];
            for (VertexID u2 = 0; u2 < numAttributes; ++u2) {
                if (std::find(lastAttrs.begin(), lastAttrs.end(), u2) != lastAttrs.end()) continue;
                if (query.getEdgeID(u2, u) != -1) continue;
                bool uU2 = false;
                for (auto rule: symmetryRules) {
                    if ((rule[0] == u2 && std::find(rule.begin(), rule.end(), u) != rule.end()) || (rule[0] == u && std::find(rule.begin(), rule.end(), u2) != rule.end()))
                        uU2 = true;
                }
                if (uU2) continue;
                subsetToCheck[i].push_back(u2);
            }
        }
    }
    else {
//        int pos2 = 0;
//        globalSymmetry.clear();
//        for (const auto &rule: symmetryRules) {
//            bool local = false;
//            for (VertexID nID = 0; nID < numNodes; ++nID) {
//                bool include = true;
//                for (VertexID u: rule) {
//                    if (std::find(nodes[nID].attributes, nodes[nID].attributes + nodes[nID].numAttributes, u) == nodes[nID].attributes + nodes[nID].numAttributes) {
//                        include = false;
//                        break;
//                    }
//                }
//                if (include) {
//                    local = true;
//                    break;
//                }
//            }
//            if (!local) globalSymmetry.push_back(rule);
//        }
//        while (pos2 < notTraversed.size()) {
//            // for each rule, check whether 'notTraversed[pos2:]' is disjoint with it or contains it
//            // if so, break
//            bool flag = true;
//            divideFactor = 1;
//            for (auto rule: globalSymmetry) {
//                int numContained = 0;
//                for (int i = pos2; i < notTraversed.size(); ++i) {
//                    VertexID u = notTraversed[i];
//                    if (std::find(rule.begin(), rule.end(), u) != rule.end()) {
//                        ++numContained;
//                    }
//                }
//                if (!(numContained == 0 || numContained == rule.size() || numContained == rule.size() - 1)) flag = false;
//                if (numContained == rule.size() || numContained == rule.size() - 1) divideFactor *= rule.size();
//            }
//            if (flag) break;
//            ++pos2;
//        }
//        if (pos2 == notTraversed.size()) --pos2;
//        currentGroups.emplace_back();
//        for (int i = pos2; i < notTraversed.size(); ++i)
//            currentGroups[0].push_back(v2n[notTraversed[i]][0]);
//        for (int i = extendLevel; i < numAttributes - currentGroups[0].size(); ++i) {
//            VertexID u = globalOrder[i];
//            VertexID nID = v2n[u][0];
//            ++currentLevels[nID];
//        }
        currentGroups.push_back({v2n[globalOrder.back()][0]});
        for (int i = extendLevel; i < numAttributes - 1; ++i) {
            VertexID u = globalOrder[i];
            VertexID nID = v2n[u][0];
            ++currentLevels[nID];
        }
    }
    for (int i = extendLevel; i < numAttributes; ++i) {
        VertexID u = globalOrder[i];
        VertexID nID = v2n[u][0];
        currentDepthToNID.push_back(nID);
    }
    adaptiveLevels[extendLevel] = levels = currentLevels;
    adaptiveGroups[extendLevel] = groups = currentGroups;
    adaptiveDepthToNID[extendLevel] = depthToNID = currentDepthToNID;
    if (!newGlobalNode && extendLevel > nodes[numNodes - 1].prefixSize) {
        trieOrder[numNodes - 1].clear();
        for (int i = 0; i < globalOrder.size(); ++i) {
            VertexID u = globalOrder[i];
            ui k = extendLevel;
            for (; k < nodes[numNodes - 1].numAttributes; ++k) {
                if (nodes[numNodes - 1].attributes[k] == u)
                    break;
            }
            if (k != nodes[numNodes - 1].numAttributes) trieOrder[numNodes - 1].push_back(k - extendLevel);
        }
    }
}

bool HyperTree::hasSubNodeOf(const HyperNode &tau) const {
    bool flag = false;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        const HyperNode &bag = nodes[nID];
        ui numInTau = 0;
        for (int i = 0; i < bag.numAttributes; ++i) {
            VertexID u = bag.attributes[i];
            for (int j = 0; j < tau.numAttributes; ++j) {
                if (tau.attributes[j] == u) {
                    ++numInTau;
                    break;
                }
            }
        }
        if (numInTau == bag.numAttributes) {
            flag = true;
            break;
        }
    }

    return flag;
}

bool HyperTree::hasSupNodeOf(const HyperNode &tau) const {
    bool flag = false;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        const HyperNode &bag = nodes[nID];
        ui numInBag = 0;
        for (int i = 0; i < tau.numAttributes; ++i) {
            VertexID u = tau.attributes[i];
            for (int j = 0; j < bag.numAttributes; ++j) {
                if (bag.attributes[j] == u) {
                    ++numInBag;
                    break;
                }
            }
        }
        if (numInBag == tau.numAttributes) {
            flag = true;
            break;
        }
    }

    return flag;
}

HyperTree::HyperTree(ui numVertices) {
    numAttributes = numVertices;
    numNodes = 0;
    v2n = new std::vector<VertexID>[numVertices];
    nodes = new HyperNode[numVertices];
    edges = std::vector<std::vector<VertexID>>(numVertices);
    extendLevel = 0;
    divideFactor = 1;
    newGlobalNode = false;
}

void HyperTree::addPeripheral(const PatternGraph &p) {
    ui coreSize;
    VertexID *coreV = p.getCoreV(coreSize);
    ui peripheralSize = p.getNumVertices() - coreSize;
    VertexID *peripheral = p.getPeripheralV();
    // build the forests by bfs. each forest has exactly one edge connects to the coreV.
    std::vector<bool> visited(p.getNumVertices(), false);
    std::queue<Edge> q;
    std::map<VertexID, VertexID> nodeID;
    for (int i = 0; i < peripheralSize; ++i) {
        VertexID u = peripheral[i];
        if (visited[u]) continue;
        for (ui j = 0; j < coreSize; ++j)
            visited[coreV[j]] = false;
        q.emplace(99, u);   // use 99 as a virtual vertex
        // each time we pop an edge, we build a node for this edge and add the link(edge) to the node in the tree
        while (!q.empty()) {
            Edge e = q.front();
            q.pop();
            VertexID u1 = e.first, u2 = e.second;
            visited[u2] = true;
            bool flag = false;
            if (u1 != 99) {
                VertexID *vertices = new VertexID[2];
                vertices[0] = u1 < u2 ? u1 : u2;
                vertices[1] = u1 < u2 ? u2 : u1;
                int id = (1 << u1) + (1 << u2);
                HyperNode *old = nodes;
                nodes = new HyperNode[numNodes + 1];
                for (int j = 0; j < numNodes; ++j) {
                    old[j].copyTo(nodes[j]);
                }
                delete[] old;
                nodes[numNodes] = HyperNode(id, vertices, 2);
                // update _v2n
//                v2n[u1].push_back(numNodes);
//                v2n[u2].push_back(numNodes);
                nodeID[u2] = numNodes;
//                if (nodeID.find(u1) != nodeID.end()) {
//                    edges[numNodes].push_back(nodeID[u1]);
//                    edges[nodeID[u1]].push_back(numNodes);
//                }
                // check whether u2 is in coreV. If so, build an edge between node containing u2 and the last node
                for (int j = 0; j < coreSize; ++j) {
                    if (coreV[j] == u2) {
                        ui k;
                        for (k = 0; k < numNodes; ++k) {
                            if (nodes[k].hasVertex(u2))
                                break;
                        }
//                        edges[numNodes].push_back(k);
//                        edges[k].push_back(numNodes);
                        flag = true;
                        break;
                    }
                }
                ++numNodes;
            }
            // add neighbors only when u2 is not coreV
            if (!flag) {
                ui nbrCnt;
                VertexID *neighbors = p.getNeighbors(u2, nbrCnt);
                for (int j = 0; j < nbrCnt; ++j) {
                    VertexID next = neighbors[j];
                    if (!visited[next])
                        q.emplace(u2, next);
                }
            }
        }
    }
}

PrefixNode *buildPrefixTree(std::vector<std::vector<VertexID>> &orders, const Graph &query,
                            const std::vector<std::vector<size_t>> &dist,
                            std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> &visitedPT,
                            bool &exist) {
    PrefixNode *root = new PrefixNode(99);
    for (VertexID nID = 0; nID < orders.size(); ++nID) {
        if (orders[nID].empty()) continue;
        const std::vector<VertexID> &order = orders[nID];
        PrefixNode *currentNode = root;
        if (order[0] == 99) {
            root -> nIDsToCall.push_back(nID);
            if (nID == orders.size() - 1) currentNode -> pathToGlobal = true;
            continue;
        }
        for (int i = 0; i < order.size(); ++i) {
            VertexID u = order[i];
            int childIndex = currentNode->findChildIndex(u);
            currentNode->nIDsToCall.push_back(nID);
            if (nID != orders.size() - 1) currentNode -> pathToGlobal = false;
            else currentNode -> pathToGlobal = true;
            if (childIndex == -1) {
                // If the child does not exist, create a new child node
                PrefixNode *newNode = new PrefixNode(u);
                currentNode->children.push_back(newNode);
                currentNode = currentNode->children.back();
            } else {
                // If the child exists, move to the child node
                currentNode = currentNode->children[childIndex];
            }
            if (i == order.size() - 1) {
                currentNode->nIDsToCall.push_back(nID);
                if (nID != orders.size() - 1) currentNode -> pathToGlobal = false;
                else currentNode -> pathToGlobal = true;
            }
        }
    }
    std::queue<PrefixNode *> Q;
    Q.push(root);
    while (!Q.empty()) {
        PrefixNode *pn = Q.front();
        Q.pop();
        for (int i = 0; i < pn -> children.size(); ++i) {
            if (i != pn -> children.size() - 1 && pn -> children[i]->pathToGlobal) {
                std::swap(pn -> children[i], pn -> children.back());
                break;
            }
        }
        for (const auto &c : pn->children)
            Q.push(c);
    }
    // refine the prefix tree
    root -> refine(orders);
    if (visitedPT.find(root) == visitedPT.end()) exist = false;
    else {
        exist = true;
        delete root;
        return nullptr;
    }
    root -> initPoses(orders, query, dist);
    return root;
}

PrefixNode *fullAttributeTree(const HyperTree &t, std::vector<PrefixNode *> &attributeOrder,
                              std::map<PrefixNode *, std::vector<VertexID>> &bagsBelow, bool share) {
    PrefixNode *root = new PrefixNode(99);
    root->pathToGlobal = true;
    for (VertexID nID = 0; nID < t.numNodes - 1; ++nID) bagsBelow[root].push_back(nID);
    attributeOrder.push_back(root);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        bool global = nID == t.numNodes - 1;
        const HyperNode &bag = t.nodes[nID];
        PrefixNode *currentNode = root;
        for (int i = 0; i < bag.numAttributes; ++i) {
            VertexID u = bag.attributes[i];
            int childIndex = -1;
            if (share) childIndex = currentNode->findChildIndex(u);
            if (childIndex == -1) {
                if (global) break;
                // If the child does not exist, create a new child node
                PrefixNode *newNode = new PrefixNode(u);
                newNode->pathToGlobal = false;
                currentNode->children.push_back(newNode);
                currentNode = currentNode->children.back();
                attributeOrder.push_back(newNode);
                bagsBelow[newNode].push_back(nID);
            } else {
                // If the child exists, move to the child node
                currentNode = currentNode->children[childIndex];
                bagsBelow[currentNode].push_back(nID);
                if (global) currentNode->pathToGlobal = true;
            }
        }
    }

    return root;
}

PrefixNode *fullAttributeTree(PrefixNode *attrTree, std::vector<PrefixNode *> &attributeOrder,
                              std::map<PrefixNode *, std::vector<VertexID>> &bagsBelow, const HyperTree &t) {
    PrefixNode *root = attrTree->clone();
    PrefixNode *pn1 = attrTree;
    PrefixNode *pn2 = root;
    ui height = pn1 -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes1, nodes2;
    int depth = 0;
    attributeOrder.push_back(pn2);
    for (VertexID nID: pn1->nIDsToCall) {
        if (depth + 1 < t.nodes[nID].numAttributes && nID != t.numNodes - 1) {
            PrefixNode *attr = pn2;
            for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
                PrefixNode *newAttr = new PrefixNode(t.nodes[nID].attributes[i]);
                newAttr->pathToGlobal = false;
                attributeOrder.push_back(newAttr);
                attr->children.push_back(newAttr);
                attr = newAttr;
            }
            attr->nIDsToCall.push_back(nID);
        }
    }
    while (depth >= 0) {
        while (childPoses[depth] < pn1 -> children.size()) {
            PrefixNode *current1 = pn1 -> children[childPoses[depth]];
            PrefixNode *current2 = pn2 -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes1.push_back(current1);
            nodes2.push_back(current2);
            bool exists = false;
            for (PrefixNode *attr : attributeOrder) {
                if (attr == current2) {
                    exists = true;
                    break;
                }
            }
            if (!exists) attributeOrder.push_back(current2);
            for (VertexID nID: current1->nIDsToCall) {
                if (depth + 1 < t.nodes[nID].numAttributes && nID != t.numNodes - 1) {
                    PrefixNode *attr = current2;
                    for (int i = depth + 1; i < t.nodes[nID].numAttributes; ++i) {
                        PrefixNode *newAttr = new PrefixNode(t.nodes[nID].attributes[i]);
                        newAttr->pathToGlobal = false;
                        attributeOrder.push_back(newAttr);
                        attr->children.push_back(newAttr);
                        attr = newAttr;
                    }
                    attr->nIDsToCall.push_back(nID);
                }
            }
            if (!current1 -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else {
                nodes1.pop_back();
                nodes2.pop_back();
            }
            if (depth > 0) {
                pn1 = nodes1[depth - 1];
                pn2 = nodes2[depth - 1];
            }
            else {
                pn1 = attrTree;
                pn2 = root;
            }
        }
        --depth;
        if (depth >= 0) {
            nodes1.pop_back();
            nodes2.pop_back();
            if (depth == 0) {
                pn1 = attrTree;
                pn2 = root;
            }
            else {
                pn1 = nodes1[depth - 1];
                pn2 = nodes2[depth - 1];
            }
        }
    }
    std::vector<PrefixNode *> path = root->locate(t.numNodes - 1);
    PrefixNode *attr = root;
    if (!path.empty()) attr = path.back();
    for (int i = path.size(); i < t.nodes[t.numNodes - 1].numAttributes; ++i) {
        PrefixNode *newAttr = new PrefixNode(t.nodes[t.numNodes - 1].attributes[i]);
        newAttr->pathToGlobal = true;
        attributeOrder.push_back(newAttr);
        attr->children.push_back(newAttr);
        attr = newAttr;
    }
    attr->nIDsToCall.push_back(t.numNodes - 1);
    attrTree->addBagsBelow(bagsBelow);
    root->addBagsBelow(bagsBelow);

    return root;
}

std::vector<std::vector<int>>
buildAttrIDMap(const std::vector<PrefixNode *> &attributes, const std::map<PrefixNode *, std::vector<VertexID>> &bagsBelow,
               const HyperTree &t) {
    std::vector<std::vector<int>> attrIDMap(t.numNodes);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) attrIDMap[nID].resize(t.nodes[nID].numAttributes);
    int depth = 0;
    const PrefixNode *pn = attributes[0];
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes.push_back(current);
            int id = 0;
            for (int i = 0; i < attributes.size(); ++i) {
                if (attributes[i] == current) {
                    id = i;
                    break;
                }
            }
            for (VertexID nID : bagsBelow.at(current))
                attrIDMap[nID][depth] = id;
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = attributes[0];
        }
        --depth;
        if (depth >= 0) {
            nodes.pop_back();
            if (depth == 0) pn = attributes[0];
            else pn = nodes[depth - 1];
        }
    }

    return attrIDMap;
}

std::vector<std::pair<int, VertexID>> buildDependent(const std::vector<PrefixNode *> &attributes) {
    std::vector<std::pair<int, VertexID>> depend(attributes.size());
    depend[0].first = -1;
    depend[0].second = -1;
    int depth = 0;
    const PrefixNode *pn = attributes[0];
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes.push_back(current);
            int id = 0, pid = 0;
            for (int i = 0; i < attributes.size(); ++i) {
                if (attributes[i] == current) {
                    id = i;
                }
                if (attributes[i] == pn) {
                    pid = i;
                }
            }
            depend[id].first = pid;
            depend[id].second = -1;
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = attributes[0];
        }
        --depth;
        if (depth >= 0) {
            nodes.pop_back();
            if (depth == 0) pn = attributes[0];
            else pn = nodes[depth - 1];
        }
    }

    return depend;
}

// Function to clone a PrefixNode tree
PrefixNode* cloneTree(const PrefixNode* root) {
    if (!root) return nullptr;
    PrefixNode* newRoot = new PrefixNode(root->u);
    for (const auto& child : root->children) {
        newRoot->children.push_back(cloneTree(child));
    }
    return newRoot;
}

int PrefixNode::findChildIndex(VertexID id) {
    for (int i = 0; i < children.size(); ++i) {
        if (children[i]->u == id) {
            return i;
        }
    }
    return -1;
}

ui PrefixNode::getHeight() const {
    std::queue<const PrefixNode *> Q;
    Q.push(this);
    ui height = 0;
    while (!Q.empty()) {
        ++height;
        std::queue<const PrefixNode *> newQ;
        while (!Q.empty()) {
            const PrefixNode *pn = Q.front();
            Q.pop();
            for (const auto &c : pn->children)
                newQ.push(c);
        }
        Q = newQ;
    }

    return height;
}

std::vector<PrefixNode *> PrefixNode::locate(VertexID nID) const {
    int depth = 0;
    const PrefixNode *pn = this;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    for (VertexID nID2: pn -> nIDsToCall) {
        if (nID2 == nID) return nodes;
    }
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes.push_back(current);
            for (VertexID nID2: current -> nIDsToCall) {
                if (nID2 == nID) return nodes;
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = this;
        }
        --depth;
        nodes.pop_back();
        if (depth >= 0) {
            if (depth == 0) pn = this;
            else pn = nodes[depth - 1];
        }
    }

    return nodes;
}

void PrefixNode::locate(PrefixNode *attribute, std::vector<PrefixNode *> &path, std::vector<ui> &childPoses) {
    int depth = 0;
    PrefixNode *pn = this;
    ui height = pn -> getHeight();
    childPoses = std::vector<ui>(height, 0);
    if (pn == attribute) return;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            path.push_back(current);
            if (current == attribute) {
                return;
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else {
                ++childPoses[depth];
                path.pop_back();
            }
            if (depth > 0) pn = path[depth - 1];
            else pn = this;
        }
        --depth;
        if (depth >= 0) {
            ++childPoses[depth];
            path.pop_back();
            if (depth == 0) pn = this;
            else pn = path[depth - 1];
        }
    }
}

int PrefixNode::locate(PrefixNode *attribute, VertexID firstID){
    int depth = 0;
    PrefixNode *pn = this;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses = std::vector<ui>(height, 0);
    std::vector<PrefixNode *> path;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            path.push_back(current);
            if (current->u == attribute->u) {
                std::vector<VertexID> below = current->getBagsBelow();
                if (std::find(below.begin(), below.end(), firstID) != below.end())
                    return depth;
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else {
                ++childPoses[depth];
                path.pop_back();
            }
            if (depth > 0) pn = path[depth - 1];
            else pn = this;
        }
        --depth;
        if (depth >= 0) {
            ++childPoses[depth];
            path.pop_back();
            if (depth == 0) pn = this;
            else pn = path[depth - 1];
        }
    }

    return -1;
}

ui PrefixNode::numAttr(const HyperTree &t) const {
    ui num = 0;
    std::queue<const PrefixNode *> Q;
    std::queue<int> heights;
    Q.push(this);
    heights.push(0);
    while (!Q.empty()) {
        const PrefixNode *pn = Q.front();
        int height = heights.front();
        Q.pop();
        heights.pop();
        ++num;
        for (VertexID nID : pn -> nIDsToCall)
            num += t.nodes[nID].numAttributes - height;
        for (const auto &c : pn->children) {
            Q.push(c);
            heights.push(height + 1);
        }
    }

    return num - 1;
}

void PrefixNode::refine(const std::vector<std::vector<VertexID>> &sharedAttrs) {
    std::queue<PrefixNode *> Q;
    Q.push(this);
    while (!Q.empty()) {
        PrefixNode *pn = Q.front();
        Q.pop();
        std::vector<VertexID> newNID;
        std::set<VertexID> childCallNID;
        std::set<VertexID> buildNID;
        std::vector<int> validPoses;
        std::vector<PrefixNode *> newChildren;
        if (pn -> pathToGlobal) {
            for (VertexID nID = 0; nID < sharedAttrs.size() - 1; ++nID) {
                if (std::find(pn->nIDsToCall.begin(), pn->nIDsToCall.end(), nID) == pn->nIDsToCall.end()
                    && std::find(sharedAttrs[nID].begin(), sharedAttrs[nID].end(), pn->u) != sharedAttrs[nID].end()) {
                    pn->nIDsToJoin.push_back(nID);
                }
            }
        }
        for (int i = 0 ; i < pn->children.size(); ++i) {
            PrefixNode *c = pn->children[i];
            if (c->nIDsToCall.size() > 1) {
                Q.push(c);
                childCallNID.insert(c->nIDsToCall.begin(), c->nIDsToCall.end());
                validPoses.push_back(i);
                if (pn -> pathToGlobal && !c -> pathToGlobal) buildNID.insert(c->nIDsToCall.begin(), c->nIDsToCall.end());
            }
            else {
                delete c;
            }
        }
        for (auto pos: validPoses) {
            newChildren.push_back(pn -> children[pos]);
        }
        for (VertexID nID: pn -> nIDsToCall) {
            if (childCallNID.find(nID) == childCallNID.end())
                newNID.push_back(nID);
        }
        pn -> nIDsToCall = newNID;
        pn -> children = newChildren;
        if (pn -> pathToGlobal) {
            for (VertexID nID: pn -> nIDsToCall) {
                if (nID != sharedAttrs.size() - 1)
                    pn -> nIDsToBuild.push_back(nID);
            }
            for (VertexID nID : buildNID)
                pn -> nIDsToBuild.push_back(nID);
        }
    }
}

void PrefixNode::initPoses(const std::vector<std::vector<VertexID>> &bagAttrs, const Graph &query,
                           const std::vector<std::vector<size_t>> &dist) {
    std::vector<VertexID> attributes(query.getNumVertices());
    std::vector<PrefixNode *> nodes(query.getNumVertices(), nullptr);
    int depth = 0;
    std::vector<ui> childPoses(query.getNumVertices(), 0);
    PrefixNode *pn = this;
    // rebuild nIDsToJoin
    std::vector<VertexID> materializedNodes = pn -> nIDsToBuild;
    while (!pn->children.empty()) {
        PrefixNode *current = pn -> children.back();
        if (!(current -> pathToGlobal)) break;
        current -> nIDsToJoin.clear();
        for (VertexID nID : materializedNodes) {
            if (std::find(bagAttrs[nID].begin(), bagAttrs[nID].end(), current->u) != bagAttrs[nID].end()) {
                current -> nIDsToJoin.push_back(nID);
            }
        }
        for (VertexID nID : current -> nIDsToBuild) {
            materializedNodes.push_back(nID);
        }
        pn = current;
    }
    depth = 0;
    pn = this;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u1 = current -> u;
            attributes[depth] = u1;
            current -> attributesBefore.clear();
            for (int i = 0; i < depth; ++i) {
                VertexID u2 = attributes[i];
                if (query.getEdgeID(u1, u2) != -1)
                    current -> attributesBefore.push_back(u2);
            }
            if (current -> pathToGlobal) {
                std::set<VertexID> coveredAttrs;
                for (VertexID nID : current -> nIDsToJoin) {
                    for (VertexID u2 : bagAttrs[nID])
                        coveredAttrs.insert(u2);
                }
                std::vector<ui> attrsToJoin;
                for (VertexID u2: current -> attributesBefore) {
                    if (coveredAttrs.find(u2) == coveredAttrs.end())
                        attrsToJoin.push_back(u2);
                }
                current -> attributesBefore = attrsToJoin;
            }
            bool cartesian = (depth > 0 && current -> attributesBefore.empty() && current -> nIDsToJoin.empty());
            if (cartesian) {
                size_t minDist = std::numeric_limits<size_t>::max();
                for (VertexID i = 0; i < depth; ++i) {
                    VertexID u2 = attributes[i];
                    if (dist[u1][u2] < minDist) {
                        minDist = dist[u1][u2];
                        current -> cartesianParent = u2;
                    }
                }
            }
            nodes[depth] = current;
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
                pn = nodes[depth - 1];
            }
        }
        --depth;
        if (depth == 0) pn = this;
        else if (depth > 0) pn = nodes[depth - 1];
    }
}

void PrefixNode::print() const {
    std::queue<const PrefixNode *> Q;
    Q.push(this);
    while (!Q.empty()) {
        const PrefixNode *pn = Q.front();
        Q.pop();
        if (!pn -> children.empty()) std::cout << pn -> u << " -- ";
        for (int i = 0; i < pn -> children.size(); ++i) {
            std::cout << " " << pn -> children[i] -> u << " ";
        }
        if (!pn -> children.empty()) std::cout << std::endl;
        for (const auto &c : pn->children)
            Q.push(c);
    }
}

bool PrefixNode::operator==(const PrefixNode &rhs) const {
    if (children.size() != rhs.children.size()) return false;
    for (int i = 0; i < children.size(); ++i) {
        if (*children[i] != *rhs.children[i]) return false;
    }
    return u == rhs.u && nIDsToCall == rhs.nIDsToCall;
}

bool PrefixNode::operator!=(const PrefixNode &rhs) const {
    return !(rhs == *this);
}

void PrefixNode::checkCallAll(ui numNodes) const {
    std::set<VertexID> nIDs;
    std::queue<const PrefixNode *> Q;
    Q.push(this);
    while (!Q.empty()) {
        const PrefixNode *pn = Q.front();
        Q.pop();
        nIDs.insert(pn -> nIDsToCall.begin(), pn->nIDsToCall.end());
        for (const auto &c : pn->children)
            Q.push(c);
    }
//    assert(nIDs.size() == numNodes);
}

std::vector<VertexID> PrefixNode::getBagsBelow() const {
    std::vector<VertexID> nIDs;
    std::stack<const PrefixNode *> S;
    S.push(this);
    while (!S.empty()) {
        const PrefixNode *pn = S.top();
        S.pop();
        for (VertexID nID : pn -> nIDsToCall)
            nIDs.push_back(nID);
        for (auto it = pn->children.rbegin(); it != pn->children.rend(); ++it)
            S.push(*it);
    }

    return nIDs;
}

PrefixNode *PrefixNode::clone() const {
    PrefixNode* newNode = new PrefixNode(u);
    newNode->attributesBefore = attributesBefore;
    newNode->cartesianParent = cartesianParent;
    newNode->nIDsToCall = nIDsToCall;
    newNode->nIDsToJoin = nIDsToJoin;
    newNode->nIDsToBuild = nIDsToBuild;
    newNode->pathToGlobal = pathToGlobal;

    for (const auto& child : children) {
        newNode->children.push_back(child->clone());
    }

    return newNode;
}

void PrefixNode::initNIDsToBuild(ui numNodes) {
    PrefixNode *pn = this;
    // find the path to global
    std::queue<PrefixNode *> q1;
    q1.push(pn);
    ui height = getHeight();
    std::vector<PrefixNode *> path, nodes(height, nullptr);
    int depth = 0;
    std::vector<ui> childPoses(height, 0);
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes[depth] = current;
            current -> pathToGlobal = false;
            current -> nIDsToBuild.clear();
            current -> nIDsToJoin.clear();
            for (VertexID nID: current -> nIDsToCall) {
                if (nID == numNodes - 1) {
                    path.assign(nodes.begin(), nodes.begin() + depth + 1);
                    break;
                }
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = this;
            else pn = nodes[depth - 1];
        }
    }
    this -> pathToGlobal = true;
    for (PrefixNode *node : path) node -> pathToGlobal = true;
    pn = this;
    std::queue<PrefixNode *> Q;
    Q.push(this);
    while (!Q.empty()) {
        pn = Q.front();
        Q.pop();
        for (int i = 0; i < pn -> children.size(); ++i) {
            if (i != pn -> children.size() - 1 && pn -> children[i]->pathToGlobal) {
                std::swap(pn -> children[i], pn -> children.back());
                break;
            }
        }
        for (const auto &c : pn->children)
            Q.push(c);
    }
    pn = this;
    while (true) {
        pn -> nIDsToBuild.clear();
        for (VertexID nID : pn -> nIDsToCall)
            if (nID != numNodes - 1)
                pn -> nIDsToBuild.push_back(nID);
        std::queue<PrefixNode *>q;
        for (PrefixNode *c : pn -> children) {
            if (!c->pathToGlobal)
                q.push(c);
        }
        while (!q.empty()) {
            PrefixNode *current = q.front();
            q.pop();
            for (VertexID nID : current -> nIDsToCall) pn -> nIDsToBuild.push_back(nID);
            for (PrefixNode *c : current -> children) q.push(c);
        }
        std::sort(pn -> nIDsToBuild.begin(), pn->nIDsToBuild.end());
        if (pn -> children.empty() || !pn->children.back()->pathToGlobal) break;
        else pn = pn -> children.back();
    }
}

void PrefixNode::refineNIDsToBuild(const HyperTree &t) {
    std::vector<VertexID> notBuild;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (t.trieOrder[nID].size() < 2)
            notBuild.push_back(nID);
    }
    PrefixNode *pn = this;
    // find the path to global
    std::queue<PrefixNode *> q1;
    q1.push(pn);
    ui height = getHeight();
    std::vector<PrefixNode *> path, nodes(height, nullptr);
    int depth = 0;
    std::vector<ui> childPoses(height, 0);
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes[depth] = current;
            std::vector<VertexID> newNIDsToBuild;
            for (VertexID nID : current -> nIDsToBuild) {
                if (std::find(notBuild.begin(), notBuild.end(), nID) == notBuild.end())
                    newNIDsToBuild.push_back(nID);
            }
            current->nIDsToBuild = newNIDsToBuild;
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = this;
            else pn = nodes[depth - 1];
        }
    }
}

// check whether the left subtree can be merged to the right path
void PrefixNode::mergeToRight(std::vector<std::vector<VertexID>> &localOrders, std::vector<VertexID> &rightCall) {
    if (this -> children.empty()) return;
    PrefixNode *pn = this;
    PrefixNode *right = this -> children.back();
    // merge bags to the right
    std::vector<VertexID> nIDsToRemove;
    std::sort(nIDsToCall.begin(), nIDsToCall.end());
    for (int i = 0; i < nIDsToCall.size(); ++i) {
        VertexID nID = nIDsToCall[i];
        if (localOrders[nID][0] == right->u) {
            nIDsToRemove.push_back(nID);
            rightCall.push_back(nID);
            std::sort(rightCall.begin(), rightCall.end());
        }
    }
    std::vector<VertexID> old = nIDsToCall;
    nIDsToCall.clear();
    std::set_difference(old.begin(), old.end(), nIDsToRemove.begin(), nIDsToRemove.end(), std::back_inserter(nIDsToCall));
    // merge subtree to the right
    int pos = -1;
    for (int i = 0; i < children.size() - 1; ++i) {
        if (children[i] -> u == right->u) {
            pos = i;
            break;
        }
    }
    if (pos == -1) return;
    PrefixNode *left = this -> children[pos];
    std::queue<PrefixNode *> q;
    q.push(left);
    while (!q.empty()) {
        PrefixNode *subtree = q.front();
        q.pop();
        for (VertexID nID: subtree->nIDsToCall) rightCall.push_back(nID);
        for (PrefixNode *c: subtree->children) q.push(c);
    }
    std::sort(rightCall.begin(), rightCall.end());
    std::vector<PrefixNode *> newChild;
    for (int i = 0; i < children.size(); ++i) {
        if (i != pos) newChild.push_back(children[i]);
    }
    children = newChild;
}

void PrefixNode::getTraverseOrder(std::vector<PrefixNode *> &attributes, std::vector<VertexID> &nIDs, const HyperTree &t) {
    int depth = 0;
    const PrefixNode *pn = this;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    for (VertexID nID2: pn -> nIDsToCall) {
        if (nID2 !=t.numNodes - 1)
            nIDs.push_back(nID2);
    }
    attributes.push_back(this);
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            attributes.push_back(current);
            ++childPoses[depth];
            nodes.push_back(current);
            for (VertexID nID2: current -> nIDsToCall) {
                if (nID2 !=t.numNodes - 1) nIDs.push_back(nID2);
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = this;
        }
        --depth;
        if (depth >= 0) {
            nodes.pop_back();
            if (depth == 0) pn = this;
            else pn = nodes[depth - 1];
        }
    }
    nIDs.push_back(t.numNodes - 1);
}

void PrefixNode::addBagsBelow(std::map<PrefixNode *, std::vector<VertexID>> &bagsBelow) {
    int depth = 0;
    PrefixNode *pn = this;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    bagsBelow[pn] = pn->getBagsBelow();
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            bagsBelow[current] = current->getBagsBelow();
            ++childPoses[depth];
            nodes.push_back(current);
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = this;
        }
        --depth;
        if (depth >= 0) {
            nodes.pop_back();
            if (depth == 0) pn = this;
            else pn = nodes[depth - 1];
        }
    }
}

void PrefixNode::addBagsBelow(std::map<const PrefixNode *, std::vector<VertexID>> &bagsBelow, VertexID lastID) const {
    int depth = 0;
    const PrefixNode *pn = this;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    bagsBelow[pn] = pn->getBagsBelow();
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            bagsBelow[current] = current->getBagsBelow();
            for (int i = 0; i < bagsBelow[current].size(); ++i) {
                if (bagsBelow[current][i] == lastID) {
                    std::swap(bagsBelow[current][i], bagsBelow[current][bagsBelow[current].size() - 1]);
                }
            }
            ++childPoses[depth];
            nodes.push_back(current);
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = this;
        }
        --depth;
        if (depth >= 0) {
            nodes.pop_back();
            if (depth == 0) pn = this;
            else pn = nodes[depth - 1];
        }
    }
}

void buildFromPrefixTree(PrefixNode *prefixTree, const std::vector<std::vector<VertexID>> &nodeOrder, HyperTree &t,
                         const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query, CandidateSpace &cs) {
    PrefixNode *pn = prefixTree;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    int depth = 0;
    std::vector<PrefixNode *> nodes(query.getNumVertices(), nullptr);
    t.globalOrder.clear();
    for (VertexID nID = 0; nID < t.numNodes; ++nID) t.nodes[nID].prefixSize = 0;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            nodes[depth] = current;
            if (current -> pathToGlobal && u != 99) t.globalOrder.push_back(u);
            for (VertexID nID: current -> nIDsToCall) {
                t.nodes[nID].numAttributes = nodeOrder[nID].size();
                t.nodes[nID].attributes = new VertexID [nodeOrder[nID].size()];
                for (int j = 0; j <= depth; ++j) {
                    if (nodes[j]->pathToGlobal) {
                        VertexID u2 = nodes[j] -> u;
                        if (std::find(t.v2n[u2].begin(), t.v2n[u2].end(), nID) != t.v2n[u2].end())
                            ++t.nodes[nID].prefixSize;
                    }
                    else break;
                }
                t.nodes[nID].prefix = new VertexID [t.nodes[nID].prefixSize];
                for (int j = 0; j < t.nodes[nID].prefixSize; ++j) {
                    t.nodes[nID].prefix[j] = nodeOrder[nID][j];
                }
                for (int j = 0; j < t.nodes[nID].numAttributes; ++j) {
                    t.nodes[nID].attributes[j] = nodeOrder[nID][j];
                }
                t.nodes[nID].initPoses(sharedAttrs, query, cs.dist, nID == t.numNodes - 1);
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = prefixTree;
            else pn = nodes[depth - 1];
        }
    }
    for (VertexID nID: prefixTree -> nIDsToCall) {
        t.nodes[nID].numAttributes = nodeOrder[nID].size();
        t.nodes[nID].attributes = new VertexID [nodeOrder[nID].size()];
        t.nodes[nID].prefixSize = 0;
        for (int j = 0; j < t.nodes[nID].numAttributes; ++j) {
            t.nodes[nID].attributes[j] = nodeOrder[nID][j];
        }
        t.nodes[nID].initPoses(sharedAttrs, query, cs.dist, nID == t.numNodes - 1);
    }
    for (int i = 0; i < nodeOrder.back().size(); ++i) {
        VertexID u = nodeOrder.back()[i];
        if (std::find(t.globalOrder.begin(), t.globalOrder.end(), u) == t.globalOrder.end())
            t.globalOrder.push_back(u);
    }
    for (VertexID nID = 0; nID < t.numNodes - 1; ++nID) {
        for (int j = 0; j < nodeOrder[nID].size(); ++j) {
            VertexID u = nodeOrder[nID][j];
            if (std::find(t.globalOrder.begin(), t.globalOrder.end(), u) == t.globalOrder.end())
                t.globalOrder.push_back(u);
        }
    }
    t.extendLevel = t.nodes[t.numNodes - 1].prefixSize;
    if (!cs.labeled && t.extendLevel < t.nodes[t.numNodes - 1].numAttributes - 1)
        t.extendLevel = t.nodes[t.numNodes - 1].numAttributes - 1;
    t.nIDs = std::vector<std::vector<VertexID>>(t.globalOrder.size());
    for (int j = 0; j < t.globalOrder.size(); ++j) {
        VertexID u = t.globalOrder[j];
        for (VertexID nID = 0; nID < nodeOrder.size(); ++nID) {
            if (std::find(nodeOrder[nID].begin(), nodeOrder[nID].end(), u) != nodeOrder[nID].end())
                t.nIDs[j].push_back(nID);
        }
    }
    t.buildTrieOrder();
}

void buildFromPrefixTree(const std::vector<PrefixNode *> &prefixTrees,
                         const std::vector<std::vector<std::vector<VertexID>>> &bestOrders,
                         std::vector<HyperTree> &trees, const HyperTree &reference,
                         const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query,
                         CandidateSpace &cs) {
    trees.resize(prefixTrees.size());
    for (int i = 0; i < prefixTrees.size(); ++i) {
        HyperTree &t = trees[i];
        t.numAttributes = reference.numAttributes;
        t.numNodes = reference.numNodes;
        t.nodes = new HyperNode[t.numNodes];
        t.v2n = new std::vector<VertexID> [t.numAttributes];
        for (VertexID u = 0; u < t.numAttributes; ++u)
            t.v2n[u] = reference.v2n[u];
        buildFromPrefixTree(prefixTrees[i], bestOrders[i], t, sharedAttrs, query, cs);
    }
}

void removeNID(VertexID nID, PrefixNode *root) {
    for (auto it = root -> nIDsToBuild.begin(); it != root -> nIDsToBuild.end();) {
        if (*it == nID) {
            root -> nIDsToBuild.erase(it);
            break;
        }
        else ++it;
    }
    for (auto it = root -> nIDsToCall.begin(); it != root -> nIDsToCall.end();) {
        if (*it == nID) {
            root -> nIDsToCall.erase(it);
            break;
        }
        else ++it;
    }
}

// remove nID in the branch: the nIDsToBuild for the global node, and the nIDsToCall for the current node
// for nodes in the path, if there is no bags to call, it is also removed
void removeNID(VertexID nID, std::vector<PrefixNode *> &nodes, int depth, PrefixNode *&root) {
    int pos = 0;
    for (; pos < depth + 1; ++pos) {
        if (!nodes[pos]->pathToGlobal)
            break;
    }
    PrefixNode *pn = root;
    if (pos != 0 ) pn = nodes[pos - 1];
    PrefixNode *current = nodes[depth];
    for (auto it = pn -> nIDsToBuild.begin(); it != pn -> nIDsToBuild.end();) {
        if (*it == nID) {
            pn -> nIDsToBuild.erase(it);
            break;
        }
        else ++it;
    }
    for (auto it = current -> nIDsToCall.begin(); it != current -> nIDsToCall.end();) {
        if (*it == nID) {
            current -> nIDsToCall.erase(it);
            break;
        }
        else ++it;
    }
    if (current -> children.empty() && current -> nIDsToCall.empty()) {
        pos = depth - 1;
        PrefixNode *child = current;
        current = root;
        if (pos != -1) current = nodes[pos];
        for (auto it = current -> children.begin(); it != current -> children.end(); ++it) {
            if (*it == child) {
                current -> children.erase(it);
                break;
            }
        }
        while (current -> children.empty() && current -> nIDsToCall.empty()) {
            child = current;
            --pos;
            current = root;
            if (pos != -1) current = nodes[pos];
            for (auto it = current -> children.begin(); it != current -> children.end(); ++it) {
                if (*it == child) {
                    current -> children.erase(it);
                    break;
                }
            }
        }
    }
}

std::vector<int> getMappingSizes(const HyperTree &t, const PrefixNode *pt) {
    std::vector<int> mappingSizes(t.numNodes);
    std::vector<VertexID> attrsInPath;
    int depth = 0;
    const PrefixNode *pn = pt;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes(height, nullptr);
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            attrsInPath.push_back(u);
            nodes[depth] = current;
            for (VertexID nID: current -> nIDsToCall) {
                for (VertexID u2: attrsInPath) {
                    if (std::find(t.v2n[u2].begin(), t.v2n[u2].end(), nID) != t.v2n[u2].end())
                        ++mappingSizes[nID];
                }
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

    return mappingSizes;
}

// Function to convert a subset to an ID
uint64_t getSubsetID(const std::vector<VertexID>& subset) {
    uint64_t id = 0;
    for (VertexID vertex : subset) {
        id |= (1ULL << vertex);
    }
    return id;
}

bool subsetConnectivity(const Graph &query, const std::vector<uint64_t> &cc, const std::vector<VertexID> &subset) {
    if (subset.empty()) return true;
    std::vector<std::vector<VertexID>> subsetCC;
    query.computeConnectedComponents(subset, subsetCC);
    if (subsetCC.size() == 1) return true;
    else return false;
}

bool orderConnectivity(const Graph &query, const std::vector<VertexID> &order, const std::vector<uint64_t> &cc) {
    if (order.empty()) return false;
    uint64_t ccID = 1 << order[0];
    for (int i = 1; i < order.size(); ++i) {
        bool connected = false;
        VertexID u1 = order[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = order[j];
            if (query.getEdgeID(u1, u2) != -1) {
                connected = true;
                ccID += 1 << u1;
                break;
            }
        }
        if (!connected) {
            bool valid = false;
            for (auto queryCC: cc) {
                if (queryCC == ccID) {
                    valid = true;
                    break;
                }
            }
            if (!valid) return false;
            ccID = 1 << u1;
        }
    }
    return true;
}

bool orderConnectivity(const Graph &query, const std::vector<VertexID> &order) {
    for (int i = 1; i < order.size(); ++i) {
        VertexID u1 = order[i];
        bool connected = false;
        for (int j = 0; j < i; ++j) {
            VertexID u2 = order[j];
            if (query.getEdgeID(u1, u2) != -1) {
                connected = true;
                break;
            }
        }
        if (!connected) return false;
    }

    return true;
}

void
genAllPrefixTree(const HyperTree &t, const Graph &query, CandidateSpace &cs, std::vector<PrefixNode *> &prefixTrees) {
    std::vector<ui> poses(query.getNumVertices(), 0);
    std::vector<std::vector<VertexID>> v2n(query.getNumVertices());
    ui numNodes = t.numNodes;
    std::vector<std::vector<VertexID>> noShareOrders(numNodes);
    std::vector<std::vector<ui>> numBackWards(numNodes);
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            v2n[u].push_back(nID);
            noShareOrders[nID].push_back(u);
        }
    }
    std::vector<std::vector<VertexID>> sharedAttrs(t.numNodes);
    std::vector<VertexID> internalAttrs;
    std::vector<std::vector<VertexID>> attrIntersections(t.numNodes * t.numNodes);
    std::vector<ui> repetitions(query.getNumVertices(), 0) ;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            ++repetitions[u];
            if (v2n[u].size() > 1) {
                sharedAttrs[nID].push_back(u);
            }
        }
        std::sort(sharedAttrs[nID].begin(), sharedAttrs[nID].end());
    }
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        if (v2n[u].size() > 1)
            internalAttrs.push_back(u);
    }
    std::vector<std::vector<VertexID>> components;
    std::vector<uint64_t> cc;
    query.computeConnectedComponents(internalAttrs, components);
    for (auto &c : components) cc.push_back(getSubsetID(c));
    for (VertexID nID1 = 0; nID1 < numNodes; ++nID1) {
        for (VertexID nID2 = nID1 + 1; nID2 < numNodes; ++nID2) {
            std::set_intersection(sharedAttrs[nID1].begin(), sharedAttrs[nID1].end(), sharedAttrs[nID2].begin(),
                                  sharedAttrs[nID2].end(), std::back_inserter(attrIntersections[nID1 * t.numNodes + nID2]));
        }
    }
    std::unordered_set<PrefixNode*, PrefixNodePtrHash, PrefixNodePtrEqual> visitedPT;
    PrefixNode *empty = new PrefixNode(99);
    for (VertexID nID = 0; nID < numNodes; ++nID) empty->nIDsToCall.push_back(nID);
    std::stack<PrefixNode *> plans;
    plans.push(empty);
    while (!plans.empty()) {
        PrefixNode *pt = plans.top();
        plans.pop();
        bool valid = true;
        for (VertexID nID: pt->nIDsToCall) {
            if (t.nodes[nID].numAttributes >= 3 && !t.nodes[nID].isTriangle(query)) valid = false;
        }
        if (!pt->children.empty()) {
            for (int i = 0; i < pt->children.size() - 1; ++i) {
                for (VertexID nID: pt->children[i]->getBagsBelow()) {
                    if (t.nodes[nID].numAttributes >= 3 && !t.nodes[nID].isTriangle(query)) valid = false;
                }
            }
        }
        if (valid) prefixTrees.push_back(pt);
        if (prefixTrees.size() >= 1000) break;
        PrefixNode *pn = pt;
        ui height = pn -> getHeight();
        std::vector<ui> childPoses(height, 0);
        std::vector<PrefixNode *> nodes(height, nullptr);
        int depth = 0;
        std::vector<VertexID> attrsInPath;
        std::vector<VertexID> siblingAttr;
        extendPrefixTree(query, t, plans, pt, pn, attrsInPath, siblingAttr, -1, childPoses,
                         visitedPT, attrIntersections, repetitions, sharedAttrs, cc);
        while (depth >= 0) {
            while (childPoses[depth] < pn -> children.size()) {
                PrefixNode *current = pn -> children[childPoses[depth]];
                VertexID u = current -> u;
                attrsInPath.push_back(u);
                nodes[depth] = current;
                siblingAttr.clear();
                if (depth == 0) {
                    for (PrefixNode *sibling: pt -> children)
                        if (sibling != current) siblingAttr.push_back(sibling -> u);
                }
                else {
                    for (PrefixNode *sibling: nodes[depth - 1] -> children)
                        if (sibling != current) siblingAttr.push_back(sibling -> u);
                }
                ++childPoses[depth];
                extendPrefixTree(query, t, plans, pt, current, attrsInPath, siblingAttr, depth, childPoses,
                                 visitedPT, attrIntersections, repetitions, sharedAttrs, cc);
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
    }
    for (PrefixNode *pt : prefixTrees) {
        pt->initNIDsToBuild(t.numNodes);
        pt ->initPoses(sharedAttrs, query, cs.dist);
    }
}

void extendPrefixTree(const Graph &query, const HyperTree &t, std::stack<PrefixNode *> &plans, const PrefixNode *pt,
                      const PrefixNode *current, const std::vector<VertexID> &attrsInPath,
                      const std::vector<VertexID> &siblingAttr, int depth, std::vector<ui> &childPoses,
                      std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> &visitedPT,
                      const std::vector<std::vector<VertexID>> &attrIntersections, const std::vector<ui> &repetitions,
                      const std::vector<std::vector<VertexID>> &sharedAttrs, const std::vector<uint64_t> &cc) {
    VertexID u = current -> u;
    if (!current -> children.empty()) {
        for (int i = 0; i < current -> children.size(); ++i) {
            PrefixNode *c = current -> children[i];
            VertexID u2 = c -> u;
            for (VertexID nID : current -> nIDsToCall) {
                if (std::find(sharedAttrs[nID].begin(), sharedAttrs[nID].end(), u2) == sharedAttrs[nID].end()) continue;
                PrefixNode *newPT = pt -> clone();
                PrefixNode *newPN = newPT;
                for (int j = 0; j <= depth; ++j) {
                    newPN = newPN -> children[childPoses[j] - 1];
                }
                newPN -> nIDsToCall.clear();
                for (VertexID nID2 : current -> nIDsToCall) {
                    if (nID != nID2) newPN -> nIDsToCall.push_back(nID2);
                }
                newPN -> children[i] -> nIDsToCall.push_back(nID);
                std::sort(newPN -> children[i] -> nIDsToCall.begin(), newPN->children[i]->nIDsToCall.end());
                if (visitedPT.find(newPT) != visitedPT.end()) {
                    delete newPT;
                    continue;
                }
                else {
                    visitedPT.insert(newPT);
                    plans.push(newPT);
                }
            }
        }
    }
    if (current -> nIDsToCall.size() >= 2) {
        std::vector<std::vector<bool>> choices = chooseK(current->nIDsToCall.size(), 2);
        for (auto choice: choices) {
            std::vector<VertexID> selected;
            for (int j = 0; j < current->nIDsToCall.size(); ++j) {
                VertexID nID = current->nIDsToCall[j];
                if (choice[j]) selected.push_back(nID);
            }
            ui pos = selected[0] * t.numNodes + selected[1];
            if (selected[0] > selected[1]) pos = selected[1] * t.numNodes + selected[0];
            std::vector<VertexID> attrIntersect = attrIntersections[pos];
            std::sort(attrIntersect.begin(), attrIntersect.end(), [&repetitions](const VertexID a, const VertexID b) {
                if (repetitions[a] < repetitions[b]) return true;
                else if (repetitions[a] == repetitions[b] && a > b) return true;
                return false;
            });
            for (VertexID u2: attrIntersect) {
                if (std::find(siblingAttr.begin(), siblingAttr.end(), u2) != siblingAttr.end()) continue;
                if (std::find(attrsInPath.begin(), attrsInPath.end(), u2) != attrsInPath.end()) continue;
                std::vector<VertexID> prefix = attrsInPath;
                bool connected = prefix.empty();
                for (VertexID up: prefix) {
                    if (query.getEdgeID(up, u2) != -1) {
                        connected = true;
                        break;
                    }
                }
                if (!connected) continue;
                prefix.push_back(u2);
                PrefixNode *newPT = pt -> clone();
                PrefixNode *newPN = newPT;
                for (int i = 0; i <= depth; ++i) {
                    newPN = newPN -> children[childPoses[i] - 1];
                }
                PrefixNode *newChild = new PrefixNode(u2);
                newChild -> nIDsToCall = selected;
                if (selected[0] == t.numNodes - 1 || selected[1] == t.numNodes - 1) newPN->children.push_back(newChild);
                else {
                    std::vector<PrefixNode *> old = newPN -> children;
                    newPN -> children = {newChild};
                    for (PrefixNode * c: old) newPN -> children.push_back(c);
                }
                newPN -> nIDsToCall.clear();
                for (int i = 0; i < current->nIDsToCall.size(); ++i) {
                    if (!choice[i]) newPN -> nIDsToCall.push_back(current -> nIDsToCall[i]);
                }
                if (visitedPT.find(newPT) != visitedPT.end()) {
                    delete newPT;
                    continue;
                }
                else {
                    visitedPT.insert(newPT);
                    plans.push(newPT);
                }
            }
        }
    }
}

std::vector<VertexID> globalOrder(const Graph &query, const HyperTree &t, const CandidateSpace &cs, const std::vector<VertexID> &prefix) {
    std::vector<VertexID> localOrder;
    std::vector<ui> repetitions(query.getNumVertices(), 0);
    std::vector<ui> sizes(query.getNumVertices(), 0);
    const HyperNode &global = t.nodes[t.numNodes - 1];
    for (int i = 0; i < global.numAttributes; ++i) {
        VertexID u = global.attributes[i];
        if (std::find(prefix.begin(), prefix.end(), u) == prefix.end())
            localOrder.push_back(u);
    }
    for (VertexID u : localOrder) {
        repetitions[u] = t.v2n[u].size();
        sizes[u] = cs.candidateSet[u].size();
    }

    std::sort(localOrder.begin(), localOrder.end(), [&repetitions, &sizes](const VertexID &a, const VertexID &b) {
        if (repetitions[a] > repetitions[b]) return true;
        if (repetitions[a] == repetitions[b] && sizes[a] < sizes[b]) return true;
        return false;
    });
    std::vector<VertexID> totalOrder = prefix;
    for (VertexID u : localOrder) totalOrder.push_back(u);
    return totalOrder;
}

void swapBags(HyperTree &t, VertexID nID1, VertexID nID2) {
    if (nID1 == nID2) return;
    HyperNode tmp;
    tmp.numAttributes = t.nodes[nID1].numAttributes;
    tmp.attributes = new VertexID[tmp.numAttributes];
    memcpy(tmp.attributes, t.nodes[nID1].attributes, sizeof(VertexID) * tmp.numAttributes);
    t.nodes[nID1].numAttributes = t.nodes[nID2].numAttributes;
    delete[] t.nodes[nID1].attributes;
    t.nodes[nID1].attributes = new VertexID [t.nodes[nID1].numAttributes];
    memcpy(t.nodes[nID1].attributes, t.nodes[nID2].attributes, sizeof(VertexID) * t.nodes[nID2].numAttributes);
    t.nodes[nID2].numAttributes = tmp.numAttributes;
    delete[] t.nodes[nID2].attributes;
    t.nodes[nID2].attributes = new VertexID [t.nodes[nID2].numAttributes];
    memcpy(t.nodes[nID2].attributes, tmp.attributes, sizeof(VertexID) * tmp.numAttributes);
    std::swap(t.nodes[nID1].id, t.nodes[nID2].id);
    std::swap(t.nodes[nID1].fw, t.nodes[nID2].fw);
    std::swap(t.nodes[nID1].canonValue, t.nodes[nID2].canonValue);
    if (!t.symmetryRules.empty()) {
        std::swap(t.nodes[nID1].smallerAttrs, t.nodes[nID2].smallerAttrs);
        std::swap(t.nodes[nID1].largerAttrs, t.nodes[nID2].largerAttrs);
    }
    for (int i = 0; i < t.nodes[nID2].numAttributes; ++i) {
        VertexID u = t.nodes[nID2].attributes[i];
        for (int j = 0; j < t.v2n[u].size(); ++j) {
            if (t.v2n[u][j] == nID1)
                t.v2n[u][j] = nID2;
        }
    }
    for (int i = 0; i < t.nodes[nID1].numAttributes; ++i) {
        VertexID u = t.nodes[nID1].attributes[i];
        for (int j = 0; j < t.v2n[u].size(); ++j) {
            if (t.v2n[u][j] == nID2)
                t.v2n[u][j] = nID1;
        }
    }
}

void addGlobalBag(const Graph &query, HyperTree &t) {
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
    std::vector<VertexID> globalCand;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (std::includes(sharedAttrs[nID].begin(), sharedAttrs[nID].end(), globalAttr.begin(), globalAttr.end())) {
            globalBag = false;
            globalCand.push_back(nID);
        }
    }
    if (!query.isConnected(globalAttr)) globalBag = true;
    if (!globalBag) {
        VertexID nID = globalCand[0];
        ui maxNum = t.nodes[nID].numAttributes;
        for (int i = 1; i < globalCand.size(); ++i) {
            ui num = t.nodes[globalCand[i]].numAttributes;
            if (num > maxNum) {
                maxNum = num;
                nID = globalCand[i];
            }
        }
        if (nID != t.numNodes - 1) {
            swapBags(t, nID, t.numNodes - 1);
            std::swap(sharedAttrs[nID], sharedAttrs[t.numNodes - 1]);
        }
    }
    else {
        sharedAttrs.push_back(globalAttr);
        t.addGlobalNode(globalAttr);
    }
    delete[] t.v2n;
    t.v2n = new std::vector<VertexID> [query.getNumVertices()];
    for (ui i = 0; i < t.numNodes; ++i) {
        for (ui j = 0; j < t.nodes[i].numAttributes; ++j) {
            t.v2n[t.nodes[i].attributes[j]].push_back(i);
        }
    }
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
}

void changeNIDsToCall(const Graph &query, HyperTree &t, PrefixNode *pt) {
    int depth = 0;
    PrefixNode *pn = pt;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    std::map<PrefixNode *, std::vector<PrefixNode *>> newChild;
    // root
    for (VertexID nID: pn ->nIDsToCall) {
        if (!t.newGlobalNode && nID == t.numNodes - 1) continue;
        VertexID u = t.nodes[nID].attributes[0];
        PrefixNode *child = new PrefixNode(u);
        child->pathToGlobal = false;
        child->attributesBefore = t.nodes[nID].attributesBefore[0];
        child->nIDsToCall = {nID};
        newChild[pn].push_back(child);
    }
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes.push_back(current);
            for (VertexID nID: current -> nIDsToCall) {
                if (!t.newGlobalNode && nID == t.numNodes - 1) continue;
                VertexID u = t.nodes[nID].attributes[depth + 1];
                PrefixNode *child = new PrefixNode(u);
                child->pathToGlobal = false;
                child->attributesBefore = t.nodes[nID].attributesBefore[depth + 1];
                child->nIDsToCall = {nID};
                newChild[current].push_back(child);
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = pt;
        }
        --depth;
        if (depth >= 0) {
            nodes.pop_back();
            if (depth == 0) pn = pt;
            else pn = nodes[depth - 1];
        }
    }
    for (auto &item : newChild) {
        PrefixNode *attr = item.first;
        std::vector<PrefixNode *> bags = item.second;
        std::vector<PrefixNode *> newChild = bags;
        for (PrefixNode *c: attr->children) newChild.push_back(c);
        attr->children = newChild;
        if (!t.newGlobalNode && std::find(attr->nIDsToCall.begin(), attr->nIDsToCall.end(),
                                          t.numNodes - 1) != attr->nIDsToCall.end()) {
            attr->nIDsToCall = {t.numNodes - 1};
        }
        else attr->nIDsToCall.clear();
    }
}

void addGlobal(const Graph &query, HyperTree &t, PrefixNode *pt, VertexID nID, VertexID maxCostID, CandidateSpace &cs,
               std::vector<std::vector<VertexID>> &nodeOrders, std::vector<ui> &prefixSizes, bool globalOrderShare) {
    t.largerAttrs.emplace_back(t.nodes[t.numNodes - 1].numAttributes);
    t.smallerAttrs.emplace_back(t.nodes[t.numNodes - 1].numAttributes);
    std::vector<PrefixNode *> attributes2 = pt->locate(maxCostID);
    std::vector<PrefixNode *> attributes;
    for (PrefixNode *attribute: attributes2) {
        if (t.v2n[attribute->u].size() >= 2) attributes.push_back(attribute);
        else break;
    }
    std::vector<VertexID> prefix;
    for (int i = 0; i < attributes.size(); ++i) prefix.push_back(attributes[i]->u);
    std::vector<VertexID> newShareAttr;
    int pos = attributes.size();
    for (; pos < nodeOrders[maxCostID].size(); ++pos) {
        VertexID u = nodeOrders[maxCostID][pos];
        if (t.v2n[u].size() > 1) {
            newShareAttr.push_back(u);
            prefix.push_back(u);
        }
        else break;
    }
    nodeOrders[nID] = globalOrder(query, t, cs, prefix);
    PrefixNode *attr = pt;
    if (!attributes.empty()) attr = attributes.back();
    for (VertexID u : newShareAttr) {
        PrefixNode *newAttr = new PrefixNode(u);
        attr -> children.push_back(newAttr);
        attr = newAttr;
    }
    if (!prefixSizes.empty()) {
        for (int i = 0; i < attributes.size(); ++i) {
            for (VertexID nID2 : attributes[i]->getBagsBelow())
                prefixSizes[nID2] = i + 1;
        }
        prefixSizes[nID] = prefixSizes[maxCostID] = prefix.size();
    }
    if (!globalOrderShare) {
        newShareAttr.clear();
        prefix.clear();
        nodeOrders[nID] = globalOrder(query, t, cs, prefix);
    }
    if (newShareAttr.empty()) attr->nIDsToCall.push_back(nID);
    else {
        PrefixNode *pn = pt;
        if (!attributes.empty()) pn = attributes.back();
        for (auto it = pn->nIDsToCall.begin(); it != pn->nIDsToCall.end(); ++it) {
            if (*it == maxCostID) {
                pn->nIDsToCall.erase(it);
                break;
            }
        }
        attr->nIDsToCall = {maxCostID, nID};
    }
}

void refineIntersectionInfo(const Graph &query, HyperTree &t, PrefixNode *pt, bool &type) {
    std::vector<ui> prefixSum(t.numNodes + 1, 0);
    for (VertexID nID = 1; nID <= t.numNodes; ++nID) {
        prefixSum[nID] = prefixSum[nID - 1] + t.nodes[nID - 1].numAttributes;
    }
    std::vector<PrefixNode *> attributeOrder;
    std::vector<VertexID> bagOrder;
    pt->getTraverseOrder(attributeOrder, bagOrder, t);
    std::vector<std::vector<VertexID>> oldBefore(prefixSum.back());
    std::vector<int> candPos(prefixSum.back());
    std::vector<bool> nIDsEmpty(prefixSum.back(), true);
    for (int i = 0; i < prefixSum.back(); ++i) {
        candPos[i] = i;
    }
    PrefixNode *pn = pt;
    int depth = -1;
    while (!pn->children.empty()) {
        pn = pn->children.back();
        ++depth;
        if (!t.nodes[t.numNodes - 1].nIDs[depth].empty()) {
            std::vector<VertexID> below = pn->getBagsBelow();
            for (VertexID nID: below) {
                nIDsEmpty[prefixSum[nID] + depth] = false;
            }
        }
    }
    for (int i = 0; i < bagOrder.size(); ++i) {
        VertexID nID = bagOrder[i];
        HyperNode &bag = t.nodes[nID];
        bag.candidatesBefore = std::vector<int>(bag.numAttributes, 0);
        for (int j = 0; j < bag.numAttributes; ++j) {
            oldBefore[prefixSum[nID] + j] = bag.attributesBefore[j];
        }
        // find the largest cover
        for (int j = 2; j < bag.numAttributes; ++j) {
            if (bag.cartesianParent[j] != 99 || (!bag.nIDs.empty() && !bag.nIDs[j].empty())) continue;
            std::vector<VertexID> oldAttrBefore = bag.attributesBefore[j];
            std::vector<VertexID> oldSmallerAttrs = bag.smallerAttrs[j];
            std::vector<VertexID> oldLargerAttrs = bag.largerAttrs[j];
            int maxCoverSize = 1, maxPos, maxPrevious = -1;
            bool equal = false;
            // previous bags
            if (std::includes(bag.attributes, bag.attributes + bag.prefixSize, oldAttrBefore.begin(),
                              oldAttrBefore.end())) {
                for (int l = 0; l < i; ++l) {
                    HyperNode &previous = t.nodes[bagOrder[l]];
                    for (int k = 2; k < previous.numAttributes; ++k) {
                        if (!nIDsEmpty[prefixSum[bagOrder[l]] + k]) continue;
                        if (orderedSubset(oldBefore[prefixSum[bagOrder[l]] + k], oldAttrBefore) &&
                            unorderedSubset(previous.largerAttrs[k], oldLargerAttrs) &&
                            unorderedSubset(previous.smallerAttrs[k], oldSmallerAttrs)) {
                            if (oldBefore[prefixSum[bagOrder[l]] + k].size() > maxCoverSize) {
                                if (oldBefore[prefixSum[bagOrder[l]] + k].size() == oldAttrBefore.size()) {
                                    equal = true;
                                }
                                else equal = false;
                                maxCoverSize = oldBefore[prefixSum[bagOrder[l]] + k].size();
                                maxPos = k;
                                maxPrevious = bagOrder[l];
                            }
                        }
                    }
                }
            }
            // current bag
            for (int k = 2; k < j; ++k) {
                if (!nIDsEmpty[prefixSum[nID] + k]) continue;
                if (orderedSubset(oldBefore[prefixSum[nID] + k], oldAttrBefore) &&
                    unorderedSubset(bag.largerAttrs[k], oldLargerAttrs) &&
                    unorderedSubset(bag.smallerAttrs[k], oldSmallerAttrs)) {
                    if (oldBefore[prefixSum[nID] + k].size() > maxCoverSize) {
                        if (oldBefore[prefixSum[nID] + k].size() == oldAttrBefore.size()) equal = true;
                        else equal = false;
                        maxCoverSize = oldBefore[prefixSum[nID] + k].size();
                        maxPos = k;
                        maxPrevious = nID;
                    }
                }
            }
            if (maxCoverSize == 1) continue;
            if (maxPrevious == -1) maxPrevious = nID;
            std::unordered_set<VertexID> exists(oldBefore[prefixSum[maxPrevious] + maxPos].begin(),
                                                oldBefore[prefixSum[maxPrevious] + maxPos].end());
            bag.attributesBefore[j].clear();
            for (VertexID u: oldAttrBefore) {
                if (exists.find(u) == exists.end())
                    bag.attributesBefore[j].push_back(u);
            }
            bag.candidatesBefore[j] = candPos[prefixSum[maxPrevious] + maxPos];
            if (equal) candPos[prefixSum[nID] + j] = candPos[prefixSum[maxPrevious] + maxPos];
        }
    }
    if (!type) {
        std::vector<std::vector<std::vector<VertexID>>> oldLarger(t.numNodes), oldSmaller(t.numNodes);
        for (VertexID nID = 0; nID < t.numNodes; ++nID) {
            oldLarger[nID] = t.nodes[nID].largerAttrs;
            oldSmaller[nID] = t.nodes[nID].smallerAttrs;
        }
        depth = 0;
        pn = pt;
        ui height = pn -> getHeight();
        std::vector<ui> childPoses(height, 0);
        std::vector<PrefixNode *> nodes;
        std::map<const PrefixNode *, std::vector<VertexID>> bagsBelow;
        std::vector<std::vector<bool>> copied(t.numNodes, std::vector<bool>(t.numAttributes, false));
        std::vector<std::vector<int>> reverse(t.numNodes, std::vector<int>(t.numAttributes, -1));
        for (VertexID nID = 0; nID < t.numNodes; ++nID) {
            for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
                VertexID u = t.nodes[nID].attributes[i];
                reverse[nID][u] = i;
            }
        }
        pt->addBagsBelow(bagsBelow, t.numNodes - 1);
        for (VertexID nID = 0; nID < t.numNodes; ++nID) {
            t.nodes[nID].candidatesAfter.resize(t.nodes[nID].numAttributes);
            t.nodes[nID].attributesAfter.resize(t.nodes[nID].numAttributes);
            t.nodes[nID].copyAfter.resize(t.nodes[nID].numAttributes);
            t.nodes[nID].attrAfterLarger.resize(t.nodes[nID].numAttributes);
            t.nodes[nID].attrAfterSmaller.resize(t.nodes[nID].numAttributes);
            t.nodes[nID].copyAfterTypes.resize(t.nodes[nID].numAttributes);
        }
        t.candIndex.resize(t.numNodes);
        int startIndex = 0;
        for (VertexID nID = 0; nID < t.numNodes; ++nID) {
            t.candIndex[nID].resize(t.nodes[nID].numAttributes);
            for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
                ui segmentSize = (t.nodes[nID].candidatesBefore[i] != 0) + t.nodes[nID].attributesBefore[i].size();
                if (segmentSize == 0) {
                    t.candIndex[nID][i] = -1;
                    continue;
                }
                startIndex += segmentSize;
                t.candIndex[nID][i] = startIndex - 1;
            }
        }
        std::vector<int> startIndexes(t.numNodes, 0);
        for (VertexID nID: bagOrder) {
            for (int i = 0; i < t.nodes[nID].numAttributes; ++i)
                if (!setAfter(t, nID, reverse, prefixSum, copied, i, startIndexes[nID], oldLarger, oldSmaller)) {
                    type = true;
                    return;
                }
            std::vector<PrefixNode *> path = pn -> locate(nID);
            for (int i = 0; i < path.size(); ++i) {
                for (VertexID nID2: bagsBelow[path[i]])
                    startIndexes[nID2] = i + 1;
            }
        }
    }
    t.refineSymmetryAttrs();
}

bool intersectType(const Graph &query, const HyperTree &t) {
    bool type = true;
    ui maxNumAttributes = 0;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        const HyperNode &bag = t.nodes[nID];
        if (bag.numAttributes > maxNumAttributes)
            maxNumAttributes = bag.numAttributes;
    }
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        const HyperNode &bag = t.nodes[nID];
        if (bag.numAttributes < 4 || bag.numAttributes < maxNumAttributes) continue;
        ui maxDegree = 0;
        bool clique = true;
        for (int i = 0; i < bag.numAttributes; ++i) {
            ui degree = 0;
            for (int j = 0; j < bag.numAttributes; ++j) {
                if (i == j) continue;
                VertexID u1 = bag.attributes[i], u2 = bag.attributes[j];
                if (query.getEdgeID(u1, u2) != -1)
                    ++degree;
            }
            if (degree > maxDegree) maxDegree = degree;
            if (degree != bag.numAttributes - 1) clique = false;
        }
        if (clique) continue;
        VertexID last = bag.attributes[bag.numAttributes - 1];
        ui degree = 0;
        for (int i = 0; i < bag.numAttributes - 1; ++i) {
            VertexID u = bag.attributes[i];
            if (query.getEdgeID(u, last) != -1)
                ++degree;
        }
        if (degree == maxDegree) type = false;
    }

    return type;
}

// 计算(nID, i)在allCandidates中的起始位置
ui getAllCandidatesStartIndex(const HyperTree &t, VertexID targetBag, int targetIndex) {
    ui startIndex = 0;
    // 累加bags 0到targetBag-1的所有segments
    for (VertexID nID = 0; nID < targetBag; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            ui segmentSize = (t.nodes[nID].candidatesBefore[i] != 0) + t.nodes[nID].attributesBefore[i].size();
            startIndex += segmentSize;
        }
    }
    // 累加targetBag中index 0到targetIndex-1的segments
    for (int i = 0; i < targetIndex; ++i) {
        ui segmentSize = (t.nodes[targetBag].candidatesBefore[i] != 0) + t.nodes[targetBag].attributesBefore[i].size();
        startIndex += segmentSize;
    }
    return startIndex;
}

bool setAfter(HyperTree &t, VertexID bag, const std::vector<std::vector<int>> &reverse,
              const std::vector<ui> &prefixSum, std::vector<std::vector<bool>> &copied, int depth, int start,
              const std::vector<std::vector<std::vector<VertexID>>> &oldLarger,
              const std::vector<std::vector<std::vector<VertexID>>> &oldSmaller) {
    VertexID u = t.nodes[bag].attributes[depth];
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        int sharedLength = 0;
        while (t.nodes[bag].attributes[sharedLength] == t.nodes[nID].attributes[sharedLength]) {
            ++sharedLength;
            if (sharedLength == t.nodes[bag].numAttributes || sharedLength == t.nodes[nID].numAttributes)
                break;
        }
        for (int i = 2; i < t.nodes[nID].numAttributes; ++i) {
            if (t.nodes[nID].candidatesBefore[i] == prefixSum[bag] + depth && prefixSum[bag] + depth != 0) {
                ui allCandStartIndex = getAllCandidatesStartIndex(t, nID, i);
                if (!t.nodes[nID].attributesBefore[i].empty() && depth >= sharedLength) {
                    return false;
                }
                else {
                    t.nodes[bag].candidatesAfter[depth].push_back(allCandStartIndex);
                    copied[nID][i] = true;
                }
                for (VertexID u3 : oldLarger[bag][depth]) {
                    auto it = std::find(t.nodes[nID].largerAttrs[i].begin(), t.nodes[nID].largerAttrs[i].end(), u3);
                    if (it != t.nodes[nID].largerAttrs[i].end()) {
                        t.nodes[nID].largerAttrs[i].erase(it);
                    }
                }
                for (VertexID u3 : oldSmaller[bag][depth]) {
                    auto it = std::find(t.nodes[nID].smallerAttrs[i].begin(), t.nodes[nID].smallerAttrs[i].end(), u3);
                    if (it != t.nodes[nID].smallerAttrs[i].end()) {
                        t.nodes[nID].smallerAttrs[i].erase(it);
                    }
                }
            }
        }
    }
    int pos1 = start > depth + 1 ? start : depth + 1;
    for (int i = pos1; i < t.nodes[bag].numAttributes; ++i) {
        if (bag == t.numNodes - 1 && !t.nodes[bag].nIDs[i].empty()) continue;
        if (std::find(t.nodes[bag].attributesBefore[i].begin(), t.nodes[bag].attributesBefore[i].end(), u) != t.nodes[bag].attributesBefore[i].end()) {
            VertexID u2 = t.nodes[bag].attributes[i];
            if (std::find(t.nodes[bag].attributesBefore[i].begin(), t.nodes[bag].attributesBefore[i].end(), u) == t.nodes[bag].attributesBefore[i].end())
                continue;
            bool copy = false;
            if (!copied[bag][i]) copy = true;
            if (copy) {
                int type = 0;
                ui allCandStartIndex = getAllCandidatesStartIndex(t, bag, i);
                auto it = std::find(t.nodes[bag].largerAttrs[i].begin(), t.nodes[bag].largerAttrs[i].end(), u);
                if (it != t.nodes[bag].largerAttrs[i].end()) {
                    t.nodes[bag].largerAttrs[i].erase(it);
                    type = 1;
                }
                it = std::find(t.nodes[bag].smallerAttrs[i].begin(), t.nodes[bag].smallerAttrs[i].end(), u);
                if (it != t.nodes[bag].smallerAttrs[i].end()) {
                    t.nodes[bag].smallerAttrs[i].erase(it);
                    type = 2;
                }
                t.nodes[bag].copyAfter[depth].push_back(allCandStartIndex);
                t.nodes[bag].copyAfterTypes[depth].push_back(type);
                copied[bag][i] = true;
            }
            else {
                ui allCandStartIndex = getAllCandidatesStartIndex(t, bag, i);
                // 找到u在u2的attributesBefore中的位置
                auto it = std::find(t.nodes[bag].attributesBefore[i].begin(), t.nodes[bag].attributesBefore[i].end(), u);
                int uPosInU2 = it - t.nodes[bag].attributesBefore[i].begin();
                ui nextPos;
                if (t.nodes[bag].candidatesBefore[i] == 0) {
                    nextPos = allCandStartIndex + uPosInU2;
                } else {
                    nextPos = allCandStartIndex + uPosInU2 + 1;
                }
                t.nodes[bag].attributesAfter[depth].push_back(nextPos);
                t.nodes[bag].attrAfterLarger[depth].emplace_back();
                t.nodes[bag].attrAfterSmaller[depth].emplace_back();
                for (VertexID u3: oldLarger[bag][i]) {
                    if (reverse[bag][u3] <= depth) {
                        t.nodes[bag].attrAfterLarger[depth].back().push_back(u3);
                        // Remove u3 from largerAttrs[depth] if it exists
                        auto it = std::find(t.nodes[bag].largerAttrs[i].begin(), t.nodes[bag].largerAttrs[i].end(), u3);
                        if (it != t.nodes[bag].largerAttrs[i].end()) {
                            t.nodes[bag].largerAttrs[i].erase(it);
                        }
                    }
                }
                for (VertexID u3: oldSmaller[bag][i]) {
                    if (reverse[bag][u3] <= depth) {
                        t.nodes[bag].attrAfterSmaller[depth].back().push_back(u3);
                        // Remove u3 from smallerAttrs[depth] if it exists
                        auto it = std::find(t.nodes[bag].smallerAttrs[i].begin(), t.nodes[bag].smallerAttrs[i].end(), u3);
                        if (it != t.nodes[bag].smallerAttrs[i].end()) {
                            t.nodes[bag].smallerAttrs[i].erase(it);
                        }
                    }
                }
            }
        }
    }

    return true;
}