//
// Created by anonymous authors on 2024/2/27.
//

#ifndef ASDMATCH_DECOMPOSITION_H
#define ASDMATCH_DECOMPOSITION_H

#include "config.h"
#include "graph.h"
#include "candidate_space.h"
#include <stack>
#include <fstream>
#include <sstream>

struct HyperNode {
    CanonType id;
    CanonType canonValue;
    VertexID *attributes;   // sorted by the matching order
    ui numAttributes;
    std::vector<std::vector<VertexID>> attributesBefore;
    std::vector<VertexID> cartesianParent;
    std::vector<std::vector<VertexID>> nIDs; // used for the global join
    VertexID *prefix;
    ui prefixSize;
    double fw;
    /******* used for symmetry-breaking in unlabeled graphs ********/
    std::vector<std::vector<VertexID>> largerAttrs;
    std::vector<std::vector<VertexID>> smallerAttrs;
    std::vector<bool> hasLargerAttrs;
    std::vector<bool> hasSmallerAttrs;
    std::vector<int> candidatesBefore;
    std::vector<std::vector<int>> attributesAfter;
    std::vector<std::vector<std::vector<VertexID>>> attrAfterLarger;
    std::vector<std::vector<std::vector<VertexID>>> attrAfterSmaller;
    std::vector<std::vector<int>> copyAfter;
    std::vector<std::vector<int>> copyAfterTypes;
    std::vector<std::vector<int>> candidatesAfter;

    HyperNode(): id(0), canonValue(0), attributes(nullptr), prefix(nullptr), numAttributes(0), prefixSize(0), fw(0.0) {}

    HyperNode(CanonType id, VertexID *attributes, ui numAttributes);

    HyperNode(const HyperNode &other);

    HyperNode& operator=(const HyperNode &other);

    ~HyperNode() {
        delete[] attributes;
        delete[] prefix;
    }

    void initPoses(const Graph &query, const std::vector<std::vector<size_t>> &dist);
    void initPoses(const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query, const std::vector<std::vector<size_t>> &dist, bool global);
    void copyTo(HyperNode &other) const;
    bool hasVertex(VertexID u) const;
    bool isTriangle(const Graph &query) const;
};

struct HyperTree {
    HyperNode *nodes;
    std::vector<VertexID>* v2n;           // for each query vertex, store the nodes that contains it.
    ui numAttributes;
    ui numNodes;
    std::vector<std::vector<VertexID>> edges;
    std::vector<ui> compressionSizes;
    std::vector<std::vector<VertexID>> nIDs;
    int extendLevel;
    std::vector<VertexID> globalOrder;
    std::vector<std::vector<VertexID>> trieOrder;
    /*******only used for scope-style join********/
    std::vector<std::vector<VertexID>> nodesAtStep;
    // dimension 1 is nID, dimension 2 is mapping size, dimension 3 are the children to call
    std::vector<VertexID> defaultPartition;
    std::vector<std::vector<VertexID>> attributesBefore;
    bool newGlobalNode;    // whether there is a virtual global node
    std::vector<VertexID> cartesianParent;
    /*******only used for global join and enumerate********/
    std::vector<VertexID> depthToNID;
    std::vector<ui> levels;
    std::vector<int> trieParents;
    std::vector<std::vector<VertexID>> groups;
    std::map<int, std::vector<VertexID>> adaptiveDepthToNID;
    std::map<int, std::vector<ui>> adaptiveLevels;
    std::map<int, std::vector<std::vector<VertexID>>> adaptiveGroups;
    /******* used for symmetry-breaking in unlabeled graphs ********/
    std::vector<std::vector<VertexID>> largerAttrs;
    std::vector<std::vector<VertexID>> smallerAttrs;
    std::vector<bool> hasLargerAttrs;
    std::vector<bool> hasSmallerAttrs;
    std::vector<std::vector<VertexID>> symmetryRules;
    int divideFactor;
    std::vector<VertexID> attributesToCheck;
    std::vector<VertexID> symmLastLevel;
    std::vector<VertexID> subsetLastLevel;
    std::vector<std::vector<VertexID>> subsetToCheck;
    std::vector<std::vector<int>> candIndex;
    /*********************************************/

    HyperTree(): nodes(nullptr), v2n(nullptr), numAttributes(0), numNodes(0), extendLevel(0), divideFactor(1), newGlobalNode(false) {};
    HyperTree(const HyperTree &other);
    HyperTree& operator=(const HyperTree &other);
    explicit HyperTree(ui numVertices);
    ~HyperTree() {
        delete[] nodes;
        delete[] v2n;
    }

    void initPoses(const Graph &query, const CandidateSpace &cs, bool handleNode);
    void buildTrieOrder();
    void buildTraverseStruct(const Graph &query);
    void buildTraverseAdaptive(const Graph &query, int currentExtendLevel);
    void buildTraverseUnlabeled(const Graph &query, CandidateSpace &cs, bool iep);
    void buildFromTD(FHD &fhd);
    void print(const Graph &query) const;
    void addGlobalNode(const std::vector<VertexID> &globalAttrs);
    void writeToStream(std::ostream &outStream);
    void selectSymmetry(const PatternGraph &p);
    void setSymmetry();
    void refineSymmetryAttrs();
    bool hasSubNodeOf(const HyperNode &tau) const;
    bool hasSupNodeOf(const HyperNode &tau) const;
    void addPeripheral(const PatternGraph &p);
};

// the attribute tree structure
struct PrefixNode {
    VertexID u;
    std::vector<VertexID> attributesBefore;
    VertexID cartesianParent;
    // hyper nodes to call in this node
    std::vector<VertexID> nIDsToCall;
    // hyper nodes to join, for the prefix nodes in the path to the global join
    std::vector<VertexID> nIDsToJoin;
    // hyper nodes to build tries
    std::vector<VertexID> nIDsToBuild;
    std::vector<PrefixNode *> children;
    bool pathToGlobal;
    PrefixNode(VertexID id) : u(id), cartesianParent(99), pathToGlobal(true) {}
    int findChildIndex(VertexID id);

    PrefixNode(const PrefixNode& other) = delete;
    PrefixNode& operator=(const PrefixNode &other) = delete;
    ~PrefixNode() {
        for (auto & c : children)
            delete c;
    }

    ui getHeight() const;
    std::vector<PrefixNode *> locate(VertexID nID) const;
    void locate(PrefixNode *attribute, std::vector<PrefixNode *> &path, std::vector<ui> &childPoses);
    int locate(PrefixNode *attribute, VertexID firstID);
    void refine(const std::vector<std::vector<VertexID>> &sharedAttrs);
    void initPoses(const std::vector<std::vector<VertexID>> &bagAttrs, const Graph &query, const std::vector<std::vector<size_t>> &dist);
    void initNIDsToBuild(ui numNodes);
    void refineNIDsToBuild(const HyperTree &t);
    void print() const;
    void checkCallAll(ui numNodes) const;
    std::vector<VertexID> getBagsBelow() const;
    ui numAttr(const HyperTree &t) const;

    bool operator==(const PrefixNode &rhs) const;

    bool operator!=(const PrefixNode &rhs) const;
    PrefixNode *clone() const;
    void mergeToRight(std::vector<std::vector<VertexID>> &localOrders, std::vector<VertexID> &rightCall);
    void getTraverseOrder(std::vector<PrefixNode *> &attributes, std::vector<VertexID> &nIDs, const HyperTree &t);
    void addBagsBelow(std::map<PrefixNode *, std::vector<VertexID>> &bagsBelow);
    void addBagsBelow(std::map<const PrefixNode *, std::vector<VertexID>> &bagsBelow, VertexID lastID) const;
};

struct PrefixNodePtrHash {
    std::size_t operator()(const PrefixNode* ptr) const {
        return std::hash<int>()(ptr->u);
    }
};

struct PrefixNodePtrEqual {
    bool operator()(const PrefixNode* lhs, const PrefixNode* rhs) const {
        return *lhs == *rhs;
    }
};

PrefixNode *buildPrefixTree(std::vector<std::vector<VertexID>> &orders, const Graph &query,
                            const std::vector<std::vector<size_t>> &dist,
                            std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> &visitedPT,
                            bool &exist);
PrefixNode *fullAttributeTree(const HyperTree &t, std::vector<PrefixNode *> &attributeOrder,
                              std::map<PrefixNode *, std::vector<VertexID>> &bagsBelow, bool share);
PrefixNode *fullAttributeTree(PrefixNode *attrTree, std::vector<PrefixNode *> &attributeOrder,
                              std::map<PrefixNode *, std::vector<VertexID>> &bagsBelow, const HyperTree &t);
std::vector<std::vector<int>>
buildAttrIDMap(const std::vector<PrefixNode *> &attributes, const std::map<PrefixNode *, std::vector<VertexID>> &bagsBelow,
               const HyperTree &t);
std::vector<std::pair<int, VertexID>> buildDependent(const std::vector<PrefixNode *> &attributes);
PrefixNode* cloneTree(const PrefixNode* root);

void buildFromPrefixTree(PrefixNode *prefixTree, const std::vector<std::vector<VertexID>> &nodeOrder, HyperTree &t,
                         const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query, CandidateSpace &cs);
void buildFromPrefixTree(const std::vector<PrefixNode *> &prefixTrees,
                         const std::vector<std::vector<std::vector<VertexID>>> &bestOrders,
                         std::vector<HyperTree> &trees, const HyperTree &reference,
                         const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query, CandidateSpace &cs);
void removeNID(VertexID nID, PrefixNode *root);
void removeNID(VertexID nID, std::vector<PrefixNode *> &nodes, int depth, PrefixNode *&root);
std::vector<int> getMappingSizes(const HyperTree &t, const PrefixNode *pt);
uint64_t getSubsetID(const std::vector<VertexID>& subset);
bool subsetConnectivity(const Graph &query, const std::vector<uint64_t> &cc, const std::vector<VertexID> &subset);
bool orderConnectivity(const Graph &query, const std::vector<VertexID> &order, const std::vector<uint64_t> &cc);
bool orderConnectivity(const Graph &query, const std::vector<VertexID> &order);
void
genAllPrefixTree(const HyperTree &t, const Graph &query, CandidateSpace &cs, std::vector<PrefixNode *> &prefixTrees);
void extendPrefixTree(const Graph &query, const HyperTree &t, std::stack<PrefixNode *> &plans, const PrefixNode *pt,
                      const PrefixNode *current, const std::vector<VertexID> &attrsInPath,
                      const std::vector<VertexID> &siblingAttr, int depth, std::vector<ui> &childPoses,
                      std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> &visitedPT,
                      const std::vector<std::vector<VertexID>> &attrIntersections, const std::vector<ui> &repetitions,
                      const std::vector<std::vector<VertexID>> &sharedAttrs, const std::vector<uint64_t> &cc);
std::vector<VertexID> globalOrder(const Graph &query, const HyperTree &t, const CandidateSpace &cs, const std::vector<VertexID> &prefix);
void swapBags(HyperTree &t, VertexID nID1, VertexID nID2);
void addGlobalBag(const Graph &query, HyperTree &t);
void changeNIDsToCall(const Graph &query, HyperTree &t, PrefixNode *pt);
void addGlobal(const Graph &query, HyperTree &t, PrefixNode *pt, VertexID nID, VertexID maxCostID, CandidateSpace &cs,
               std::vector<std::vector<VertexID>> &nodeOrders, std::vector<ui> &prefixSizes, bool globalOrderShare);
void refineIntersectionInfo(const Graph &query, HyperTree &t, PrefixNode *pt, bool &type);
std::vector<int>
getPrefixAttrID(PrefixNode *pt, const std::vector<VertexID> &prefix, const std::vector<std::vector<int>> &attrIDMap,
                const std::vector<PrefixNode *> &dynamicPartition);
bool intersectType(const Graph &query, const HyperTree &t);
bool setAfter(HyperTree &t, VertexID bag, const std::vector<std::vector<int>> &reverse,
              const std::vector<ui> &prefixSum, std::vector<std::vector<bool>> &copied, int depth, int start,
              const std::vector<std::vector<std::vector<VertexID>>> &oldLarger,
              const std::vector<std::vector<std::vector<VertexID>>> &oldSmaller);


#endif //ASDMATCH_DECOMPOSITION_H
