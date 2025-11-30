//
// Created by on 25-6-8.
//

#ifndef GRAPH_MINING_SYSTEM_ALL_TD_H
#define GRAPH_MINING_SYSTEM_ALL_TD_H

#include "decomposition.h"
#include "../utility/td/td.h"
#include "../utility/td/SCsolver.h"

std::vector<HyperNode> getAllNodes(const PatternGraph &p);
std::vector<HyperNode> getCandidateNodes(const HyperTree &t, const std::vector<HyperNode> &allNodes, const PatternGraph &p);
void findCliquesRecursive(const PatternGraph &graph,
                          std::vector<VertexID> &currentClique,
                          std::vector<VertexID> &potentialClique,
                          std::vector<VertexID> &processedVertices,
                          std::queue<HyperNode> &cliques);
std::queue<HyperNode> findMaximalCliques(const PatternGraph &graph);
std::vector<HyperTree> getAllTree(const PatternGraph &p);
std::vector<HyperTree> getAllTree(const std::vector<HyperNode> &allNode, const PatternGraph &p);
std::vector<HyperTree> getMinWidthTrees(const PatternGraph &p);
std::vector<HyperTree> addNode(HyperNode &tau, const HyperTree &tree, const PatternGraph &p);
std::vector<HyperNode> getCandidateNodes(const HyperTree &t, const std::vector<HyperNode> &allNodes, const PatternGraph &p);
static std::map<CanonType, double> canonToFW;
void fractionalWidth(HyperNode &tau, const PatternGraph &p);
std::vector<VertexID> getUncovered(const HyperTree &t, const PatternGraph &p);
bool isValid(const HyperTree &t, const PatternGraph &p);

// New functions for enumeration-based FHD computation
double computeFractionalHypertreeWidth(const PatternGraph &p, const std::vector<size_t> &eliminationOrder, FHD &fhd);
std::vector<HyperTree> buildOptimalFHDsByEnumeration(const PatternGraph &p, std::vector<std::vector<size_t>> &optimalOrders, double &minWidth);

#endif //GRAPH_MINING_SYSTEM_ALL_TD_H
