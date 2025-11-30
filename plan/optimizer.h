//
// Created by anonymous authors on 2024/5/22.
//

#ifndef ASDMATCH_OPTIMIZER_H
#define ASDMATCH_OPTIMIZER_H

#include "subset_structure.h"
#include "all_td.h"

void setTDExtention(HyperTree &t, const Graph &query, bool labeled);
std::vector<VertexID> globalCandidates(const Graph &query, HyperTree &t);
void optCostPlan(const PatternGraph &p, const Graph &g, CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount,
            VertexID **totalCandidates, ui *totalCandCount, std::vector<ui> &poses, std::vector<VertexID> &tmpCand,
                      HyperTree &t, PrefixNode *&pt);
void genRemainingOrder(const PatternGraph &p, CandidateSpace &cs, std::vector<std::vector<VertexID>> &remainingOrders,
        const std::vector<VertexID> &prefix, const HyperNode &bag);
double computeOrderCost(const PatternGraph &p, const Graph &g, HyperTree &t, VertexID nID, CandidateSpace &cs, bool *visited,
                        VertexID *partMatch, VertexID **candidates, ui *candCount, const std::vector<VertexID> &order,
                        VertexID **totalCandidates, ui *totalCandCount, std::vector<ui> &poses, std::vector<VertexID> &tmpCand,
                        std::vector<std::vector<VertexID>> &symmetryRules, std::vector<VertexID> &visitedBag,
                        std::vector<std::vector<std::vector<VertexID>>> &attributesBefore,
                        std::vector<std::vector<std::vector<VertexID>>> &smallerAttrs,
                        std::vector<std::vector<std::vector<VertexID>>> &largerAttrs,
                        std::vector<std::vector<int>> &candidatesBefore, std::vector<std::vector<VertexID>> &cartesianParent,
                        std::vector<PrefixNode *> &path, std::set<PrefixNode *> &computed);
void optCostOrder(const PatternGraph &p, const Graph &g, HyperTree &t, CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount,
                VertexID **totalCandidates, ui *totalCandCount, std::vector<ui> &poses, std::vector<VertexID> &tmpCand,
                std::vector<std::vector<VertexID>> &symmetryRules, PrefixNode *&bestPT, std::vector<std::vector<VertexID>> &bestOrders, double &minCost,
                int type = 0);
void optCostOrder(const PatternGraph &p, const Graph &g, HyperTree &t, CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount,
                  VertexID **totalCandidates, ui *totalCandCount, std::vector<ui> &poses, std::vector<VertexID> &tmpCand,
                  std::vector<std::vector<VertexID>> &symmetryRules, PrefixNode *&bestPT, std::vector<std::vector<VertexID>> &bestOrders, double &minCost,
                  double &intersectCost, double &materializeCost, int type = 0);
void buildTreesFromOrderFile(const std::string &filename, const Graph &query, CandidateSpace &cs,
                             HyperTree &t, PrefixNode *&pt);

#endif //ASDMATCH_OPTIMIZER_H
