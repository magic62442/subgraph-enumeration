//
// Created by anonymous authors on 2024/5/22.
//

#ifndef ASDMATCH_OPTIMIZER_H
#define ASDMATCH_OPTIMIZER_H

#include "estimator.h"
#include "all_td.h"

void setTDExtention(HyperTree &t, const Graph &query);
std::vector<VertexID> globalCandidates(const Graph &query, HyperTree &t);
void optCostPlan(const PatternGraph &p, const Graph &g, CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount,
            VertexID **totalCandidates, ui *totalCandCount, std::vector<ui> &poses, std::vector<VertexID> &tmpCand,
                      HyperTree &t, PrefixNode *&pt);
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

#endif //ASDMATCH_OPTIMIZER_H
