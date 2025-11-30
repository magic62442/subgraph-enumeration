#ifndef PREPROCESSING_HZY_H
#define PREPROCESSING_HZY_H

#include "hypergraph.h"

bool InducedHyperG(HyperG& H, std::vector<size_t>& remain_vertex, std::map<size_t, size_t>& Vres_map);
bool DelDegree1Vertex(HyperG& H, Order& prefix_o, std::map<size_t, size_t>& Vres_map);
bool DelCoveredEdge(HyperG& H);
bool DelISOVertex(HyperG& H, Order& prefix_o, std::map<size_t, size_t>& Vres_map);
void Preprocessing(HyperG& H, Order& prefix_o, std::map<size_t, size_t>& Vres_map);

#endif // PREPROCESSING_HZY_H
