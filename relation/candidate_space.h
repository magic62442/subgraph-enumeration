//
// Created by anonymous authors on 2024/2/29.
//

#ifndef ASDMATCH_CANDIDATE_SPACE_H
#define ASDMATCH_CANDIDATE_SPACE_H

#include "graph.h"
#include "compute_set_intersection.h"
#include <set>

class CandidateSpace {
private:
    const Graph &_query;
    const DataGraph &_data;

public:
    std::vector<std::vector<VertexID>> candidateSet;
    std::vector<std::vector<VertexID>> next; // Next vertex on the shortest path, used for cartesian product
    std::vector<std::vector<size_t>> dist;

    CandidateSpace(const Graph &query, const DataGraph &data);

    ~CandidateSpace() = default;

    bool buildCandCFL();
    void setQueryGraphWeights(Graph &query);
    size_t getDist(VertexID u1, VertexID u2) const;
    std::vector<VertexID> reconstructPath(VertexID i, VertexID j) const;
    ui getMaxSize() const;
    ui getMaxDegree() const;
    ui candSize(VertexID u) const;
    const VertexID * candSet(VertexID u) const;
    const VertexID * getNeighbors(VertexID u1, VertexID v1, VertexID u2, ui &count) const;
    const VertexID * getNeighbors(VertexID v1, ui &count) const;
    bool checkExists(VertexID v1, VertexID v2) const {
        return _data.getEdgeID(v1, v2) != -1;
    };
    std::vector<std::pair<VertexID, VertexID>> getEdges(int type) const;
};

#endif //ASDMATCH_CANDIDATE_SPACE_H
