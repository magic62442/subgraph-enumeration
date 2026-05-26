//
// Created by anonymous authors on 2024/2/29.
//

#include "candidate_space.h"

CandidateSpace::CandidateSpace(const Graph &query, const DataGraph &data)
        : _query(query), _data(data) {
    candidateSet.resize(query.getNumVertices());
}

bool CandidateSpace::buildCandCFL() {
    candidateSet[0].resize(_data.getNumVertices());
    for (VertexID v = 0; v < _data.getNumVertices(); ++v) {
        candidateSet[0][v] = v;
    }
    return true;
}

void CandidateSpace::setQueryGraphWeights(Graph &query) {
    query.initWeights();
    const EdgeID *queryOffset = query.getOffsets();
    const VertexID *queryNbr = query.getNbors();
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        query.setVertexWeight(u, _data.getNumVertices());
        for (EdgeID queryEdge = queryOffset[u]; queryEdge < queryOffset[u + 1]; ++queryEdge) {
            query.setEdgeWeight(queryEdge, 2);
        }
    }
    query.allSourcesShortestPaths(dist, next);
}

size_t CandidateSpace::getDist(VertexID u1, VertexID u2) const {
    return dist[u1][u2];
}

std::vector<VertexID> CandidateSpace::reconstructPath(VertexID i, VertexID j) const {
    if (next[i][j] == std::numeric_limits<VertexID>::max()) return {};
    std::vector<VertexID> path = {i};
    while (i != j) {
        i = next[i][j];
        path.push_back(i);
    }
    return path;
}

ui CandidateSpace::getMaxSize() const {
    return _data.getNumVertices();
}

ui CandidateSpace::getMaxDegree() const {
    return _data.getMaxDegree();
}

ui CandidateSpace::candSize(VertexID) const {
    return _data.getNumVertices();
}

const VertexID * CandidateSpace::candSet(VertexID) const {
    return candidateSet[0].data();
}

const VertexID * CandidateSpace::getNeighbors(VertexID, VertexID v1, VertexID, ui &count) const {
    return _data.getNeighbors(v1, count);
}

const VertexID * CandidateSpace::getNeighbors(VertexID v1, ui &count) const {
    return _data.getNeighbors(v1, count);
}

std::vector<std::pair<VertexID, VertexID>> CandidateSpace::getEdges(int type) const {
    std::vector<std::pair<VertexID, VertexID>> edges;
    for (VertexID u = 0; u < _data.getNumVertices(); ++u) {
        ui degree;
        const VertexID *neighbors = _data.getNeighbors(u, degree);
        for (int i = 0; i < degree; ++i) {
            VertexID v = neighbors[i];
            if (type == 0) edges.emplace_back(u, v);
            else if (type == 1 && v > u) edges.emplace_back(u, v);
            else if (type == 2 && v < u) edges.emplace_back(u, v);
        }
    }
    return edges;
}
