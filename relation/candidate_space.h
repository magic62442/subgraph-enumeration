//
// Created by anonymous authors on 2024/2/29.
//

#ifndef ASDMATCH_CANDIDATE_SPACE_H
#define ASDMATCH_CANDIDATE_SPACE_H

#include "labeled_graph.h"
#include "compute_set_intersection.h"
#include <queue>
#include <set>

struct Dag {
    VertexID root;
    std::vector<std::vector<VertexID>> edges;
    std::vector<std::vector<VertexID>> children;
    std::vector<std::vector<VertexID>> parents;
    std::vector<VertexID> bfsSequence;

    Dag() = default;
    explicit Dag(const LabeledGraph &query);
    virtual ~Dag() = default;
    Dag &operator=(const Dag &rhs);
    void buildDAG(const LabeledGraph &query, const LabeledGraph &data);
    void selectRoot(const LabeledGraph &query, const LabeledGraph &data);
};

class CandidateSpace {
private:
    bool **_bitsetVertex;
    bool **_bitsetEdge;

    const LabeledGraph &_query;
    const LabeledGraph &_data;
    Dag _dag;

public:
    bool labeled;
    std::vector<std::vector<VertexID>> candidateSet;
    // index candidate edge by u1, v(a candidate of u1), u2
    std::vector<std::vector<std::vector<std::vector<VertexID>>>> candidateEdge;
    // all source shortest path for the query graph
    std::vector<std::vector<VertexID>> next; // Next vertex on the shortest path, used for cartesian product
    std::vector<std::vector<size_t>> dist;

    CandidateSpace(const LabeledGraph &query, const LabeledGraph &data, bool labeled);

    ~CandidateSpace() = default;

    bool buildCandVEQ();
    bool init();
    bool filter();
    bool filter(VertexID u, VertexID v);
    void construct();
    void buildCandidateEdge(bool checkEdge);
    bool edgeBipartiteSafety(VertexID u, VertexID v);
    bool buildCandCFL();
    VertexID selectCFLRoot();
    void buildCFLBFSTree(VertexID root, std::vector<std::vector<VertexID>> &levels, std::vector<VertexID> &order,
                         std::vector<std::vector<VertexID>> &before, std::vector<std::vector<VertexID>> &lowerLevelAfter,
                         std::vector<std::vector<VertexID>> &sameLevelAfter,
                         std::vector<std::vector<VertexID>> &children);
    void writeToStream(std::ofstream &outStream);
    void setQueryGraphWeights(LabeledGraph &query);
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

struct FilterVertex {
    double penalty;
    int filteredTime;
    VertexID u;

    FilterVertex(double penalty, int filteredTimes, VertexID u) : penalty(penalty), filteredTime(filteredTimes),
                                                                  u(u) {}

    bool operator<(const FilterVertex &rhs) const {
        return penalty > rhs.penalty;
    }
    bool operator>(const FilterVertex &rhs) const {
        return penalty < rhs.penalty;
    }
};

#endif //ASDMATCH_CANDIDATE_SPACE_H
