//
// Created by anonymous author on 2022/7/19.
//

#ifndef SCOPE_GRAPH_H
#define SCOPE_GRAPH_H

#include "config.h"
#include "utils.h"
#include "clique.h"
extern "C" {
    #include "automorphism/nauty.h"
    #include "automorphism/gtools.h"
    #include "automorphism/naugroup.h"
};

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <map>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <functional>

class Graph {
protected:
    ui _numVertices;           // number of vertices in the graph
    ui _numEdges;              // number of directed edges in the graph.
                               // undirected graph has _numEdges / 2 undirected edges.
    EdgeID *_offsets;          // vertex v's edges are [offset[v], offset[v + 1])
    VertexID *_nbrs;           // incoming or outgoing neighbors
    EdgeID *_degree;
    ui *_nodeWeight;
    ui *_edgeWeight;
    ui _maxDegree;             // maximum degree in the graph

public:
    Graph();

    Graph(const Graph& rhs);

    Graph(ui numVertices, ui numEdges);

    Graph& operator=(const Graph& rhs);

    virtual ~Graph();

    inline ui getNumVertices() const {
        return _numVertices;
    }

    inline ui getNumEdges() const {
        return _numEdges;
    }

    inline VertexID *getNeighbors(VertexID v, ui &count) const {
        count = _offsets[v + 1] - _offsets[v];
        return _nbrs + _offsets[v];
    }

    inline ui getDegree(VertexID v) const {
        return _offsets[v + 1] - _offsets[v];
    }

    inline ui getMaxDegree() const {
        return _maxDegree;
    }

    ui getVertexWeight(const VertexID u) const {
        return _nodeWeight[u];
    }

    ui getEdgeWeight(const EdgeID e) const {
        return _edgeWeight[e];
    }

    void setVertexWeight(VertexID u, ui weight) {
        _nodeWeight[u] = weight;
    }

    void setEdgeWeight(EdgeID e, ui weight) {
        _edgeWeight[e] = weight;
    }

    void initWeights();

    int getEdgeID(VertexID v, VertexID w) const;
    int getUndirectedEID(VertexID v, VertexID w) const;
    void addDirectedEdges(const Edge *edgeList, ui num);

    EdgeID *getOffsets() const {
        return _offsets;
    }

    VertexID *getNbors() const {
        return _nbrs;
    }
    void allSourcesShortestPaths(std::vector<std::vector<size_t>> &dist, std::vector<std::vector<VertexID>> &next) const;
    void computeConnectedComponents(const std::vector<VertexID> &vertices, std::vector<std::vector<VertexID>> &components) const;
    bool isConnected(const std::vector<VertexID> &vertices) const;
    bool isClique() const {
        return _numEdges == _numVertices * (_numVertices - 1);
    };
    std::vector<VertexID> bfsForCC(VertexID start, std::set<VertexID> &visited, const std::vector<VertexID> &vertices) const;
    LabelID getVertexLabel(const VertexID v) const {
        return 0;
    }
};

class DataGraph : public Graph {
private:
    EdgeID *_largerOffsets;    // used for symmetry breaking, the offset of neighbors larger than itself

public:
    DataGraph();

    DataGraph(ui numVertices, ui numEdges);

    DataGraph& operator=(const DataGraph& rhs);

    ~DataGraph() override;
    inline VertexID *getNeighborsLargerID(VertexID u, ui &count) const {
        count = _offsets[u + 1] - _offsets[u] - _largerOffsets[u];
        return _nbrs + _offsets[u] + _largerOffsets[u];
    }
    inline VertexID *getNeighborsSmallerID(VertexID u, ui &count) const {
        count = _largerOffsets[u];
        return _nbrs + _offsets[u];
    }
    void loadDataGraph(const std::string &file);
    void saveBinaryGraph(const std::string &file) const;
    void loadBinaryGraph(const std::string &file);
    void buildLargerOffset();
    void initSpecialSparse(specialsparse *sg) const;
};

class PatternGraph : public Graph {
private:
    bool **_adjMatrix;
    int _orbitType;                 // 0: global counting. 1: vertex orbit. 2: edge orbit
    VertexID *_coreV;               // vertices with core number \geq 2
    VertexID *_peripheralV;         // other vertices
    ui _coreSize;
    CanonType _canonValue;          // the canonical value of the induced subgraph
    int _autoSize;                  // number of automorphisms
    int _divideFactor;              // _autoSize / auto group size of vertex 0 / edge (0,1)
    int *_v2o;                      // the canonical orbits of vertices in this node
    VertexID *_v2l;                 // the canonical labeling of vertices in this node
    std::vector<std::vector<std::vector<VertexID>>> _candidateRules;

public:
    PatternGraph();

    PatternGraph(ui numVertices, ui numEdges, const Edge *edgeList);

    PatternGraph(const PatternGraph & rhs);

    PatternGraph& operator=(const PatternGraph &rhs);

    ~PatternGraph() override;

    inline bool isEdge(VertexID u1, VertexID u2) const {
        return _adjMatrix[u1][u2];
    }

    inline VertexID *getCoreV(ui &coreSize) const {
        coreSize = _coreSize;
        return _coreV;
    }

    inline VertexID *getPeripheralV() const {
        return _peripheralV;
    }

    inline VertexID *getV2L() const {
        return _v2l;
    }

    inline int getOrbitType() const {
        return _orbitType;
    }

    inline int getOrbit(VertexID u) const {
        return _v2o[u];
    }

    inline int getAutoSize() const {
        return _autoSize;
    }

    inline int getDivideFactor() const {
        return _divideFactor;
    }

    inline CanonType getCanonValue() const {
        return _canonValue;
    }

    const std::vector<std::vector<std::vector<VertexID>>> &getCandidateRules() const {
        return _candidateRules;
    }

    void loadPatternGraph(const std::string &file);
    void addEdgeList(const Edge *edgeList, ui num);
    void computeAutoGroup();
    void computeCandidateRules();
    void buildCorePeripheral();
    Edge *coreUndirectedEdges(ui &num) const;
    void printGraph(bool direction = true) const;
    bool isClique() const;
};

CanonType subgraphCanonValue(const PatternGraph &p, const std::vector<VertexID> &v, int *v2o = nullptr);

#endif //SCOPE_GRAPH_H
