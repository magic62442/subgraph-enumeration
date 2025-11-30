//
// Created by anonymous authors on 2024/2/27.
//

#ifndef ASDMATCH_GRAPH_H
#define ASDMATCH_GRAPH_H


#include "config.h"
#include "utils.h"
#include "graph.h"
#include <map>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <queue>
#ifndef __APPLE__
#include <immintrin.h>
#include <x86intrin.h>
#endif

class LabeledGraph : public Graph {
private:
    ui _numLabels;
    LabelID *_labels;
    // undirected graph has _numEdges / 2 undirected edges.
    EdgeID *_reverseID;
    std::map<LabelID, ui> *_nlc;
    std::map<LabelID, EdgeID> *_labelOffsets; // vertex v's edges with label l starts from [offset[labelOffset[l]]
    std::map<LabelID, ui> _labelCount;
    ui *_offset2;   // offset of 'verticesByLabel'
    VertexID *_verticesByLabel; // vertices grouped by label

public:
    LabeledGraph();
    LabeledGraph(const Graph & baseGraph);
    ~LabeledGraph() override;

    const EdgeID *getOffsets() const {
        return _offsets;
    }

    const VertexID *getNbrs() const {
        return _nbrs;
    }

    const ui * getVerticesByLabel(const LabelID l, ui& count) const {
        count = _offset2[l + 1] - _offset2[l];
        return _verticesByLabel + _offset2[l];
    }

    VertexID *const getNeighbors(const VertexID v, ui& count) const {
        count = _offsets[v + 1] - _offsets[v];
        return _nbrs + _offsets[v];
    }

    LabelID getVertexLabel(const VertexID v) const {
        if (_labels != nullptr) return _labels[v];
        else return 0;
    }

    LabelID *getLabels() const {
        return _labels;
    }

    EdgeID getReverseID(const EdgeID e) const {
        return _reverseID[e];
    }

    EdgeID getEdgesByLabel(const VertexID v, const LabelID label, ui& count) const {
        count = _nlc[v][label];
        return _labelOffsets[v][label];
    }

    const std::map<LabelID, ui> &getNLC(const VertexID u) const {
        return _nlc[u];
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

    void setLabel(VertexID u, LabelID l) {
        _labels[u] = l;
    }

    void loadLabeldGraphFromText(const std::string& file);
    void buildReverseID();
    void buildNLC();
    EdgeID getEdgeID(VertexID v, VertexID w) const;
    EdgeID getUndirectedEID(VertexID v, VertexID w) const;
    void readFromStream(std::ifstream &inFile);
    void writeToStream(std::ofstream &outFile) const;
    void initWeights();
    void buildHyperGraph(HyperG &h) const;
    void buildFHD(FHD &fhd) const;
    void sortNeighborsUnlabeled() {
        for (ui i = 0; i < _numVertices; ++i) {
            std::sort(_nbrs + _offsets[i], _nbrs + _offsets[i + 1]);
        }
    }
    void generateRandomGraph(ui numVertices, ui numEdges);
    LabeledGraph generateRandomSubgraph(ui numVertices, ui numEdges) const;
    void writeToTextFile(const std::string &filename) const;
    void loadFromGraph(const Graph &baseGraph);
};

bool canAddVertex(const LabeledGraph& query, const std::vector<VertexID>& currentSet, VertexID newVertex);
void maxIndependentSet(const LabeledGraph& graph, const std::vector<VertexID>& subset, std::vector<VertexID>& maxSet);

#endif //ASDMATCH_GRAPH_H
