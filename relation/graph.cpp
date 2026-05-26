//
// Created by anonymous author on 2022/7/19.
//

#include "graph.h"

Graph::Graph() {
    _numVertices = 0;
    _numEdges = 0;
    _offsets = nullptr;
    _nbrs = nullptr;
    _degree = nullptr;
    _nodeWeight = nullptr;
    _edgeWeight = nullptr;
    _maxDegree = 0;
}

Graph::Graph(ui numVertices, ui numEdges) : _numVertices(numVertices), _numEdges(numEdges) {
    _offsets = new VertexID[_numVertices + 1];
    _nbrs = new VertexID[_numEdges];
    _degree = new VertexID[_numVertices];
    memset(_degree, 0, sizeof(VertexID) * _numVertices);
    _nodeWeight = nullptr;
    _edgeWeight = nullptr;
    _maxDegree = 0;
}

Graph::Graph(const Graph& rhs) {
    _numVertices = rhs._numVertices;
    _numEdges = rhs._numEdges;
    _maxDegree = rhs._maxDegree;
    if (rhs._offsets == nullptr) _offsets = nullptr;
    else {
        _offsets = new EdgeID[_numVertices + 1];
        memcpy(_offsets, rhs._offsets, sizeof(EdgeID) * (_numVertices + 1));
    }
    if (rhs._nbrs == nullptr) _nbrs = nullptr;
    else {
        _nbrs = new VertexID[_numEdges];
        memcpy(_nbrs, rhs._nbrs, sizeof(VertexID) * _numEdges);
    }
    if (rhs._degree == nullptr) _degree = nullptr;
    else {
        _degree = new VertexID[_numVertices];
        memcpy(_degree, rhs._degree, sizeof(VertexID) * _numVertices);
    }
    _nodeWeight = nullptr;
    _edgeWeight = nullptr;
}

Graph &Graph::operator=(const Graph &rhs) {
    if (this == &rhs)
        return *this;
    _numVertices = rhs._numVertices;
    _numEdges = rhs._numEdges;
    _maxDegree = rhs._maxDegree;
    if (rhs._offsets == nullptr) _offsets = nullptr;
    else {
        delete[] _offsets;
        _offsets = new EdgeID[_numVertices + 1];
        memcpy(_offsets, rhs._offsets, sizeof(EdgeID) * (_numVertices + 1));
    }
    if (rhs._nbrs == nullptr) _nbrs = nullptr;
    else {
        delete[] _nbrs;
        _nbrs = new VertexID[_numEdges];
        memcpy(_nbrs, rhs._nbrs, sizeof(VertexID) * _numEdges);
    }
    if (rhs._degree == nullptr) _degree = nullptr;
    else {
        delete[] _degree;
        _degree = new VertexID[_numVertices];
        memcpy(_degree, rhs._degree, sizeof(VertexID) * _numVertices);
    }
    if (rhs._nodeWeight == nullptr) _nodeWeight = nullptr;
    else {
        delete[] _nodeWeight;
        _nodeWeight = new ui[_numVertices];
        memcpy(_nodeWeight, rhs._nodeWeight, sizeof(ui) * _numVertices);
    }
    if (rhs._edgeWeight == nullptr) _edgeWeight = nullptr;
    else {
        delete[] _edgeWeight;
        _edgeWeight = new ui[_numEdges];
        memcpy(_edgeWeight, rhs._edgeWeight, sizeof(ui) * _numEdges);
    }
    return *this;
}

Graph::~Graph() {
    delete[] _offsets;
    delete[] _nbrs;
    delete[] _degree;
    delete[] _nodeWeight;
    delete[] _edgeWeight;
}

void Graph::initWeights() {
    delete[] _nodeWeight;
    delete[] _edgeWeight;
    _nodeWeight = new ui[_numVertices];
    _edgeWeight = new ui[_numEdges];
    std::fill(_nodeWeight, _nodeWeight + _numVertices, 1);
    std::fill(_edgeWeight, _edgeWeight + _numEdges, 1);
}

int Graph::getEdgeID(VertexID v, VertexID w) const {
    if (v > w) std::swap(v, w);
    int low = (int)_offsets[v];
    int high = (int)_offsets[v + 1] - 1;
    int mid;

    while (high - low >= 16) {
        mid = low + ((high - low) >> 1);
#ifndef __APPLE__
        _mm_prefetch((char *) &_nbrs[(mid + 1 + high) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &_nbrs[(mid - 1 + low) / 2], _MM_HINT_T0);
#endif
        if (_nbrs[mid] == w) return mid;
        if (_nbrs[mid] > w) high = mid - 1;
        else low = mid + 1;
    }

    for (int i = low; i <= high; ++i) {
        if (_nbrs[i] == w) return i;
    }

    return -1;
}

int Graph::getUndirectedEID(VertexID v, VertexID w) const {
    int low = (int)_offsets[v];
    int high = (int)_offsets[v + 1] - 1;
    int mid;

    while (high - low >= 16) {
        mid = low + ((high - low) >> 1);
#ifndef __APPLE__
        _mm_prefetch((char *) &_nbrs[(mid + 1 + high) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &_nbrs[(mid - 1 + low) / 2], _MM_HINT_T0);
#endif
        if (_nbrs[mid] == w) return mid;
        if (_nbrs[mid] > w) high = mid - 1;
        else low = mid + 1;
    }

    for (int i = low; i <= high; ++i) {
        if (_nbrs[i] == w) return i;
    }

    return -1;
}

void Graph::addDirectedEdges(const Edge *edgeList, ui num) {
    // build CSR graph from the edge list
    delete[] _offsets;
    _offsets = new EdgeID[_numVertices + 1];
    _offsets[0] = 0;
    delete[] _nbrs;
    _nbrs = new VertexID[num];
    delete[] _degree;
    _degree = new VertexID[_numVertices];
    memset(_degree, 0, sizeof(VertexID) * _numVertices);
    VertexID *nborsOffset = new VertexID[_numVertices + 1];
    memset(nborsOffset, 0, sizeof(VertexID) * (_numVertices + 1));
    for (int i = 0; i < num; ++i) {
        VertexID u1 = edgeList[i].first;
        ++_degree[u1];
    }
    for (int i = 0; i < _numVertices; ++i) {
        _offsets[i + 1] = _offsets[i] + _degree[i];
    }
    for (int i = 0; i < num; ++i) {
        VertexID u1 = edgeList[i].first, u2 = edgeList[i].second;
        EdgeID offset1 = _offsets[u1] + nborsOffset[u1];
        _nbrs[offset1] = u2;
        ++nborsOffset[u1];
    }
    // sort neighbors in ascending order of vertex id
    for (int i = 0; i < _numVertices; ++i)
        std::sort(_nbrs + _offsets[i], _nbrs + _offsets[i + 1]);

    delete[] nborsOffset;
}

void Graph::allSourcesShortestPaths(std::vector<std::vector<size_t>> &dist, std::vector<std::vector<VertexID>> &next) const {
    ui n = _numVertices;
    dist = std::vector<std::vector<size_t>>(n, std::vector<size_t>(n, std::numeric_limits<size_t>::max()));
    next = std::vector<std::vector<VertexID>>(n, std::vector<VertexID>(n, std::numeric_limits<VertexID>::max()));

    for (ui i = 0; i < n; ++i) {
        dist[i][i] = std::numeric_limits<size_t>::max();
        next[i][i] = i; // Self-loops are initialized to the vertex itself
    }

    // Initialize distances and paths for direct edges
    for (VertexID v = 0; v < n; ++v) {
        for (EdgeID e = _offsets[v]; e < _offsets[v + 1]; ++e) {
            VertexID u = _nbrs[e];
            dist[v][u] = getEdgeWeight(e);
            next[v][u] = u;
        }
    }

    // Floyd-Warshall algorithm to update distances and paths
    for (ui k = 0; k < n; ++k) {
        for (ui i = 0; i < n; ++i) {
            for (ui j = 0; j < n; ++j) {
                if (dist[i][k] == std::numeric_limits<size_t>::max() || dist[k][j] == std::numeric_limits<size_t>::max()) continue;
                if (dist[i][k] * dist[k][j] < dist[i][j]) {
                    dist[i][j] = dist[i][k] * dist[k][j];
                    next[i][j] = next[i][k];
                }
            }
        }
    }
}

void Graph::computeConnectedComponents(const std::vector<VertexID> &vertices, std::vector<std::vector<VertexID>> &components) const {
    std::set<VertexID> visited;
    for (VertexID v : vertices) {
        if (visited.find(v) == visited.end()) {
            std::vector<VertexID> component = bfsForCC(v, visited, vertices);
            components.push_back(component);
        }
    }
}

bool Graph::isConnected(const std::vector<VertexID> &vertices) const {
    if (isClique()) return true;
    std::set<VertexID> visited;
    return bfsForCC(vertices[0], visited, vertices).size() == vertices.size();
}

std::vector<VertexID> Graph::bfsForCC(VertexID start, std::set<VertexID> &visited, const std::vector<VertexID>& vertices) const {
    std::vector<VertexID> component;
    std::queue<VertexID> queue;
    queue.push(start);
    while (!queue.empty()) {
        VertexID current = queue.front();
        queue.pop();

        if (visited.find(current) == visited.end()) {
            visited.insert(current);
            component.push_back(current);

            // Iterate through all neighbors of current vertex
            for (EdgeID i = _offsets[current]; i < _offsets[current + 1]; ++i) {
                VertexID nbr = _nbrs[i];
                // Check if this neighbor is part of the subgraph we are interested in
                if (visited.find(nbr) == visited.end() && std::find(vertices.begin(), vertices.end(), nbr) != vertices.end()) {
                    queue.push(nbr);
                }
            }
        }
        for (VertexID u: vertices)
            if (visited.find(u) == visited.end()) {
                queue.push(u);
                break;
            }
    }

    return component;
}

void DataGraph::loadDataGraph(const std::string &file) {
    // Check if file is binary by reading magic number
    std::ifstream checkFile(file, std::ios::binary);
    if (!checkFile.is_open()) {
        std::cout << "Can not open file " << file << "!" << std::endl;
        exit(1);
    }

    char magic[8] = {0};
    checkFile.read(magic, 7);
    checkFile.close();

    // If magic number matches "BINGRPH", load as binary
    if (std::string(magic) == "BINGRPH") {
        loadBinaryGraph(file);
        return;
    }

    // Otherwise load as text file
    std::ifstream inFile(file);
    // get the number of vertices and edges from the graph file
    inFile >> _numVertices >> _numEdges;
    _numEdges *= 2;
    _offsets = new EdgeID[_numVertices + 1];
    _offsets[0] = 0;
    _nbrs = new VertexID[_numEdges];
    _degree = new VertexID[_numVertices];
    memset(_degree, 0, sizeof(VertexID) * _numVertices);
    // read edge list
    EdgeID edgeCnt = 0;
    Edge *edgeList = new Edge[_numEdges];
    for (ui i = 0; i < _numEdges / 2; i++) {
        VertexID _v1, _v2;
        if (inFile >> _v1 >> _v2) {
            edgeList[edgeCnt] = std::make_pair(_v1, _v2);
            ++edgeCnt;
            edgeList[edgeCnt] = std::make_pair(_v2, _v1);
            ++edgeCnt;
        }
        else {
            printf("error: expecting %d edges, getting %d edges.\n", _numEdges / 2, edgeCnt / 2);
            exit(1);
        }
    }
    inFile.close();
    addDirectedEdges(edgeList, _numEdges);

    // compute maximum degree
    _maxDegree = 0;
    for (VertexID i = 0; i < _numVertices; ++i) {
        if (_degree[i] > _maxDegree) {
            _maxDegree = _degree[i];
        }
    }

    // build larger offset
    buildLargerOffset();

    delete[] edgeList;
}

void DataGraph::saveBinaryGraph(const std::string &file) const {
    std::ofstream outFile(file, std::ios::binary);
    if (!outFile.is_open()) {
        std::cout << "Can not create file " << file << "!" << std::endl;
        exit(1);
    }

    // Write binary format marker
    const char magic[] = "BINGRPH";
    outFile.write(magic, 7);

    // Write number of vertices and edges (original _numEdges is doubled)
    ui originalNumEdges = _numEdges / 2;
    outFile.write(reinterpret_cast<const char*>(&_numVertices), sizeof(ui));
    outFile.write(reinterpret_cast<const char*>(&originalNumEdges), sizeof(ui));

    // Reconstruct original edge list from graph structure
    std::vector<Edge> edgeList;
    edgeList.reserve(originalNumEdges);

    for (VertexID u = 0; u < _numVertices; ++u) {
        ui count;
        VertexID *neighbors = getNeighbors(u, count);
        for (ui i = 0; i < count; ++i) {
            VertexID v = neighbors[i];
            // Only add each undirected edge once (avoid duplicates)
            if (u < v) {
                edgeList.push_back(std::make_pair(u, v));
            }
        }
    }

    // Write edges
    for (const auto &edge : edgeList) {
        outFile.write(reinterpret_cast<const char*>(&edge.first), sizeof(VertexID));
        outFile.write(reinterpret_cast<const char*>(&edge.second), sizeof(VertexID));
    }

    outFile.close();
}

void DataGraph::loadBinaryGraph(const std::string &file) {
    std::ifstream inFile(file, std::ios::binary);
    if (!inFile.is_open()) {
        std::cout << "Can not open file " << file << "!" << std::endl;
        exit(1);
    }

    // Read and verify binary format marker
    char magic[8] = {0};
    inFile.read(magic, 7);
    if (std::string(magic) != "BINGRPH") {
        std::cout << "Invalid binary graph file format!" << std::endl;
        exit(1);
    }

    // Read number of vertices and edges
    ui originalNumEdges;
    inFile.read(reinterpret_cast<char*>(&_numVertices), sizeof(ui));
    inFile.read(reinterpret_cast<char*>(&originalNumEdges), sizeof(ui));

    _numEdges = originalNumEdges * 2;  // Double for directed representation

    // Allocate memory
    _offsets = new EdgeID[_numVertices + 1];
    _offsets[0] = 0;
    _nbrs = new VertexID[_numEdges];
    _degree = new VertexID[_numVertices];
    memset(_degree, 0, sizeof(VertexID) * _numVertices);

    // Read edge list
    Edge *edgeList = new Edge[_numEdges];
    EdgeID edgeCnt = 0;

    for (ui i = 0; i < originalNumEdges; i++) {
        VertexID u, v;
        inFile.read(reinterpret_cast<char*>(&u), sizeof(VertexID));
        inFile.read(reinterpret_cast<char*>(&v), sizeof(VertexID));

        // Add both directions for undirected graph
        edgeList[edgeCnt] = std::make_pair(u, v);
        ++edgeCnt;
        edgeList[edgeCnt] = std::make_pair(v, u);
        ++edgeCnt;
    }

    inFile.close();

    // Build graph structure using existing function
    addDirectedEdges(edgeList, _numEdges);

    // Compute maximum degree
    _maxDegree = 0;
    for (VertexID i = 0; i < _numVertices; ++i) {
        if (_degree[i] > _maxDegree) {
            _maxDegree = _degree[i];
        }
    }

    // Build larger offset
    buildLargerOffset();

    delete[] edgeList;
}

DataGraph::DataGraph() {
    _largerOffsets = nullptr;
}

DataGraph::~DataGraph() {
    delete[] _largerOffsets;
}

DataGraph::DataGraph(ui numVertices, ui numEdges) : Graph(numVertices, numEdges) {
    _largerOffsets = nullptr;
}

DataGraph& DataGraph::operator=(const DataGraph& rhs) {
    if (this == &rhs)
        return *this;
    Graph::operator=(rhs);
    delete[] _largerOffsets;
    _largerOffsets = new EdgeID[_numVertices + 1];
    memcpy(_largerOffsets, rhs._offsets, sizeof(EdgeID) * (_numVertices + 1));

    return *this;
}

void DataGraph::buildLargerOffset() {
    _largerOffsets = new EdgeID[_numVertices];
    for (VertexID i = 0; i < _numVertices; ++i) {
        ui temp;
        VertexID *neighbors = getNeighbors(i, temp);
        int begin = 0, end = (int)temp - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] > i) {
                end = mid - 1;
            }
            else {
                begin = mid + 1;
            }
        }
        _largerOffsets[i] = (EdgeID)begin;
    }
}

void DataGraph::initSpecialSparse(specialsparse *sg) const {
    sg->n = _numVertices;
    sg->e = _numEdges;
    sg->edges = (edge *)malloc(sg->e * sizeof(edge));
    for (VertexID u = 0; u < _numVertices; ++u) {
        for (EdgeID e = _offsets[u]; e < _offsets[u + 1]; ++e) {
            VertexID w = _nbrs[e];
            sg->edges[e].s = u;
            sg->edges[e].t = w;
        }
    }
}

PatternGraph::PatternGraph() {
    _adjMatrix = nullptr;
    _orbitType = 0;
    _coreV = nullptr;
    _peripheralV = nullptr;
    _coreSize = 0;
    _canonValue = 0;
    _autoSize = 0;
    _divideFactor = 1;
    _v2o = nullptr;
    _v2l = nullptr;
}

PatternGraph::PatternGraph(ui numVertices, ui numEdges, const Edge *edgeList) : Graph(
        numVertices, numEdges) {
    _adjMatrix = new bool*[_numVertices];
    _degree = new VertexID[_numVertices];
    memset(_degree, 0, sizeof(VertexID) * _numVertices);
    for (int i = 0; i < _numVertices; ++i) {
        _adjMatrix[i] = new bool[_numVertices];
        memset(_adjMatrix[i], false, sizeof(bool) * (_numVertices));
    }
    _orbitType = 0;
    addEdgeList(edgeList, _numEdges);
    _coreV = nullptr;
    _peripheralV = nullptr;
    _coreSize = 0;
    _canonValue = 0;
    _autoSize = 0;
    _divideFactor = 1;
    _v2o = nullptr;
    _v2l = nullptr;
}

PatternGraph::PatternGraph(const PatternGraph & rhs) : Graph(rhs) {
    _orbitType = rhs._orbitType;
    _coreSize = rhs._coreSize;
    _canonValue = rhs._canonValue;
    _autoSize = rhs._autoSize;
    _divideFactor = rhs._divideFactor;
    _candidateRules = rhs._candidateRules;
    if (rhs._adjMatrix == nullptr)
        _adjMatrix = nullptr;
    else {
        _adjMatrix = new bool*[_numVertices];
        for (int i = 0; i < _numVertices; ++i) {
            _adjMatrix[i] = new bool[_numVertices];
            memcpy(_adjMatrix[i], rhs._adjMatrix[i], sizeof(bool) * (_numVertices));
        }
    }
    if (rhs._coreV == nullptr)
        _coreV = nullptr;
    else {
        _coreV = new VertexID[_coreSize];
        memcpy(_coreV, rhs._coreV, sizeof(VertexID) * _coreSize);
    }
    if (rhs._peripheralV == nullptr)
        _peripheralV = nullptr;
    else {
        ui peripheralSize = _numVertices - _coreSize;
        _peripheralV = new VertexID[peripheralSize];
        memcpy(_peripheralV, rhs._peripheralV, sizeof(VertexID) * peripheralSize);
    }
    if (rhs._v2o == nullptr)
        _v2o = nullptr;
    else {
        _v2o = new int[_numVertices];
        memcpy(_v2o, rhs._v2o, sizeof(int) * _numVertices);
    }
    if (rhs._v2l == nullptr)
        _v2l = nullptr;
    else {
        _v2l = new VertexID[_numVertices];
        memcpy(_v2l, rhs._v2l, sizeof(int) * _numVertices);
    }
}

PatternGraph& PatternGraph::operator=(const PatternGraph &rhs) {
    if (this == &rhs)
        return *this;
    Graph::operator=(rhs);
    _orbitType = rhs._orbitType;
    _coreSize = rhs._coreSize;
    _canonValue = rhs._canonValue;
    _autoSize = rhs._autoSize;
    _divideFactor = rhs._divideFactor;
    _candidateRules = rhs._candidateRules;
    if (rhs._adjMatrix == nullptr)
        _adjMatrix = nullptr;
    else {
        if (_adjMatrix) {
            for (int i = 0; i < _numVertices; ++i) {
                delete[] _adjMatrix[i];
            }
            delete[] _adjMatrix;
        }
        _adjMatrix = new bool*[_numVertices];
        for (int i = 0; i < _numVertices; ++i) {
            _adjMatrix[i] = new bool[_numVertices];
            memcpy(_adjMatrix[i], rhs._adjMatrix[i], sizeof(bool) * (_numVertices));
        }
    }
    if (rhs._coreV == nullptr)
        _coreV = nullptr;
    else {
        delete _coreV;
        _coreV = new VertexID[_coreSize];
        memcpy(_coreV, rhs._coreV, sizeof(VertexID) * _coreSize);
    }
    if (rhs._peripheralV == nullptr)
        _peripheralV = nullptr;
    else {
        delete _peripheralV;
        ui peripheralSize = _numVertices - _coreSize;
        _peripheralV = new VertexID[peripheralSize];
        memcpy(_peripheralV, rhs._peripheralV, sizeof(VertexID) * peripheralSize);
    }
    if (rhs._v2o == nullptr)
        _v2o = nullptr;
    else {
        delete[] _v2o;
        _v2o = new int[_numVertices];
        memcpy(_v2o, rhs._v2o, sizeof(int) * _numVertices);
    }
    if (rhs._v2l == nullptr)
        _v2l = nullptr;
    else {
        delete[] _v2l;
        _v2l = new VertexID[_numVertices];
        memcpy(_v2l, rhs._v2l, sizeof(int) * _numVertices);
    }
    return *this;
}

PatternGraph::~PatternGraph() {
    for (int i = 0; i < _numVertices; ++i) {
        delete[] _adjMatrix[i];
    }
    delete[] _adjMatrix;
    delete[] _coreV;
    delete[] _peripheralV;
    delete[] _v2o;
    delete[] _v2l;
}

void PatternGraph::loadPatternGraph(const std::string &file) {
    std::ifstream inFile(file);
    if (!inFile.is_open()) {
        std::cout << "Can not open file " << file << "!" << std::endl;
        exit(1);
    }
    // get the number of vertices and edges from the graph file
    inFile >> _numVertices >> _numEdges;
    _numEdges *= 2;
    // read edge list
    EdgeID edgeCnt = 0;
    Edge *edgeList = new Edge[_numEdges];
    for (ui i = 0; i < _numEdges / 2; i++) {
        VertexID u1, u2;
        if (inFile >> u1 >> u2) {
            edgeList[edgeCnt] = std::make_pair(u1, u2);
            ++edgeCnt;
            edgeList[edgeCnt] = std::make_pair(u2, u1);
            ++edgeCnt;
        }
        else {
            printf("error: expecting %d edges, getting %d edges.\n", _numEdges / 2, edgeCnt / 2);
            exit(1);
        }
    }

    // read the orbit
    VertexID *mapping = new VertexID[_numVertices];
    for (int i = 0; i < _numVertices; ++i) {
        mapping[i] = i;
    }
    VertexID vOrbit, eOrbit1, eOrbit2;
    if (inFile >> _orbitType) {
        // vertex orbit
        if (_orbitType == 1) {
            if (inFile >> vOrbit) {
                // relabel vertices: let the orbit vertex be 0.
                for (int i = 0; i < _numVertices; ++i) {
                    if (i == 0) mapping[i] = vOrbit;
                    else if (i == vOrbit) mapping[i] = 0;
                    else mapping[i] = i;
                }
            }
            else {
                printf("error: expecting a vertex orbit.\n");
                exit(1);
            }
        }
            // edge orbit
        else if (_orbitType == 2) {
            if (inFile >> eOrbit1 >> eOrbit2) {
                // relabel vertices: let the orbit edge be 0, 1.
                for (int i = 0; i < _numVertices; ++i) {
                    if (i == eOrbit1) mapping[i] = 0;
                    else if (i == eOrbit2) mapping[i] = 1;
                }
                int i = 2, j = 0;
                while (i < _numVertices && j < _numVertices) {
                    if (j != eOrbit1 && j != eOrbit2) {
                        mapping[j] = i;
                        ++i;
                    }
                    ++j;
                }
            }
            else {
                printf("error: expecting an edge orbit.\n");
                exit(1);
            }
        }
        else {
            printf("Scope currently supports vertex or edge orbit.\n");
            exit(1);
        }
    }
    // else, the file does not give the orbit, so we do the global counting.
    // relabel the edge list according to mapping
    for (ui i = 0; i < _numEdges; ++i) {
        VertexID u1 = edgeList[i].first, u2 = edgeList[i].second;
        edgeList[i].first = mapping[u1];
        edgeList[i].second = mapping[u2];
    }
    addEdgeList(edgeList,  _numEdges);
    if (_orbitType == 2 && !(isEdge(mapping[eOrbit1], mapping[eOrbit2]))) {
        printf("error: expecting an edge orbit.\n");
        exit(1);
    }
    computeAutoGroup();
    computeCandidateRules();
    buildCorePeripheral();
    inFile.close();
    delete[] mapping;
    delete[] edgeList;
}

// should call the constructor with numVertices and numEdges before using this function.
// add a set of directed edges to an empty pattern graph.
void PatternGraph::addEdgeList(const Edge *edgeList, ui num) {
    addDirectedEdges(edgeList, num);
    _adjMatrix = new bool*[_numVertices];
    for (int i = 0; i < _numVertices; ++i) {
        _adjMatrix[i] = new bool[_numVertices];
        memset(_adjMatrix[i], false, sizeof(bool) * (_numVertices));
    }
    for (int i = 0; i < num; ++i) {
        VertexID u1 = edgeList[i].first, u2 = edgeList[i].second;
        _adjMatrix[u1][u2] = true;
    }
}

/*
 * Compute canonical orbits, canonical value and automorphism group size
 * The canonical value corresponds to the graph colored by the query orbit
 * */
void PatternGraph::computeAutoGroup() {
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(graph,cg,cg_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    options.digraph = TRUE;
    options.defaultptn = FALSE;
    statsblk stats;
    int n = int(_numVertices);
    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC2(graph,cg,cg_sz,m,n,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    EMPTYGRAPH(g,m,n);
    // treat pattern as directed
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (isEdge(i, j)) {
                        ADDONEARC(g, i, j, m);
            }
        }
    }
    for (int i = 0; i < n; ++i)
        lab[i] = i;
    for (int i = 0; i < n - 1; ++i)
        ptn[i] = 1;
    ptn[n - 1] = 0;
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
    std::map<VertexID, int> canonV2O;
    delete _v2l;
    _v2l = new VertexID[_numVertices];
    for (int i = 0; i < _numVertices; ++i) {
        _v2l[lab[i]] = i;
    }
    // group canonical vertices according to orbits and reorder groups according to the min canonical vertex id
    int o = 0;
    std::vector<std::vector<int>> canonGroup(_numVertices);
    std::vector<int> canonV2Group(_numVertices);
    std::vector<int> group2O(_numVertices, -1);
    for (int i = 0; i < _numVertices; ++i) {
        canonGroup[orbits[i]].push_back(_v2l[i]);
        canonV2Group[_v2l[i]] = orbits[i];
    }
    for (int i = 0; i < _numVertices; ++i) {
        int group = canonV2Group[i];
        if (group2O[group] == -1) {
            group2O[group] = o;
            ++o;
        }
    }
    for (int i = 0; i < _numVertices; ++i) {
        for (int u: canonGroup[i]) {
            canonV2O[u] = group2O[i];
        }
    }
    delete _v2o;
    _v2o = new int[_numVertices];
    for (int i = 0; i < _numVertices; ++i) {
        _v2o[i] = canonV2O[_v2l[i]];
    }
    _autoSize = (int)round(stats.grpsize1) * quickPow10(stats.grpsize2);
    int autoGroupSize = 0;
    if (_orbitType == 0) _divideFactor = _autoSize;
    else if (_orbitType == 1) {
        o = _v2o[0];
        for (VertexID i = 0; i < _numVertices; ++i) {
            if (_v2o[i] == o)
                ++autoGroupSize;
        }
        _divideFactor = _autoSize / autoGroupSize;
    }
    else if (_orbitType == 2) {
        int orbit0 = _v2o[0], orbit1 = _v2o[1];
        for (int u1 = 0; u1 < _numVertices; ++u1) {
            for (int u2 = u1 + 1; u2 < _numVertices; ++u2) {
                if (isEdge(u1, u2) && ((_v2o[u1] == orbit0 && _v2o[u2] == orbit1) ||
                                       (_v2o[u2] == orbit0 && _v2o[u1] == orbit1))) {
                    ++autoGroupSize;
                }
            }
        }
        _divideFactor = _autoSize / autoGroupSize;
    }
    for (int i = 0; i < n; ++i)
        lab[i] = i;
    for (int i = 0; i < n - 1; ++i)
        ptn[i] = 1;
    ptn[n - 1] = 0;
    if (_orbitType == 1)
        ptn[0] = 0;
    else if (_orbitType == 2) {
        ptn[0] = 0;
        ptn[1] = 0;
    }
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
    // save the canonical value of the pattern
    int pos = 0;
    _canonValue = 0;
    for (int i = 0; i < m * n; ++i) {
        unsigned long row = cg[i];
        bool current;
        for (unsigned long j = 0; j < n; ++j) {
            if (i != j) {
                if (1 & row >> (sizeof(unsigned long) * 8 - j - 1)) current = true;
                else current = false;
                if (current)
                    _canonValue += 1 << pos;
                ++pos;
            }
        }
    }
}

void PatternGraph::computeCandidateRules() {
    _candidateRules.clear();
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    options.digraph = TRUE;
    options.defaultptn = FALSE;
    statsblk stats;
    int n = _numVertices;
    // m should be 1
    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    EMPTYGRAPH(g,m,n);
    for (VertexID i = 0; i < n; ++i) {
        for (VertexID j = 0; j < n; ++j) {
            if (isEdge(i, j))
                        ADDONEARC(g, i, j, m);
        }
    }
    std::vector<int> initialV;
    std::queue<std::vector<int>> coloredQ;
    std::queue<std::vector<std::vector<VertexID>>> groupQ;
    groupQ.push(std::vector<std::vector<VertexID>>(0));
    coloredQ.push(initialV);
    while (!coloredQ.empty()) {
        std::vector<int> coloredV = coloredQ.front();
        std::vector<std::vector<VertexID>> currentGroup = groupQ.front();
        coloredQ.pop();
        groupQ.pop();
        // reorganize lab and ptn, put colored vertices before uncolored vertices in lab
        for (int i = 0; i < coloredV.size(); ++i) {
            lab[i] = coloredV[i];
            ptn[i] = 0;
        }
        int pos = (int)coloredV.size();
        for (int i = 0; i < n; ++i) {
            if (std::find(coloredV.begin(), coloredV.end(), i) == coloredV.end()) {
                lab[pos] = i;
                ptn[pos] = 1;
                ++pos;
            }
        }
        ptn[n - 1] = 0;
        densenauty(g, lab, ptn, orbits, &options, &stats, m, n, NULL);
        int autoSize = (int)round(stats.grpsize1) * quickPow10(stats.grpsize2);
        if (autoSize != 1) {
            std::vector<std::vector<VertexID>> o2v = std::vector<std::vector<VertexID>>(n);
            for (int i = 0; i < n; ++i)
                o2v[orbits[i]].push_back(i);
            for (const auto &group: o2v) {
                if (group.size() > 1) {
                    // group is a valid rule
                    if (group.size() > 2) {
                        for (int u: group) {
                            std::vector<int> newColoredV = coloredV;
                            std::vector<std::vector<VertexID>> newGroup = currentGroup;
                            newColoredV.push_back(u);
                            std::vector<VertexID> rule;
                            rule.push_back(u);
                            for (int w: group) {
                                if (w != u) rule.push_back(w);
                            }
                            newGroup.push_back(rule);
                            coloredQ.push(newColoredV);
                            groupQ.push(newGroup);
                        }
                    }
                    else {
                        std::vector<int> newColoredV = coloredV;
                        std::vector<std::vector<VertexID>> newGroup = currentGroup;
                        newColoredV.push_back(group[0]);
                        newGroup.push_back(group);
                        coloredQ.push(newColoredV);
                        groupQ.push(newGroup);
                    }
                }
            }
        }
        else {
//            // Refine currentGroup to keep only minimal elements in each rule
//            if (!currentGroup.empty()) {
//                // Build partial order from the rules in currentGroup
//                std::map<VertexID, std::set<VertexID>> partialOrder;
//
//                // Initialize empty sets for all vertices that appear in rules
//                for (const auto &rule : currentGroup) {
//                    for (VertexID u : rule) {
//                        if (partialOrder.find(u) == partialOrder.end()) {
//                            partialOrder[u] = std::set<VertexID>();
//                        }
//                    }
//                }
//
//                // Build transitive closure of partial order: u < v means u is smaller than v
//                for (const auto &rule : currentGroup) {
//                    for (int j = 1; j < rule.size(); ++j) {
//                        partialOrder[rule[0]].insert(rule[j]);
//                    }
//                }
//
//                // Compute transitive closure
//                bool changed = true;
//                while (changed) {
//                    changed = false;
//                    for (auto &[u, larger] : partialOrder) {
//                        std::set<VertexID> newLarger = larger;
//                        for (VertexID v : larger) {
//                            for (VertexID w : partialOrder[v]) {
//                                if (newLarger.find(w) == newLarger.end()) {
//                                    newLarger.insert(w);
//                                    changed = true;
//                                }
//                            }
//                        }
//                        larger = newLarger;
//                    }
//                }
//
//                // Helper function to get minimal elements
//                auto getMinimalElements = [&](const std::vector<VertexID> &vertices) -> std::vector<VertexID> {
//                    std::vector<VertexID> minimal;
//                    for (VertexID u : vertices) {
//                        bool isMinimal = true;
//                        for (VertexID v : vertices) {
//                            if (u != v && partialOrder[v].count(u) > 0) {
//                                isMinimal = false;
//                                break;
//                            }
//                        }
//                        if (isMinimal) {
//                            minimal.push_back(u);
//                        }
//                    }
//                    return minimal;
//                };
//
//                // Refine each rule: keep rule[0] and minimal elements from rule[1] to rule[last]
//                for (auto &rule : currentGroup) {
//                    if (rule.size() > 1) {
//                        VertexID firstElement = rule[0];
//                        std::vector<VertexID> remainingElements(rule.begin() + 1, rule.end());
//                        std::vector<VertexID> minimalRemaining = getMinimalElements(remainingElements);
//
//                        rule.clear();
//                        rule.push_back(firstElement);
//                        rule.insert(rule.end(), minimalRemaining.begin(), minimalRemaining.end());
//                    }
//                }
//            }

            std::sort(currentGroup.begin(), currentGroup.end());
            bool exists = false;
            for (const auto &group: _candidateRules) {
                if (group == currentGroup) {
                    exists = true;
                    break;
                }
            }
            if (!exists) _candidateRules.push_back(currentGroup);
        }
    }
}

Edge * PatternGraph::coreUndirectedEdges(ui &num) const {
    Edge *edgeList = new Edge[_numEdges / 2];
    num = 0;
    for (ui i = 0; i < _coreSize; ++i) {
        for (ui j = i + 1; j < _coreSize; ++j) {
            VertexID u1 = _coreV[i], u2 = _coreV[j];
            if (isEdge(u1, u2)) {
                edgeList[num] = std::make_pair(u1, u2);
                ++num;
            }
        }
    }

    return edgeList;
}

void PatternGraph::printGraph(bool direction) const {
    printf("graph in CSR:\n");
    printf("vertices: %d, edges: %d\n", _numVertices, _numEdges);
    printf("from   to\n");
    for (VertexID i = 0; i < _numVertices; i++)
        for (EdgeID e = _offsets[i]; e < _offsets[i + 1]; e++) {
            VertexID j = _nbrs[e];
            if (direction)
                printf("%4d %4d\n", i, j);
            else
                printf("%4d %4d\n", j, i);
        }
}

bool PatternGraph::isClique() const {
    return _numVertices > 2 && _numVertices * (_numVertices - 1) == _numEdges;
}

void PatternGraph::buildCorePeripheral() {
    _coreV = new VertexID[_numVertices];
    _peripheralV = new VertexID[_numVertices];
    ui peripheralSize = 0;
    bool *visited = new bool[_numVertices + 1];
    int *pos = new int[_numVertices + 1];
    int *bin = new int[_numVertices + 1];
    int *deg = new int[_numVertices + 1];
    int *vert = new int[_numVertices + 1];
    memset(visited, false, sizeof(bool) * (_numVertices + 1));
    memset(pos, 0, sizeof(int) * (_numVertices + 1));
    memset(bin, 0, sizeof(int) * (_numVertices + 1));
    memset(deg, 0, sizeof(int) * (_numVertices + 1));
    memset(vert, 0, sizeof(int) * (_numVertices + 1));
    int maxDeg = 0;
    for (int v = 1; v <= _numVertices; ++v) {
        int d = (int) getDegree(v - 1);
        deg[v] = d;
        if (d > maxDeg) maxDeg = d;
    }
    for (int v = 1; v <= _numVertices; ++v)
        ++bin[deg[v]];
    int start = 1;
    for (int d = 0; d <= maxDeg; ++d) {
        int num = bin[d];
        bin[d] = start;
        start += num;
    }
    for (int v = 1; v <= _numVertices; ++v) {
        pos[v] = bin[deg[v]];
        vert[pos[v]] = v;
        ++bin[deg[v]];
    }
    for (int d = maxDeg; d > 0; --d)
        bin[d] = bin[d-1];
    bin[0] = 1;
    for (int i = 1; i <= _numVertices; ++i) {
        VertexID v = vert[i];
        if (deg[v] < 2){
            // v is a peripheral vertex
            _peripheralV[peripheralSize] = v - 1;
            visited[v - 1] = true;
            ++peripheralSize;
            ui nbrCnt = 0;
            VertexID *neighbors = getNeighbors(v - 1, nbrCnt);
            for (ui j = 0; j < nbrCnt; ++j) {
                VertexID u = neighbors[j] + 1;
                if (deg[u] > deg[v]) {
                    int du = deg[u], pu = pos[u], pw = bin[du], w = vert[pw];
                    if(u != w) {
                        pos[u] = pw;
                        vert[pu] = w;
                        pos[w] = pu;
                        vert[pw] = int(u);
                    }
                    ++bin[du];
                    --deg[u];
                }
            }
        }
        else break;
    }

    _coreSize = _numVertices - peripheralSize;
    ui coreIdx = 0;
    for (VertexID i = 0; i < _numVertices; ++i) {
        if (!visited[i]) {
            _coreV[coreIdx] = i;
            ++coreIdx;
        }
    }
    delete[] visited;
    delete[] pos;
    delete[] bin;
    delete[] deg;
    delete[] vert;
}

CanonType subgraphCanonValue(const PatternGraph &p, const std::vector<VertexID> &v, int *v2o) {
    CanonType canonValue = 0;
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(graph,cg,cg_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    options.digraph = TRUE;
    statsblk stats;
    int n = (int)v.size();
    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC2(graph,cg,cg_sz,m,n,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    EMPTYGRAPH(g,m,n);
    for (int i = 0; i < v.size(); ++i) {
        for (int j = 0; j < v.size(); ++j) {
            if (p.isEdge(v[i], v[j]))
                        ADDONEARC(g, i, j, m);
        }
    }
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
    // save the canonical value of the pattern
    int pos = 0;
    for (int i = 0; i < m * n; ++i) {
        unsigned long row = cg[i];
        bool current;
        for (unsigned long j = 0; j < n; ++j) {
            if (i != j) {
                if (1 & row >> (sizeof(unsigned long) * 8 - j - 1)) current = true;
                else current = false;
                if (current)
                    canonValue += 1 << pos;
                ++pos;
            }
        }
    }
    if (v2o != nullptr) {
        std::map<VertexID, int> canonV2O;
        int *v2l = new int[n];
        for (int i = 0; i < n; ++i) {
            v2l[lab[i]] = i;
        }
        // group canonical vertices according to orbits and reorder groups according to the min canonical vertex id
        int o = 0;
        std::vector<std::vector<int>> canonGroup(n);
        std::vector<int> canonV2Group(n);
        std::vector<int> group2O(n, -1);
        for (int i = 0; i < n; ++i) {
            canonGroup[orbits[i]].push_back(v2l[i]);
            canonV2Group[v2l[i]] = orbits[i];
        }
        for (int i = 0; i < n; ++i) {
            int group = canonV2Group[i];
            if (group2O[group] == -1) {
                group2O[group] = o;
                ++o;
            }
        }
        for (int i = 0; i < n; ++i) {
            for (int u: canonGroup[i]) {
                canonV2O[u] = group2O[i];
            }
        }
        for (int i = 0; i < n; ++i) {
            v2o[v[i]] = canonV2O[v2l[i]];
        }
        delete[] v2l;
    }
    return canonValue;
}
