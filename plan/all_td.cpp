//
// Created by 李祺彦 on 25-6-8.
//

#include "all_td.h"
#include <unordered_set>
#include "subset_structure.h"

void findCliquesRecursive(const PatternGraph &graph,
                          std::vector<VertexID> &currentClique,
                          std::vector<VertexID> &potentialClique,
                          std::vector<VertexID> &processedVertices,
                          std::queue<HyperNode> &cliques) {
    if (potentialClique.empty() && processedVertices.empty()) {
        HyperNode node;
        node.numAttributes = currentClique.size();
        node.attributes = new VertexID[node.numAttributes];
        for (int i = 0; i < currentClique.size(); ++i) {
            node.attributes[i] = currentClique[i];
            node.id += 1 << currentClique[i];
        }
        std::copy(currentClique.begin(), currentClique.end(), node.attributes);
        cliques.push(node);
        return;
    }
    for (int i = 0; i < potentialClique.size(); ++i) {
        std::vector<VertexID> newCurrentClique = currentClique;
        VertexID newVertex = potentialClique[i];
        newCurrentClique.push_back(newVertex);
        std::vector<VertexID> newPotentialClique, newProcessedClique;
        // refine potentialClique and processedVertices
        for (auto & u: potentialClique) {
            if (u > newVertex && graph.isEdge(u, newVertex))
                newPotentialClique.push_back(u);
        }
        for (auto & u: processedVertices) {
            if (graph.isEdge(u, newVertex))
                newProcessedClique.push_back(u);
        }
        findCliquesRecursive(graph, newCurrentClique, newPotentialClique, newProcessedClique, cliques);
        processedVertices.push_back(newVertex);
    }
}

std::queue<HyperNode> findMaximalCliques(const PatternGraph &graph) {
    std::queue<HyperNode> cliques;
    std::vector<VertexID> currentClique, potentialClique(graph.getNumVertices()), processedVertices;

    // Initialize potentialClique with all vertices
    for (ui i = 0; i < graph.getNumVertices(); ++i) {
        potentialClique[i] = i;
    }

    findCliquesRecursive(graph, currentClique, potentialClique, processedVertices, cliques);
    return cliques;
}

std::vector<HyperNode> getAllNodes(const PatternGraph &p) {
    std::vector<HyperNode> nodes;
    // fill the initial q with cliques
    std::queue<HyperNode> q = findMaximalCliques(p);
    ui coreSize;
    ui maxNodeID = 1 << p.getNumVertices();
    bool *visited = new bool[maxNodeID];
    memset(visited, false, sizeof(bool) * maxNodeID);
    VertexID *coreV = p.getCoreV(coreSize);
    // get core graph edges and push them into q.
    ui numEdge = 0;
    Edge* edgeList = p.coreUndirectedEdges(numEdge);
    while (!q.empty()) {
        HyperNode tau = q.front();
        q.pop();
        ui tauSize = tau.numAttributes;
        int id = tau.id;
        VertexID *vertices = tau.attributes;
        if (tauSize > 2) {
            visited[id] = true;
            nodes.push_back(tau);
        }
        for (int i = 0; i < coreSize; ++i) {
            VertexID u1 = coreV[i];
            if (tau.hasVertex(u1)) continue;
            int newID = id + (1 << u1);
            if (visited[newID]) continue;
            // loop over vertices in tau, check whether u1 connects to one of them
            for (int j = 0; j < tauSize; ++j) {
                VertexID u2 = vertices[j];
                if (p.isEdge(u1, u2)) {
                    visited[newID] = true;
                    // vertices + u2 is weakly connected, thus a new node
                    // use insertion sort to keep newV sorted in ascending order
                    VertexID *newV = new VertexID[tauSize + 1];
                    int k;
                    for (k = 0; k < tauSize; ++k) {
                        if (vertices[k] < u1)
                            newV[k] = vertices[k];
                        else
                            break;
                    }
                    newV[k] = u1;
                    for (; k < tauSize; ++k)
                        newV[k + 1] = vertices[k];
                    q.emplace(newID, newV, tauSize + 1);
                    break;
                }
            }
        }
    }

    delete[] edgeList;
    delete[] visited;
    return nodes;
}

std::vector<HyperTree> getAllTree(const PatternGraph &p) {
    std::vector<HyperNode> allNodes = getAllNodes(p);
    return getAllTree(allNodes, p);
}

std::vector<HyperTree> getMinWidthTrees(const PatternGraph &p) {
    std::vector<HyperNode> allNodes = getAllNodes(p);
    std::vector<HyperTree> minWidthTrees;
    double fhtw = p.getNumEdges();
    ui tw = p.getNumVertices();
    std::vector<HyperTree> allTrees = getAllTree(allNodes, p);
    std::vector<double> fws(allTrees.size());
    std::vector<ui> widths(allTrees.size());
    for (int i = 0; i < allTrees.size(); ++i) {
        const HyperTree &t = allTrees[i];
        ui width = 0;
        double fw = 0.0;
        for (VertexID nID = 0; nID < t.numNodes; ++nID) {
            if (t.nodes[nID].numAttributes > width)
                width = t.nodes[nID].numAttributes;
            if (t.nodes[nID].fw > fw)
                fw = t.nodes[nID].fw;
        }
        fws[i] = fw;
        widths[i] = width;
        if (width < tw) tw = width;
        if (fw < fhtw) fhtw = fw;
    }
    for (int i = 0; i < allTrees.size(); ++i) {
        if (fws[i] == fhtw && widths[i] == tw) {
            allTrees[i].addPeripheral(p);
            minWidthTrees.push_back(allTrees[i]);
        }
    }
    return minWidthTrees;
}

std::vector<HyperTree> getAllTree(const std::vector<HyperNode> &allNode, const PatternGraph &p) {
    std::vector<HyperTree> trees;
    std::queue<HyperTree> q;
    std::vector<std::set<CanonType>> treeID;
    q.emplace(p.getNumVertices());
    while (!q.empty()) {
        HyperTree t = q.front();
        q.pop();
        if (isValid(t, p)) {
            std::set<CanonType> nodeIDs;
            for (VertexID nID = 0; nID < t.numNodes; ++nID) {
                HyperNode &bag = t.nodes[nID];
                std::vector<VertexID> vertices(bag.attributes, bag.attributes + bag.numAttributes);
                bag.canonValue = subgraphCanonValue(p, vertices, nullptr);
                if (canonToFW.find(bag.canonValue) == canonToFW.end()) {
                    fractionalWidth(bag, p);
                    canonToFW[bag.canonValue] = bag.fw;
                }
                else bag.fw = canonToFW[bag.canonValue];
                nodeIDs.insert(bag.canonValue);
            }
            bool idExists = false;
            for (int i = 0; i < treeID.size(); ++i) {
                if (treeID[i] == nodeIDs) {
                    idExists = true;
                    break;
                }
            }
            if (!idExists) {
                trees.push_back(t);
                treeID.push_back(nodeIDs);
                // early termination
                if (p.getNumVertices() > 8 && trees.size() >= 500) return trees;
            }
        }
        else {
            std::vector<HyperNode> candidateNodes = getCandidateNodes(t, allNode, p);
            for (auto &tau: candidateNodes) {
                std::vector<HyperTree> newTrees = addNode(tau, t, p);
                for (const auto& newT: newTrees) {
                    q.push(newT);
                }
            }
        }
    }

    return trees;
}

std::vector<HyperTree> addNode(HyperNode &tau, const HyperTree &tree, const PatternGraph &p) {
    std::vector<HyperTree> result;
    if (tree.numNodes == 0) {
        HyperTree t(tree);
        delete[] t.nodes;
        t.nodes = new HyperNode[t.numAttributes];
        t.v2n = new std::vector<VertexID>[t.numAttributes];
        t.nodes[0] = tau;
        for (int i = 0; i < tau.numAttributes; ++i) {
            VertexID u = tau.attributes[i];
            t.v2n[u].push_back(0);
        }
        t.numNodes = 1;
        result.push_back(t);
        return result;
    }
    // for each existing node that shares some vertices with the new node
    // we add an edge and check whether it is a valid TD
    for (ui i = 0; i < tree.numNodes; ++i) {
        HyperTree t(tree);
        HyperNode *old = t.nodes;
        t.nodes = new HyperNode[t.numAttributes];
        for (int j = 0; j < t.numNodes; ++j) {
            old[j].copyTo(t.nodes[j]);
        }
        delete[] old;
        t.nodes[t.numNodes] = tau;
        // for each vertex u in tau, add tau to _v2n[u]
        for (int j = 0; j < tau.numAttributes; ++j) {
            VertexID u = tau.attributes[j];
            t.v2n[u].push_back(t.numNodes);
        }
        ++t.numNodes;
        bool flag;
        for (int j = 0; j < tau.numAttributes; ++j) {
            VertexID u = tau.attributes[j];
            if (t.nodes[i].hasVertex(u)) {
                // add (i, _numNodes - 1) to edges
                flag = true;
                t.edges[i].push_back(t.numNodes - 1);
                t.edges[t.numNodes - 1].push_back(i);
                // check whether t satisfies the TD constraint
                // i.e. for any vertex v, tree nodes containing v are connected
                int k;
                for (k = 0; k < tau.numAttributes; ++k) {
                    VertexID _u = tau.attributes[k];
                    std::vector<ui> pos = t.v2n[_u];
                    if (pos.size() == 1) continue;
                    // check whether the subgraph of t._edges induced by pos is connected
                    // start a bfs at pos[0]
                    std::vector<bool> visited(t.numAttributes, true);
                    for (auto nd : pos)
                        visited[nd] = false;

                    std::queue<VertexID> q;
                    q.push(pos[0]);
                    visited[pos[0]] = true;
                    while (!q.empty()) {
                        VertexID nd1 = q.front();
                        q.pop();
                        for (auto nd2: t.edges[nd1]) {
                            if (!visited[nd2]) {
                                visited[nd2] = true;
                                q.push(nd2);
                            }
                        }
                    }
                    // if there are unvisited nodes, then nodes containing _u is unconnected, and t is invalid
                    for (auto f: visited) {
                        if (!f) flag = false;
                    }
                }
                if (flag)
                    result.push_back(t);
                break;
            }
        }
    }

    return result;
}

std::vector<HyperNode> getCandidateNodes(const HyperTree &t, const std::vector<HyperNode> &allNodes, const PatternGraph &p) {
    std::vector<HyperNode> result;
    std::vector<VertexID> unCovered = getUncovered(t, p);
    for (const auto& tau: allNodes) {
        // if tau is neither a supernode nor a subnode of some nodes in t,
        // then tau can be added to tree t
        if (t.hasSubNodeOf(tau) || t.hasSupNodeOf(tau)) continue;
        for (auto u: unCovered)
            if (tau.hasVertex(u)) {
                result.push_back(tau);
                break;
            }
    }
    return result;
}

void fractionalWidth(HyperNode &tau, const PatternGraph &p) {
    glp_term_out(GLP_OFF);

    ui n = 0, m = 0;
    Edge e[100];
    n = tau.numAttributes;
    for(ui i = 0; i < n; ++i)
        for(ui j = i + 1; j < n; ++j)
            if(p.isEdge(tau.attributes[i], tau.attributes[j]))
                e[m++] = std::make_pair(i, j);

    int ia[n * (n + m) + 1], ja[n * (n + m) + 1];
    double ar[n * (n + m) + 1], z;

    glp_prob *lp;

    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MIN);
    glp_add_rows(lp, n);

    for(size_t i = 0; i < n; ++i)
        glp_set_row_bnds(lp, i + 1, GLP_LO, 1., 1e8);

    glp_add_cols(lp, m + n);
    for(size_t i = 0; i < n + m; ++i) {
        glp_set_col_bnds(lp, i + 1, GLP_LO, 0., 0.);
        glp_set_obj_coef(lp, i + 1, 1.);
    }

    int cnt = 0;

    for(size_t i = 0; i < n; ++i)
        for(size_t j = 0; j < n; ++j)
            ia[++cnt] = i + 1, ja[cnt] = j + 1, ar[cnt] = (i == j ? 1. : 0);

    for(size_t j = 0; j < m; ++j) {
        size_t u = e[j].first, v = e[j].second;
        for(size_t i = 0; i < n; ++i)
            ia[++cnt] = i + 1, ja[cnt] = n + j + 1, ar[cnt] = ((i == u || i == v) ? 1. : 0);
    }

    glp_load_matrix(lp, cnt, ia, ja, ar);
    glp_simplex(lp, NULL);

    z = glp_get_obj_val(lp);
    glp_delete_prob(lp);

    tau.fw = z;
}

std::vector<VertexID> getUncovered(const HyperTree &t, const PatternGraph &p) {
    std::vector<VertexID> result;
    ui coreSize;
    VertexID *coreV = p.getCoreV(coreSize);
    for (ui i = 0; i < coreSize; ++i) {
        VertexID u = coreV[i];
        if (t.v2n[u].empty())
            result.push_back(u);
    }

    return result;
}

bool isValid(const HyperTree &t, const PatternGraph &p) {

    if (!getUncovered(t, p).empty()) return false;
    ui numEdge = 0;
    Edge* edgeList = p.coreUndirectedEdges(numEdge);
    std::set<ui> edgeInCore, edgesInTree;
    for (ui i = 0; i < numEdge; ++i) {
        VertexID u1 = edgeList[i].first, u2 = edgeList[i].second;
        // encode an edge into an integer
        edgeInCore.insert(100 * u1 + u2);
    }
    delete[] edgeList;
    for (ui i = 0; i < t.numNodes; ++i) {
        const HyperNode &tau = t.nodes[i];
        // loop over vertex pairs in tau
        for (int j = 0; j < tau.numAttributes; ++j) {
            VertexID u1 = tau.attributes[j];
            for (int k = j + 1; k < tau.numAttributes; ++k) {
                VertexID u2 = tau.attributes[k];
                if (p.isEdge(u1, u2)) {
                    edgesInTree.insert(100 * u1 + u2);
                }
            }
        }
    }
    if (edgeInCore == edgesInTree)
        return true;
    else
        return false;
}

bool nodeContains(const HyperNode &node1, const HyperNode &node2) {
    // Check if node1's attributes contain all of node2's attributes
    for (ui i = 0; i < node2.numAttributes; ++i) {
        VertexID v = node2.attributes[i];
        bool found = false;
        for (ui j = 0; j < node1.numAttributes; ++j) {
            if (node1.attributes[j] == v) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }
    return true;
}

bool findPerfectMatching(const HyperTree &t1, const HyperTree &t2,
                         std::vector<int> &matching, std::vector<bool> &used, int depth) {
    if (depth == t2.numNodes) {
        return true; // Found perfect matching
    }

    // Try to match t2.nodes[depth] with some unused node in t1
    for (VertexID i = 0; i < t1.numNodes; ++i) {
        if (used[i]) continue;

        if (nodeContains(t1.nodes[i], t2.nodes[depth])) {
            // Try this matching
            matching[depth] = i;
            used[i] = true;

            if (findPerfectMatching(t1, t2, matching, used, depth + 1)) {
                return true;
            }

            // Backtrack
            used[i] = false;
        }
    }

    return false;
}

bool containsHyperTree(const HyperTree &t1, const HyperTree &t2) {
    // Check if number of nodes is the same
    if (t1.numNodes != t2.numNodes) {
        return false;
    }

    // Try to find a perfect matching where each t1 node contains the corresponding t2 node
    std::vector<int> matching(t2.numNodes, -1);
    std::vector<bool> used(t1.numNodes, false);

    return findPerfectMatching(t1, t2, matching, used, 0);
}

double computeFractionalHypertreeWidth(const PatternGraph &p, const std::vector<size_t> &eliminationOrder, FHD &fhd) {
    // Convert PatternGraph to HyperG format
    HyperG H;
    H.N = p.getNumVertices();
    H.M = p.getNumEdges() / 2;
    H.e.clear();
    const EdgeID *offsets = p.getOffsets();
    const VertexID *nbrs = p.getNbors();
    for (VertexID u = 0; u < p.getNumVertices(); ++u) {
        for (EdgeID e = offsets[u]; e < offsets[u + 1]; ++e) {
            VertexID u2 = nbrs[e];
            if (u2 < u) continue;
            VertexSet temp;
            temp.Set(u);
            temp.Set(u2);
            H.e.push_back(temp);
        }
    }

    // Create FHD from elimination order
    fhd = FHD(H, eliminationOrder);
    fhd.Refine();

    // Calculate width using calc_width function
    std::vector<VertexSet> E = H.e;
    std::vector<size_t> order = eliminationOrder;
    double *w = new double[order.size()];

    double width = calc_width(H, order, E, w);

    delete[] w;
    return width;
}

std::vector<HyperTree> buildOptimalFHDsByEnumeration(const PatternGraph &p, std::vector<std::vector<size_t>> &optimalOrders, double &minWidth) {
    std::vector<HyperTree> optimalTrees;
    optimalOrders.clear();
    minWidth = 1e9; // Initialize with a very large value

    // Get canonical labeling from PatternGraph

    // Set to store visited canonical permutations
    std::unordered_set<Permutation> visitedPermutations;

    // Vector to store tree IDs for deduplication (similar to getAllTree)
    std::vector<std::set<CanonType>> treeID;

    // Generate all vertices
    std::vector<size_t> vertices;
    ui coreSize;
    VertexID *coreV = p.getCoreV(coreSize);
    for (size_t i = 0; i < coreSize; ++i) {
        vertices.push_back(coreV[i]);
    }

    // Enumerate all permutations
    do {
        // Convert current permutation to canonical labels
        std::vector<uint32_t> canonicalOrbits;
        for (size_t vertex : vertices) {
            canonicalOrbits.push_back(static_cast<uint32_t>(p.getOrbit(vertex)));
        }

        // Encode the canonical permutation
        Permutation encodedPermutation = encodePermutation(canonicalOrbits);

        // Check if this canonical permutation has been visited
        if (visitedPermutations.find(encodedPermutation) != visitedPermutations.end()) {
            continue; // Skip this permutation
        }

        // Mark this canonical permutation as visited
        visitedPermutations.insert(encodedPermutation);

        FHD currentFHD;
        double currentWidth = computeFractionalHypertreeWidth(p, vertices, currentFHD);

        // Keep track of all optimal solutions
        if (currentWidth < minWidth) {
            // Found a better solution, clear previous results
            minWidth = currentWidth;
            optimalTrees.clear();
            optimalOrders.clear();
            treeID.clear(); // Clear tree IDs for deduplication

            // Build HyperTree from current FHD
            HyperTree currentTree;
            currentTree.buildFromTD(currentFHD);

            // Check if each bag forms a connected subgraph
            bool allBagsConnected = true;
            for (VertexID nID = 0; nID < currentTree.numNodes; ++nID) {
                HyperNode &bag = currentTree.nodes[nID];
                std::vector<VertexID> nodeVertices(bag.attributes, bag.attributes + bag.numAttributes);
                if (!p.isConnected(nodeVertices)) {
                    allBagsConnected = false;
                    break;
                }
            }

            if (!allBagsConnected) {
                continue; // Skip this tree if any bag is not connected
            }

            // Compute canonical values and fractional widths for all nodes, and build nodeIDs
            std::set<CanonType> nodeIDs;
            for (VertexID nID = 0; nID < currentTree.numNodes; ++nID) {
                HyperNode &bag = currentTree.nodes[nID];
                std::vector<VertexID> nodeVertices(bag.attributes, bag.attributes + bag.numAttributes);
                bag.canonValue = subgraphCanonValue(p, nodeVertices, nullptr);
                if (canonToFW.find(bag.canonValue) == canonToFW.end()) {
                    fractionalWidth(bag, p);
                    canonToFW[bag.canonValue] = bag.fw;
                }
                else bag.fw = canonToFW[bag.canonValue];
                nodeIDs.insert(bag.canonValue);
            }

            optimalTrees.push_back(currentTree);
            optimalOrders.push_back(vertices);
            treeID.push_back(nodeIDs);
        } else if (currentWidth == minWidth) {
            // Found another optimal solution
            HyperTree currentTree;
            currentTree.buildFromTD(currentFHD);

            // Check if each bag forms a connected subgraph
            bool allBagsConnected = true;
            for (VertexID nID = 0; nID < currentTree.numNodes; ++nID) {
                HyperNode &bag = currentTree.nodes[nID];
                std::vector<VertexID> nodeVertices(bag.attributes, bag.attributes + bag.numAttributes);
                if (!p.isConnected(nodeVertices)) {
                    allBagsConnected = false;
                    break;
                }
            }

            if (!allBagsConnected) {
                continue; // Skip this tree if any bag is not connected
            }

            // Compute canonical values and fractional widths for all nodes, and build nodeIDs
            std::set<CanonType> nodeIDs;
            for (VertexID nID = 0; nID < currentTree.numNodes; ++nID) {
                HyperNode &bag = currentTree.nodes[nID];
                std::vector<VertexID> nodeVertices(bag.attributes, bag.attributes + bag.numAttributes);
                bag.canonValue = subgraphCanonValue(p, nodeVertices, nullptr);
                if (canonToFW.find(bag.canonValue) == canonToFW.end()) {
                    fractionalWidth(bag, p);
                    canonToFW[bag.canonValue] = bag.fw;
                }
                else bag.fw = canonToFW[bag.canonValue];
                nodeIDs.insert(bag.canonValue);
            }

            // Check if this tree already exists (deduplication)
            bool idExists = false;
            for (int i = 0; i < treeID.size(); ++i) {
                if (treeID[i] == nodeIDs) {
                    idExists = true;
                    break;
                }
            }

            // Only add if it's not a duplicate
            if (!idExists) {
                optimalTrees.push_back(currentTree);
                optimalOrders.push_back(vertices);
                treeID.push_back(nodeIDs);
            }
        }

    } while (std::next_permutation(vertices.begin(), vertices.end()));

    // 为每个optimal tree的每个bag计算id
    for (auto &tree : optimalTrees) {
        for (VertexID nID = 0; nID < tree.numNodes; ++nID) {
            HyperNode &bag = tree.nodes[nID];
            CanonType id = 0;
            for (ui i = 0; i < bag.numAttributes; ++i) {
                id |= (1ULL << bag.attributes[i]);
            }
            bag.id = id;
        }
        tree.addPeripheral(p);
    }

    // Filter by treeWidth (maximum numAttributes among all HyperNodes in each tree)
    if (!optimalTrees.empty()) {
        // Calculate treeWidth for each tree and find minimum
        int minTreeWidth = INT_MAX;

        for (const auto &tree : optimalTrees) {
            int maxNumAttributes = 0;
            for (VertexID nID = 0; nID < tree.numNodes; ++nID) {
                maxNumAttributes = std::max(maxNumAttributes, (int)tree.nodes[nID].numAttributes);
            }
            minTreeWidth = std::min(minTreeWidth, maxNumAttributes);
        }

        // Filter trees to keep only those with minimum treeWidth
        std::vector<HyperTree> filteredTrees;
        std::vector<std::vector<size_t>> filteredOrders;
        for (size_t i = 0; i < optimalTrees.size(); ++i) {
            int maxNumAttributes = 0;
            for (VertexID nID = 0; nID < optimalTrees[i].numNodes; ++nID) {
                maxNumAttributes = std::max(maxNumAttributes, (int)optimalTrees[i].nodes[nID].numAttributes);
            }
            if (maxNumAttributes == minTreeWidth) {
                filteredTrees.push_back(optimalTrees[i]);
                filteredOrders.push_back(optimalOrders[i]);
            }
        }

        optimalTrees = filteredTrees;
        optimalOrders = filteredOrders;
    }

    // Filter to keep only minimal trees (not contained by any other optimal tree)
    if (!optimalTrees.empty()) {
        std::vector<HyperTree> minimalTrees;
        std::vector<std::vector<size_t>> minimalOrders;

        for (size_t i = 0; i < optimalTrees.size(); ++i) {
            bool isMinimal = true;
            for (size_t j = 0; j < optimalTrees.size(); ++j) {
                if (i != j && containsHyperTree(optimalTrees[i], optimalTrees[j])) {
                    isMinimal = false;
                    break;
                }
            }
            if (isMinimal) {
                minimalTrees.push_back(optimalTrees[i]);
                minimalOrders.push_back(optimalOrders[i]);
            }
        }

        optimalTrees = minimalTrees;
        optimalOrders = minimalOrders;
    }

    return optimalTrees;
}
