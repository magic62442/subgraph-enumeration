//
// Created by anonymous authors on 2024/3/2.
//

#include "utils.h"

double getMemoryUsageGB() {
    std::ifstream statm("/proc/self/statm");
    long pageSize = sysconf(_SC_PAGESIZE); // Get system page size in bytes

    std::string line;
    double totalProgramSizeGB = 0.0;
    double residentSetSizeGB = 0.0;

    if (std::getline(statm, line)) {
        std::istringstream iss(line);
        long totalProgramPages, residentSetPages;
        if (iss >> totalProgramPages >> residentSetPages) {
            totalProgramSizeGB = totalProgramPages * pageSize / 1024.0 / 1024.0 / 1024.0;
            residentSetSizeGB = residentSetPages * pageSize / 1024.0 / 1024.0 / 1024.0;
        }
    }
    statm.close();

//    std::cout << "Total Program Size: " << totalProgramSizeGB << " GB" << std::endl;
//    std::cout << "Resident Set Size: " << residentSetSizeGB << " GB" << std::endl;
    return residentSetSizeGB;
}

bool nlcValid(const std::map<LabelID, ui> &queryNLC, const std::map<LabelID, ui> &dataNLC) {
    if (queryNLC.size() > dataNLC.size()) return false;
    for (auto it : queryNLC) {
        VertexID uL = it.first;
        ui uCount = it.second;
        if (dataNLC.find(uL) == dataNLC.end() || dataNLC.at(uL) < uCount)
            return false;
    }

    return true;
}

void sampleKElements(VertexID *array, ui length, ui k) {
    std::uniform_int_distribution<ui> dis(0, length - 1);
    std::set<ui> sampledIndices;
    while (sampledIndices.size() < k) {
        ui r = dis(gen);
        sampledIndices.insert(r);
    }
    ui i = 0;
    for (ui index: sampledIndices) {
        array[i] = array[index];
        ++i;
    }
}

void sampleKElements(VertexID *array, ui length, VertexID *samples, ui k) {
    std::uniform_int_distribution<ui> dis(0, length - 1);
    std::set<ui> sampledIndices;
    while (sampledIndices.size() < k) {
        ui r = dis(gen);
        sampledIndices.insert(r);
    }
    ui i = 0;
    for (ui index: sampledIndices) {
        samples[i] = array[index];
        ++i;
    }
}

int quickPow10(int n)
{
    static int pow10[10] = {
            1, 10, 100, 1000, 10000,
            100000, 1000000, 10000000, 100000000, 1000000000
    };

    return pow10[n];
}

Count choosec(Count n, int k)
{
    if (k > n)
        return 0;
    if (k == 0)
        return 1;
    Count nom = n;
    Count denom = 1;
    for (int i = 1; i < k; i++)
    {
        nom *= (n-i);
        denom *= (i+1);
    }
    return nom / denom;
}

void leftShit(std::vector<bool> &arr, int k) {
    int count = 0;
    for (int j = k; j >= 0; --j)
        if ( arr[j] == 1)  ++count;
    if (count == 0)  return;
    for (int j = k; j >= 0; --j)
        arr[j] = false;
    for ( int j=0; j<count;++j)
        arr[j] = true;
}

std::vector<std::vector<bool>> chooseK(ui n, int k) {
    if (k > n) {
        std::vector<std::vector<bool>> result;
        result.emplace_back(n, true);
        return result;
    }
    std::vector<bool> arr(n, false);
    for (int i = 0; i < k; ++i)
        arr[i] = true;
    std::vector<std::vector<bool>> result;
    result.push_back(arr);
    bool flag = true;
    while (flag) {
        flag = false;
        for (int i = 0; i < n - 1; ++i) {
            if (arr[i] && !arr[i + 1]){
                flag = true;
                arr[i] = false;
                arr[i + 1] = true;
                leftShit(arr,i);
                result.push_back(arr);
                break;
            }
        }
    }

    return result;
}

void sortEOrbitAggrePos(std::vector<int> &aggrePos) {
    std::vector<std::pair<int, int>> tmp(aggrePos.size() / 2);
    for (int i = 0; i < aggrePos.size(); i = i + 2) {
        if (aggrePos[i] < aggrePos[i + 1])
            tmp[i / 2] = std::make_pair(aggrePos[i], aggrePos[i + 1]);
        else
            tmp[i / 2] = std::make_pair(aggrePos[i + 1], aggrePos[i]);
    }
    std::sort(tmp.begin(), tmp.end());
    for (int i = 0; i < aggrePos.size(); i = i + 2) {
        aggrePos[i] = tmp[i / 2].first;
        aggrePos[i + 1] = tmp[i / 2].second;
    }
}

ui firstPosGreaterThan(const ui *candidate, ui begin, ui end, ui target) {
    int low = (int)begin;
    int high = (int)end - 1;
    int mid;
    while (low <= high) {
        mid = low + ((high - low) >> 1);
#ifndef __APPLE__
        _mm_prefetch((char *) &candidate[(mid + 1 + high) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &candidate[(mid - 1 + low) / 2], _MM_HINT_T0);
#endif
        if (candidate[mid] == target) {
            return mid + 1;
        } else if (candidate[mid] < target) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }

    return (ui)low;
}

// Backtracking function to generate all permutations
void generatePermutations(ui n, ui** automorphisms, ui* current, bool* used, size_t& index, size_t currentSize) {
    if (currentSize == n) {
        for (ui j = 0; j < n; j++) {
            automorphisms[index][j] = current[j];
        }
        index++;
        return;
    }

    for (ui i = 0; i < n; i++) {
        if (!used[i]) {
            used[i] = true;
            current[currentSize] = i;
            generatePermutations(n, automorphisms, current, used, index, currentSize + 1);
            used[i] = false;
        }
    }
}

void generateCliqueRulesHelper(const std::vector<VertexID> &vertices,
                               std::vector<std::vector<std::vector<VertexID>>> &rules,
                               std::vector<std::vector<VertexID>> current) {
    if (vertices.size() == 2) {
        current.push_back(vertices);
        rules.push_back(current);
        return;
    }
    for (int i = 0; i < vertices.size(); ++i) {
        std::vector<VertexID> newVertices;
        std::vector<VertexID> verticesRule;
        newVertices.reserve(vertices.size());
        verticesRule.reserve(vertices.size());
        verticesRule.push_back(vertices[i]);
        for (int j = 0; j < vertices.size(); ++j) {
            if (j != i) {
                newVertices.push_back(vertices[j]);
                verticesRule.push_back(vertices[j]);
            }
        }
        current.push_back(verticesRule);
        generateCliqueRulesHelper(newVertices, rules, current);
        current.pop_back();
    }
}

std::vector<std::vector<std::vector<VertexID>>> generateCliqueRules(const std::vector<VertexID> &vertices,
                                                                    const std::vector<std::vector<std::vector<VertexID>>> &candidateRules) {
    if (vertices.size() < 2) return candidateRules;
    std::vector<std::vector<VertexID>> current;
    std::vector<std::vector<std::vector<VertexID>>> rules;
    std::vector<std::vector<std::vector<VertexID>>> newCandRules;
    generateCliqueRulesHelper(vertices, rules, current);
    if (candidateRules.empty()) newCandRules = rules;
    else {
        for (int i = 0; i < candidateRules.size(); ++i) {
            for (int j = 0; j < rules.size(); ++j) {
                std::vector<std::vector<VertexID>> newRule = candidateRules[i];
                for (const auto &rule : rules[j])
                    newRule.push_back(rule);
                newCandRules.push_back(newRule);
            }
        }
    }

    return newCandRules;
}

int getPosition(uint64_t id, ui k, const std::vector<uint64_t> &allSets) {
    auto it = std::lower_bound(allSets.begin(), allSets.end(), id);
    if (it != allSets.end() && *it == id) {
        return std::distance(allSets.begin(), it);
    }
    return -1; // Not found
}

void subsetSupersetRelationships(ui n, const std::vector<VertexID> &elements,
                                 const std::vector<std::vector<uint64_t>> &subsets, std::vector<std::vector<std::vector<int>>> &subsetOf,
                                 std::vector<std::vector<std::vector<int>>> &supersetOf) {
    subsetOf.resize(n + 1);
    supersetOf.resize(n + 1);
    for (int k = 1; k <= n; ++k) {
        subsetOf[k].resize(subsets[k].size());
        if (k > 1) {
            for (size_t i = 0; i < subsets[k].size(); ++i) {
                uint64_t subsetID = subsets[k][i];
                for (const auto& element : elements) {
                    uint64_t elementBit = 1ULL << (element % 64);
                    if (subsetID & elementBit) {
                        uint64_t smallerSubsetID = subsetID & ~elementBit;
                        int pos = getPosition(smallerSubsetID, k - 1, subsets[k - 1]);
                        if (pos != -1) {
                            subsetOf[k][i].push_back(pos);
                        }
                    }
                }
            }
        }
        supersetOf[k].resize(subsets[k].size());
        if (k < n) {
            for (size_t i = 0; i < subsets[k].size(); ++i) {
                uint64_t subsetID = subsets[k][i];
                for (const auto& element : elements) {
                    uint64_t elementBit = 1ULL << (element % 64);
                    if (!(subsetID & elementBit)) {
                        uint64_t superSetID = subsetID | elementBit;
                        int pos = getPosition(superSetID, k + 1, subsets[k + 1]);
                        if (pos != -1) {
                            supersetOf[k][i].push_back(pos);
                        }
                    }
                }
            }
        }
    }
}

size_t numUniqueCombine(std::vector<std::vector<VertexID>> &vsets) {
    size_t num = 1;
    for (auto &vset : vsets) num *= vset.size();
    std::map<uint64_t, std::vector<VertexID>> id2Intersection;
    std::vector<std::pair<int, int>> pairs;
    // pairwise intersection
    for (int i = 0; i < vsets.size(); ++i) {
        uint64_t id = 1 << i;
        id2Intersection[id] = vsets[i];
        for (int j = i + 1; j < vsets.size(); ++j) {
            id  = (1 << i) + (1 << j);
            id2Intersection[id] = std::vector<VertexID>();
            std::vector<VertexID> &intersect = id2Intersection[id];
            if (vsets[i].size() < vsets[j].size()) intersect.reserve(vsets[i].size());
            else intersect.reserve(vsets[j].size());
            std::set_intersection(vsets[i].begin(), vsets[i].end(), vsets[j].begin(),
                                  vsets[j].end(), std::back_inserter(intersect));
            pairs.emplace_back(i, j);
        }
    }
    int n = vsets.size();
    for (int k = 1; k <= pairs.size(); ++k) {
        // choose k sets in n(n-1)/2 pairwises
        std::vector<std::vector<bool>> choices = chooseK(pairs.size(), k);
        for (auto &choice : choices) {
            UnionFind uf(n);
            for (int i = 0; i < pairs.size(); ++i) {
                if (choice[i])
                    uf.unite(pairs[i].first, pairs[i].second);
            }
            std::vector<std::vector<int>> disjointSets = uf.getDisjointSets();
            size_t term = 1;
            for (auto &disjointSet : disjointSets) {
                uint64_t id = 0;
                for (int pos : disjointSet) id |= (1 << pos);
                if (id2Intersection.find(id) != id2Intersection.end())
                    term *= id2Intersection[id].size();
                else {
                    for (int pos : disjointSet) {
                        uint64_t subsetID = id - (1 << pos);
                        if (id2Intersection.find(subsetID) != id2Intersection.end()) {
                            const std::vector<VertexID> &oldIntersect = id2Intersection[subsetID];
                            id2Intersection[id] = std::vector<VertexID>();
                            std::vector<VertexID> &intersect = id2Intersection[id];
                            if (oldIntersect.size() < vsets[pos].size()) intersect.reserve(oldIntersect.size());
                            else intersect.reserve(vsets[pos].size());
                            std::set_intersection(oldIntersect.begin(), oldIntersect.end(), vsets[pos].begin(),
                                                  vsets[pos].end(), std::back_inserter(intersect));
                            term *= intersect.size();
                            break;
                        }
                    }
                }
            }
            if (k % 2 == 0) num += term;
            else num -= term;
        }
    }

    return num;
}

bool unorderedSubset(const std::vector<VertexID> &v1, const std::vector<VertexID> &v2) {
    std::unordered_set<VertexID> set2(v2.begin(), v2.end());

    int intersectionCount = 0;
    for (VertexID id : v1) {
        if (set2.find(id) != set2.end()) {
            intersectionCount++;
        }
    }

    if (intersectionCount == v1.size()) return true;
    return false;
}

bool orderedSubset(const std::vector<VertexID> &v1, const std::vector<VertexID> &v2) {
    if (v1.size() > v2.size()) return false;
    for (int i = 0; i < v1.size(); ++i) {
        if (v1[i] != v2[i]) return false;
    }

    return true;
}

size_t bfUniqueCombine(std::vector<std::vector<VertexID>> &vsets) {
    size_t num = 0;
    ui n = vsets.size();
    std::vector<int> poses(n, 0);
    int depth = 0;
    std::vector<VertexID> tuple(n);
    while (depth >= 0) {
        while (poses[depth] < vsets[depth].size()) {
            VertexID v = vsets[depth][poses[depth]];
            tuple[depth] = v;
            ++poses[depth];
            bool exists = false;
            for (int i = 0; i < depth; ++i) {
                if (tuple[i] == v) {
                    exists = true;
                    break;
                }
            }
            if (exists) continue;
            if (depth == n - 1)
                ++num;
            else {
                ++depth;
                poses[depth] = 0;
            }
        }
        --depth;
    }

    return num;
}

int findFirstExtension(const std::vector<std::vector<uint32_t>>& vec, const std::vector<uint32_t>& prefix) {
    int left = 0;
    int right = vec.size();

    while (left < right) {
        int mid = left + (right - left) / 2;
        bool vecLarger = true;
        for (int i = 0; i < prefix.size(); ++i) {
            if (vec[mid][i] < prefix[i]) {
                vecLarger = false;
                break;
            }
            else if(vec[mid][i] > prefix[i]) break;
        }
        if (vecLarger) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }

    return left;
}

int findFirstGreater(const std::vector<std::vector<uint32_t>>& vec, const std::vector<uint32_t>& prefix) {
    int left = 0;
    int right = vec.size();

    while (left < right) {
        int mid = left + (right - left) / 2;
        bool vecSmaller = true;
        for (int i = 0; i < prefix.size(); ++i) {
            if (vec[mid][i] > prefix[i]) {
                vecSmaller = false;
                break;
            }
            else if (vec[mid][i] < prefix[i]) break;
        }
        if (vecSmaller) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    return left;
}

void makeContinuous(std::vector<VertexID>& nums) {
    // Copy and sort the original vector
    std::vector<VertexID> sorted_nums = nums;
    std::sort(sorted_nums.begin(), sorted_nums.end());

    // Map each unique value to its new continuous value
    std::map<int, int> value_map;
    int new_value = 0;
    for (int num : sorted_nums) {
        if (value_map.find(num) == value_map.end()) {
            value_map[num] = new_value++;
        }
    }

    // Replace each value in the original vector with its new value
    for (auto& num : nums) {
        num = value_map[num];
    }
}

bool checkBudget(size_t oldBudget, size_t budget, const std::vector<size_t> &tupleSizes, const std::vector<std::vector<std::vector<VertexID>>> &tuples,
                 const std::vector<std::vector<std::vector<VertexID>>> &backup) {
    for (size_t sz: tupleSizes) budget += sz;
    for (int i = 0; i < tuples.size(); ++i) {
        if (!tuples[i].empty()) budget += tuples[i].size() * tuples[i][0].size() * sizeof(VertexID);
        if (!backup[i].empty()) budget += backup[i].size() * backup[i][0].size() * sizeof(VertexID);
    }

    return oldBudget == budget;
}

std::string pathJoin(const std::string& part1, const std::string& part2) {
    if (part1.empty()) {
        return part2;
    }
    if (part2.empty()) {
        return part1;
    }

    // Define the directory separator
#ifdef _WIN32
    const char sep = '\\';
#else
    const char sep = '/';
#endif

    std::string result = part1;

    // Add separator if needed
    if (result.back() != sep) {
        result += sep;
    }

    // Append the second part, avoiding double separator if it starts with one
    if (part2.front() == sep) {
        result += part2.substr(1);
    } else {
        result += part2;
    }

    return result;
}

std::vector<int> buildMaxSubsetChain(const std::vector<std::pair<int, int>>& subsetRelation, int totalElements) {
    std::vector<std::vector<int>> graph(totalElements);
    std::vector<int> inDegree(totalElements, 0);

    for (const auto& relation : subsetRelation) {
        int i = relation.first;
        int j = relation.second;
        graph[i].push_back(j);
        inDegree[j]++;
    }
    std::queue<int> queue;
    std::vector<int> distance(totalElements, -1);
    std::vector<int> parent(totalElements, -1);

    for (int i = 0; i < totalElements; i++) {
        if (inDegree[i] == 0) {
            queue.push(i);
            distance[i] = 1;
        }
    }

    int maxDistance = 0;
    int maxNode = -1;

    while (!queue.empty()) {
        int current = queue.front();
        queue.pop();

        if (distance[current] > maxDistance) {
            maxDistance = distance[current];
            maxNode = current;
        }

        for (int neighbor : graph[current]) {
            if (distance[current] + 1 > distance[neighbor]) {
                distance[neighbor] = distance[current] + 1;
                parent[neighbor] = current;
            }

            inDegree[neighbor]--;
            if (inDegree[neighbor] == 0) {
                queue.push(neighbor);
            }
        }
    }

    std::vector<int> result;
    if (maxNode != -1) {
        int current = maxNode;
        std::vector<int> path;
        while (current != -1) {
            path.push_back(current);
            current = parent[current];
        }

        for (int i = path.size() - 1; i >= 0; i--) {
            result.push_back(path[i]);
        }
    }

    return result;
}

ui setBeginPos(const VertexID* array, ui count, VertexID maxCompared) {
    ui beginPos = 0;
    if (count <= BINARY_SEARCH_THRESHOLD) {
        for (ui j = 0; j < count; ++j) {
            if (array[j] > maxCompared) {
                beginPos = j;
                break;
            }
        }
        if (beginPos == 0 && count > 0 && array[count - 1] <= maxCompared) {
            beginPos = count;
        }
    } else {
        beginPos = std::upper_bound(array, array + count, maxCompared) - array;
    }
    return beginPos;
}

ui setEndPos(const VertexID* array, ui count, VertexID minCompared) {
    ui endPos = count;
    if (count <= BINARY_SEARCH_THRESHOLD) {
        for (ui j = 0; j < count; ++j) {
            if (array[j] >= minCompared) {
                endPos = j;
                break;
            }
        }
    } else {
        endPos = std::lower_bound(array, array + count, minCompared) - array;
    }
    return endPos;
}
