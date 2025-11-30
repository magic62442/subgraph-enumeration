//
// Created by anonymous authors on 2024/2/27.
//

#include "relation.h"

std::vector<std::vector<VertexID>> gResult1 = {};

size_t TrieNode::memoryCost() const {
    size_t mem = sizeof(*this);
    for (int i = 0; i < nodeChild.size(); ++i) {
        mem += nodeChild[i].memoryCost();
    }
    return mem;
}

void TrieNode::addMatch(VertexID *match, const VertexID *order, ui matchSize, ui startPos, VertexID **data, VertexID *length,
                        ui num) {
    TrieNode &current = *this;
    while (true) {
        if (startPos + num == matchSize) {
            break;
        } else {
            VertexID v = match[order[startPos]];
            if (current.nodeChild.empty() || current.nodeChild.back().value < v) {
                current.nodeChild.emplace_back(v);
            }
            current = current.nodeChild.back();
            startPos += 1; // Prepare for the next iteration
        }
    }
}

void TrieNode::addMatch(const std::vector<VertexID> &match, ui startPos) {
    TrieNode &current = *this;
    while (startPos < match.size()) {
        VertexID v = match[startPos];
        if (current.nodeChild.empty() || current.nodeChild.back().value < v) {
            current.nodeChild.emplace_back(v);
        }
        current = current.nodeChild.back();
        startPos++; // Move to the next position
    }
}

void TrieNode::buildTrieFromSortedMatchesBatch(
    const std::vector<std::vector<VertexID>> &matches,
    ui start, ui end, ui depth) {

    if (start >= end) return;
    if (depth == matches[start].size() - 1) {
        nodeChild.resize(end - start);
        for (ui idx = start; idx < end; ++idx)
            nodeChild[idx - start].value = matches[idx][depth];
        return;
    }

    std::vector<std::tuple<VertexID, ui, ui>> groups;
    ui i = start;

    while (i < end) {
        VertexID currentKey = matches[i][depth];
        ui groupStart = i;

        while (i < end && matches[i][depth] == currentKey) {
            ++i;
        }

        groups.emplace_back(currentKey, groupStart, i);
    }

    ui groupCount = groups.size();
    nodeChild.resize(groupCount);

    for (ui idx = 0; idx < groupCount; ++idx) {
        auto [key, groupStart, groupEnd] = groups[idx];
        nodeChild[idx].value = key;
        nodeChild[idx].buildTrieFromSortedMatchesBatch(matches, groupStart, groupEnd, depth + 1);
    }
}


size_t TrieNode::numTuples(bool *visited) const {
    if (nodeChild.empty()) return 1;
    size_t num = 0;
    for (int i = 0; i < nodeChild.size(); ++i) {
        if (!visited[nodeChild[i].value])
            num += nodeChild[i].numTuples(visited);
    }
    return num;
}

size_t TrieNode::numResults(VertexID *partMatch, const std::vector<VertexID> &toCheck) const {
    ui beginPos = 0, endPos = nodeChild.size();
    size_t num = endPos - beginPos;
    for (VertexID u: toCheck) {
        VertexID v = partMatch[u];
        if (checkExists(nodeChild.data(), beginPos, endPos, v) != endPos)
            --num;
    }

    return num;
}

size_t TrieNode::numResults(VertexID *partMatch, const std::vector<VertexID> &toCheck,
                            const std::vector<VertexID> &largerAttrs, const std::vector<VertexID> &smallerAttrs) const {
    VertexID maxCompared = 0, minCompared = (VertexID) - 1;
    for (VertexID u: largerAttrs) {
        if (partMatch[u] > maxCompared)
            maxCompared = partMatch[u];
    }
    for (VertexID u: smallerAttrs) {
        if (partMatch[u] < minCompared)
            minCompared = partMatch[u];
    }
    ui beginPos = 0, endPos = nodeChild.size();
    if (!largerAttrs.empty()) {
        beginPos = upperBound(nodeChild.data(), 0, nodeChild.size(), maxCompared);
    }
    if (!smallerAttrs.empty()) {
        endPos = lowerBound(nodeChild.data(), 0, nodeChild.size(), minCompared);
    }
    size_t num = endPos - beginPos;
    for (VertexID u: toCheck) {
        VertexID v = partMatch[u];
        if (checkExists(nodeChild.data(), beginPos, endPos, v) != endPos)
            --num;
    }

    return num;
}

int binarySearch(const std::vector<TrieNode> &array, VertexID v) {
    if (array.size() == 0) return -1;
    int left = 0;
    int right = array.size() - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (array[mid].value == v) {
            return mid;
        } else if (array[mid].value < v) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return -1; // Element not found
}

ui binarySearch(const TrieNode *child, const ui begin, const ui end, const ui target) {
    ui offset_begin = begin;
    ui offset_end = end;
    while (offset_end - offset_begin >= BINARY_SEARCH_THRESHOLD) {
        auto mid = (offset_begin + offset_end) / 2;
#ifndef __APPLE__
        _mm_prefetch((char *) &child[(mid + 1 + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &child[(mid - 1 + offset_begin) / 2], _MM_HINT_T0);
#endif
        if (child[mid].value == target) {
            return mid;
        } else if (child[mid].value < target) {
            offset_begin = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback
    for (auto offset = offset_begin; offset < offset_end; ++offset) {
        if (child[offset].value >= target) {
            return offset;
        }
    }

    return offset_end;
}

ui checkExists(const TrieNode *child, const ui begin, const ui end, const ui target) {
    ui offset_begin = begin;
    ui offset_end = end;
    while (offset_end - offset_begin >= BINARY_SEARCH_THRESHOLD) {
        auto mid = (offset_begin + offset_end) / 2;
#ifndef __APPLE__
        _mm_prefetch((char *) &child[(mid + 1 + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &child[(mid - 1 + offset_begin) / 2], _MM_HINT_T0);
#endif
        if (child[mid].value == target) {
            return mid;
        } else if (child[mid].value < target) {
            offset_begin = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback
    for (auto offset = offset_begin; offset < offset_end; ++offset) {
        if (child[offset].value == target) {
            return offset;
        }
    }

    return end;
}

ui checkExists(const VertexID *array, const ui begin, const ui end, const ui target) {
    ui offset_begin = begin;
    ui offset_end = end;
    while (offset_end - offset_begin >= BINARY_SEARCH_THRESHOLD) {
        auto mid = (offset_begin + offset_end) / 2;
#ifndef __APPLE__
        _mm_prefetch((char *) &array[(mid + 1 + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &array[(mid - 1 + offset_begin) / 2], _MM_HINT_T0);
#endif
        if (array[mid] == target) {
            return mid;
        } else if (array[mid] < target) {
            offset_begin = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback
    for (auto offset = offset_begin; offset < offset_end; ++offset) {
        if (array[offset] == target) {
            return offset;
        }
    }

    return end;
}

ui upperBound(const TrieNode *array, ui begin, ui end, VertexID v) {
    while (end - begin >= BINARY_SEARCH_THRESHOLD) {
        int mid = begin + (end - begin) / 2;
        if (array[mid].value <= v) {
            begin = mid + 1;
        } else {
            end = mid;
        }
    }

    // Linear search fallback
    for (int i = begin; i < end; ++i) {
        if (array[i].value > v) {
            return i;
        }
    }
    return end;
}

ui lowerBound(const TrieNode *array, ui begin, ui end, VertexID v) {
    while (end - begin >= BINARY_SEARCH_THRESHOLD) {
        int mid = begin + (end - begin) / 2;
        if (array[mid].value >= v) {
            end = mid;
        } else {
            begin = mid + 1;
        }
    }

    // Linear search fallback
    for (int i = begin; i < end; ++i) {
        if (array[i].value >= v) {
            return i;
        }
    }
    return end;
}

size_t TrieLevel::memoryCost() const {
    size_t mem = sizeof(*this);
    for (int i = 0; i < length; ++i) {
        mem += children[i].memoryCost();
    }
    return mem;
}

size_t TrieLevel::numTuples(bool *visited) const {
    size_t num = 0;
    for (int i = 0; i < length; ++i) {
        if (!visited[values[i]]) ++num;
        num += children[i].numTuples(visited);
    }

    return num;
}

size_t TrieLevel::numResults(VertexID *partMatch, const std::vector<VertexID> &toCheck) const {
    ui beginPos = 0, endPos = length;
    size_t num = endPos - beginPos;
    for (VertexID u: toCheck) {
        VertexID v = partMatch[u];
        if (checkExists(values, beginPos, endPos, v) != endPos)
            --num;
    }
    return num;
}

size_t TrieLevel::numResults(VertexID *partMatch, const std::vector<VertexID> &toCheck,
                             const std::vector<VertexID> &largerAttrs,
                             const std::vector<VertexID> &smallerAttrs) const {
    VertexID maxCompared = 0, minCompared = (VertexID) - 1;
    for (int i = 0; i < largerAttrs.size(); ++i) {
        if (partMatch[largerAttrs[i]] > maxCompared)
            maxCompared = partMatch[largerAttrs[i]];
    }
    for (int i = 0; i < smallerAttrs.size(); ++i) {
        if (partMatch[smallerAttrs[i]] < minCompared)
            minCompared = partMatch[smallerAttrs[i]];
    }
    ui beginPos = 0, endPos = length;
    if (!largerAttrs.empty()) {
        beginPos = setBeginPos(values, length, maxCompared);
    }
    if (!smallerAttrs.empty()) {
        endPos = setEndPos(values, length, minCompared);
    }
    size_t num = endPos - beginPos;
    for (VertexID u: toCheck) {
        VertexID v = partMatch[u];
        if (checkExists(values, beginPos, endPos, v) != endPos)
            --num;
    }
//    for (int i = beginPos; i < endPos; ++i) {
//        bool exists = false;
//        for (VertexID u: toCheck) {
//            if (values[i] == partMatch[u])
//                exists = true;
//        }
//        if (!exists) {
//            partMatch[1] = values[i];
//            gResult1.emplace_back(partMatch, partMatch + 6);
//        }
//    }
    return num;
}

void TrieLevel::buildTrieFromSortedMatchesBatch(const std::vector<std::vector<VertexID>> &matches, ui start, ui end,
                                                ui depth) {
    if (matches.empty() || start >= end) {
        length = 0;
        return;
    }
    if (depth == matches[start].size() - 1) {
        length = end - start;
        delete[] values;
        values = new VertexID[length];
        for (ui idx = start; idx < end; ++idx)
            values[idx - start] = matches[idx][depth];
        return;
    }
    std::vector<std::tuple<VertexID, ui, ui>> groups;
    ui i = start;

    while (i < end) {
        VertexID currentKey = matches[i][depth];
        ui groupStart = i;

        while (i < end && matches[i][depth] == currentKey) {
            ++i;
        }

        groups.emplace_back(currentKey, groupStart, i);
    }

    length = groups.size();
    delete[] values;
    values = new VertexID[length];

    delete[] children;
    children = new TrieLevel[length];
    for (ui idx = 0; idx < length; ++idx) {
        auto [key, groupStart, groupEnd] = groups[idx];
        values[idx] = key;
        children[idx].buildTrieFromSortedMatchesBatch(matches, groupStart, groupEnd, depth + 1);
    }
}
