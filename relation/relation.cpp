//
// Created by anonymous authors on 2024/2/27.
//

#include "relation.h"
#include "utils.h"
#include <tuple>

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
