//
// Created by anonymous authors on 2024/2/27.
//

#ifndef ASDMATCH_RELATION_H
#define ASDMATCH_RELATION_H

#include "labeled_graph.h"
#include "compute_set_intersection.h"

extern std::vector<std::vector<VertexID>> gResult1;

struct TrieNode {
    VertexID value;
    std::vector<TrieNode> nodeChild;
//    Compression *compression;

//    TrieNode(): value(0), nodeChild(), compression(nullptr) {}
    TrieNode(): value(0), nodeChild() {}
    TrieNode(VertexID value): value(value), nodeChild() {}
//    TrieNode(const TrieNode& other) = delete;
//    TrieNode& operator=(const TrieNode& other) = delete;
    ~TrieNode() {
        nodeChild.clear();
    }

    size_t memoryCost() const;
    size_t numTuples(bool *visited) const;
    size_t numResults(VertexID *partMatch, const std::vector<VertexID> &toCheck) const;
    size_t numResults(VertexID *partMatch, const std::vector<VertexID> &toCheck,
                      const std::vector<VertexID> &largerAttrs, const std::vector<VertexID> &smallerAttrs) const;
    void addMatch(VertexID *match, const VertexID *order, ui matchSize, ui startPos, VertexID **data, VertexID *length,
                  ui num);
    void addMatch(const std::vector<VertexID> &match, ui startPos);
    void buildTrieFromSortedMatchesBatch(const std::vector<std::vector<VertexID>> &matches,
                                    ui start, ui end, ui depth);
    friend bool operator==(const TrieNode& lhs, const TrieNode& rhs) {
        if (lhs.value != rhs.value) return false;
        // If no compression, compare children
        if (lhs.nodeChild.size() != rhs.nodeChild.size()) return false;
        for (int i = 0; i < lhs.nodeChild.size(); ++i) {
            if (!(lhs.nodeChild[i] == rhs.nodeChild[i])) return false;
        }

        return true;
    }
};

struct TrieLevel {
    VertexID *values;
    size_t length;
    bool oneLevel;
    TrieLevel *children;
    TrieLevel(): values(nullptr), length(0), oneLevel(false), children(nullptr) {}
    ~TrieLevel() {
        if (!oneLevel) {
            delete[] values;
            delete[] children;
        }
    }

    size_t memoryCost() const;
    size_t numTuples(bool *visited) const;
    size_t numResults(VertexID *partMatch, const std::vector<VertexID> &toCheck) const;
    size_t numResults(VertexID *partMatch, const std::vector<VertexID> &toCheck,
                      const std::vector<VertexID> &largerAttrs, const std::vector<VertexID> &smallerAttrs) const;
    void buildTrieFromSortedMatchesBatch(const std::vector<std::vector<VertexID>> &matches,
                                         ui start, ui end, ui depth);
};

int binarySearch(const std::vector<TrieNode> &array, VertexID v);
ui binarySearch(const TrieNode *child, const ui begin, const ui end, const ui target);
ui checkExists(const TrieNode *child, const ui begin, const ui end, const ui target);
ui checkExists(const VertexID *array, const ui begin, const ui end, const ui target);
ui upperBound(const TrieNode *array, ui begin, ui end, VertexID v);
ui lowerBound(const TrieNode *array, ui begin, ui end, VertexID v);

#endif //ASDMATCH_RELATION_H
