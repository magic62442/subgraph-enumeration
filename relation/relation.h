//
// Created by anonymous authors on 2024/2/27.
//

#ifndef ASDMATCH_RELATION_H
#define ASDMATCH_RELATION_H

#include "compute_set_intersection.h"
#include <vector>

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

ui checkExists(const VertexID *array, const ui begin, const ui end, const ui target);

#endif //ASDMATCH_RELATION_H
