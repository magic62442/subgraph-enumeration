//
// Created by Qiyan on 25-3-10.
//

#ifndef GRAPH_MINING_SYSTEM_CONFIG_H
#define GRAPH_MINING_SYSTEM_CONFIG_H

#endif //GRAPH_MINING_SYSTEM_CONFIG_H

#include <cstdlib>
#include <utility>
#include <cstdint>

/**
 * Set intersection method.
 * 0: Hybrid method; 1: Merge based set intersections.
 */
#define HYBRID 0

/**
 * Accelerate set intersection with SIMD instructions.
 * 0: AVX2; 1: AVX512; 2: Basic;
 * If building in mac AMD chip, can only use 2.
 * For machines using Intel CPU, use 0.
 */

#ifdef __APPLE__
#define SI 2
#else
#define SI 0
#endif

#define MAX_PATTERN_SIZE 32
#define SAMPLE_SIZE 40000
#define SAMPLE_PORTION 4
#define COEEFICIENT 20.0
#define BINARY_SEARCH_THRESHOLD 32
//#define COLLECT_STATISTICS
//#define ALL_LEVEL
//#define COLLECT_RESULT

typedef uint32_t ui;
typedef uint32_t VertexID;
typedef uint32_t EdgeID;
typedef uint32_t LabelID;
typedef uint64_t Count;
typedef Count *HashTable;
typedef uint64_t CanonType;
typedef std::pair<VertexID, VertexID> Edge;
