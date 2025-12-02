# Subgraph Enumeration

## Background

This project provides efficient algorithms for enumerating subgraphs in large networks. It is based on the new minimal fractional hypertree decompositions, a new application of symmetry-breaking, a new cost model to select attribute orders, and a new join algorithm. It can handle both single and multi-threaded execution.

## Compile

### Prerequisites

This project combines dependencies from both ASDMatch and SCOPE projects:

1. **Compile and link to the [nauty](https://pallini.di.uniroma1.it) library** (for automorphism computation):

```shell
cd utility/automorphism/
./configure
make
mv nauty.a libnauty.a
```

If you encounter relocation errors, add `-fPIC` flag:

```shell
cd utility/automorphism/
vim makefile
# add -fPIC to the end of line 6
make
mv nauty.a libnauty.a
```

2. **Compile and link to the [GLPK](https://www.gnu.org/software/glpk/) library** (for fractional hypertree decomposition):

```shell
cd utility/td/glpk-4.55
./configure --prefix={ABSOLUTE_PATH_TO_GLPK_LIB}  # relative path: ../glpk_lib
make -j
make install
cd .. # back to 'td' directory
```

3. **Build the project**:

```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```

## Input format

The pattern graph and data graph should start with 'n, m' where n is the number of nodes and m is the number of undirected edges, followed by the edge list. The node id should be consecutive and should start from 0. Each edge only appears once, in the format of 'smaller_vertex_id larger_vertex_id'. The separator is a ' '.

Example:

```
3 2
0 1
1 2
```

## Execution and Output

### Main Executable (graph_mining_system):

| Option | Description |
|--------|-------------|
| -q | Query graph path |
| -d | Data graph path |
| -r | Result path (optional) |
| -i | Enable exploration-style intersection. The default is WCOJ-style. (optinal) |
| -u | Enable IEP (Intersection-based Execution), optinal. |
| -t | Number of threads for parallel execution (default value is 1, optinal) |
| --help | Show help information |

Example:

```shell
./graph_mining_system -q query.graph -d data.graph -r results.txt -t 4
```

### Specialized Pattern Executables:

The project includes optimized executables for specific patterns  to compare with hand-optimized algorithms:

- `diamond.out`
- `house.out`
- `triple_triangle.out`
- `quad_triangle.out`
- `four_square.out`
- `near_five_clique.out` 

### Utility Tools:

- `txt2bin.out` - Convert text graph format to binary format for faster loading

## Performance Notes

- Enable exploration-style intersection for `yt` and `up`graphs
- Use appropriate thread count (-t) based on your system cores
- Consider using pattern-specific executables for better performance on known patterns

## Dataset

The project includes sample datasets in the `dataset/` directory:
- `light_query/` - Small query patterns for testing
- `huge_query/` - Large query patterns for benchmarking
- `6vertex/` - 6-vertex query patterns. The experiments use 16 patterns with 11-13 edges.
- `7vertex/` - 7-vertex query patterns. The experiments use 78 patterns with 16-19 edges.
