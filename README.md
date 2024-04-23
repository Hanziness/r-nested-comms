# Nested Community Detection Algorithms
This R package contains algorithms to do hierarchical clustering for nestedness and find overlapping fully nested subgraphs in both bipartite and non-bipartite networks. This version of the package (extended with hierarchical clustering algorithms) is part of the work submitted in the following article:

```
TBA
```

## Installation
Make sure you have tools to build C++ source files. Download the package's source code (a zip file) and use the `devtools` R package to install it:

```{r}
devtools::install_local("sourcecode.zip")
```

This will build the package and install it as `nested.comms`.

**Note.** The community detection algorithm is written in C++ using containers from `std`. It uses an edge-list to perform its operations, then returns a community list. This makes it possible to use the algorithm without R, without requiring too many modifications (remove the `Rcpp` wrapper calls and return structures from the function `get_lk_all_topsort` in `src/node_comm.cpp`).

## Usage
The package operates on `igraph` networks, so make sure to prepare your networks as `igraph` graphs. Vertices may be named, edge weights and directionality are ignored.

### Hierarchical clustering (new in v0.3)
To perform hierarchical clustering, one can use two approaches. The _bottom-up_ method is done using R's built-in `hclust` method, wrapped by `bottom_up_clustering`. This returns a list of cluster memberships, one vector for each level:

```{r}
results <- nested.comms::bottom_up_clustering(graph, method = "complete", level_per_node = TRUE)
```

_Top-down_ hierarchical clustering is implemented as a modified version of the Girvan-Newman algorithm. In order to perform top-down clustering, we first create a _nestedness graph_ using the function `nestedness_graph` (this returns a full graph where the edge weights are the nestedness values), and then run `girvan_newman_nestedness` on it. The result here is once again a list of cluster memberships, each vector corresponding to a clustering level.

```{r}
nested_g <- nested.comms::nestedness_graph(g)
results <- nested.comms::girvan_newman_nestedness(nested_g, method = "full")
```

For more information, check the documentation of each function!

### Community detection
To perform community detection on an `igraph` network, use the `nested_node_comms` function:

```{r}
graph <- igraph::read_graph("network.gml", format = "gml")
communities <- nested.comms::nested_node_comms(graph)
```

The result will be a list of vectors. Each item (vector) is a fully nested community, where the vertices are in ascending order of their neighborhood - ie. the neighborhood of `comm[i]` is a subset of the neighborhood of `comm[i + 1]`. Vertices that are not nested with any other vertex are put into their own community.

### Building a community graph
Using the results of `nested_node_comms` one can build the community graph of the network with the `nested_node_community_graph` function.

```{r}
community_graph <- nested.comms::nested_node_community_graph(communities)
```

If one wants the bidirectional edges between nodes that have equal neighborhoods, the original network must be supplied, too:

```{r}
community_graph <- nested.comms::nested_node_community_graph(communities, graph)
```

The resulting graph may contain multiple edges between a given pair of vertices (an edge is added for each community). If only a single edge is desired, the graph can be simplified afterwards:

```{r}
igraph::simplify(community_graph, remove.multiple = TRUE)
```

### Generating test graphs
To generate graphs with a ground-truth overlapping nested community structure, use `generate_benchmark_random`:

```{r}
random_nestedness_graph <- nested.comms::generate_benchmark_random(c(4, 6))
```

This generates a bipartite graph with two components (one with 4 nodes, and another with 6). The result is a list with three items:

* `final_g`: the bipartite graph with overlapping nested communities
* `base_g`: the underlying random tree (forest) from which `final_g` was generated. This is essentially the community graph of the network.
* `truth`: the nested communities in the same structure as with `nested_node_comms`.

**Note** that the resulting graph has `2 * sum(n)` vertices, the ground-truth community structure only applies to the first `sum(n)` vertices!