#' Generates a random tree/forest (depending on the value of `blocks`) and calls `generate_benchmark_from_graph` on that.
#' @param n A vector of integer numbers containing the number of vertices each component has.
#' @param seed A random seed to set before generating the random base forest.
#' @param blocks The number of components to generate. Defaults to `length(n)`.
#' @export
generate_benchmark_random <- function(n, seed = as.integer(runif(1, max = 10000)), blocks = length(n)) {
    set.seed(seed)
    n <- rep_len(n, blocks)

    base_g <- igraph::make_empty_graph(directed = TRUE)

    for (i in seq_along(n)) {
        base_g_component <- igraph::make_full_graph(n = n[i], directed = TRUE)
        base_g <- igraph::disjoint_union(base_g, igraph::subgraph.edges(base_g_component, igraph::sample_spanning_tree(base_g_component)))
    }

    generate_benchmark_from_graph(base_g)
}

#' Expand a directed base graph (tree or forest) into a bipartite graph with overlapping nested structure. An `a -> b` directed edge means the neighborhood of `a` is a subset of the neighborhood of `b` (i.e. `a` will be nested with `b`).
#' It is suggested to use `generate_benchmark_random` instead of directly calling this function.
#' @param base_g An `igraph` graph object. Should be directed and acyclic. Does not need to have only one component.
#' @return A bipartite graph where the nodes of `base_g` are nested in the way their edges are directed. This can create overlapping nested communities.
#' @export
generate_benchmark_from_graph <- function(base_g) {
    n <- igraph::components(base_g)$csize

    add_neighs <- lapply(seq_len(sum(n)), function(i) c())
    nodes_added <- 0L

    block_sizes <- c(0, cumsum(n))
    nodes_visited <- rep_len(FALSE, length.out = sum(n))

    neighs_out <- igraph::adjacent_vertices(base_g, igraph::V(base_g), mode = "out")
    neighs_in <- igraph::adjacent_vertices(base_g, igraph::V(base_g), mode = "in")

    # we generate the nested neighborhood for all blocks, separately
    for (i in seq_along(n)) {
        # start with the first vertex that has no in-neighbors
        current_vertex <- which.min(vapply(neighs_in, function(v) length(v) == 0, logical(1)))

        logger::log_debug("Starting block {i} with vertex {current_vertex}")

        # while we haven't visited all the vertices in the block, do the generation
        while (!all(nodes_visited[seq(from = block_sizes[i] + 1, to = block_sizes[i + 1])])) {
            in_current <- neighs_in[[current_vertex]]
            out_current <- neighs_out[[current_vertex]]

            # if all in-neighbors are visited, we can proceed
            if (length(in_current) == 0 || all(nodes_visited[in_current])) {
                # if this node was not visited before, generate neighbors for it and jump to the next vertex
                if (!nodes_visited[current_vertex]) {
                    if (length(in_current) > 0) {
                        add_neighs[[current_vertex]] <- c(add_neighs[[current_vertex]], unlist(add_neighs[in_current]))
                    }

                    # add an extra neighbor
                    add_neighs[[current_vertex]] <- c(add_neighs[[current_vertex]], nodes_added + 1L)

                    # mark as visited and increase the added node count
                    nodes_visited[current_vertex] <- TRUE
                    nodes_added <- nodes_added + 1L
                    logger::log_debug("✓ vertex {current_vertex} [{paste0(add_neighs[[current_vertex]], collapse = ', ')}]")
                }

                # the next vertex to visit
                current_vertex <- if (length(out_current) > 0 && !all(nodes_visited[out_current])) {
                    # pick an out-neighbor that hasn't been visited before
                    logger::log_debug("Picked outbound neighbor {current_vertex} -> {out_current[which.min(nodes_visited[out_current])]} (choices: [{paste0(nodes_visited[out_current], collapse = ', ')}])")
                    out_current[which.min(nodes_visited[out_current])]
                } else {
                    # otherwise check if there are any in-neighbors
                    if (length(in_current) > 0) {
                        # pick the first in-neighbor
                        logger::log_debug("Picked in-neighbor {current_vertex} <- {in_current[which.min(nodes_visited[in_current])]} (choices: [{paste0(!nodes_visited[in_current], collapse = ', ')}])")
                        in_current[which.min(nodes_visited[in_current])]
                    } else {
                        # if there are no unvisited out-neighbors and no in-neighbors, pick the first unvisited node
                        logger::log_debug("Picked unvisited node {current_vertex} ~> {which.min(nodes_visited)}")
                        which.min(nodes_visited)
                    }
                }
            } else {
                # if there are unvisited in-neighbors, go visit them first
                logger::log_debug("Picked an in_neighbor (not all are visited): {current_vertex} <- {in_current[which.min(nodes_visited[in_current])]} (choices: [{paste0(nodes_visited[in_current], collapse = ', ')}])")

                current_vertex <- in_current[which.min(nodes_visited[in_current])]
            }
        }
    }

    # add the neighbors
    final_g <- igraph::make_empty_graph(n = 2 * sum(n), directed = FALSE)
    for (i in seq_along(add_neighs)) {
        neighs_to_add <- unique.default(add_neighs[[i]])
        final_g <- igraph::add_edges(final_g, as.vector(t(matrix(
            data = c(
                rep_len(i, length(neighs_to_add)),
                neighs_to_add + sum(n)
            ),
            ncol = 2
        ))))
    }

    # make the final graph bipartite
    final_g_type <- rep_len(FALSE, igraph::vcount(final_g))
    final_g_type[seq_len(igraph::vcount(base_g))] <- TRUE
    final_g <- igraph::set_vertex_attr(
        graph = final_g,
        name = "type",
        value = final_g_type
    )

    # extract nested communities
    leaves <- which(igraph::degree(base_g, mode = "out") == 0)
    truth <- lapply(
        which(igraph::degree(base_g, mode = "in") == 0),
        function(i) {
            lapply(igraph::all_simple_paths(base_g, from = i, to = leaves), igraph::as_ids)
        }
    )

    truth <- unlist(truth, recursive = FALSE)

    return(list(
        final_g = final_g,
        base_g = base_g,
        truth = truth
    ))
}