#' Find maximal sets from the results of `nested_node_comm`
#' @param list_of_lk The result list of `nested_node_comm`
#' @return A list of maximal sets
#' @export
find_maximal_sets <- function(list_of_lk) {
    found_sets <- list()
    n <- length(list_of_lk)

    for (i in rev(seq_len(n))) {
        current_lk <- list_of_lk[[i]]
        for (row_index in seq_along(current_lk)) {
            current_row <- current_lk[[row_index]]
            if (!is_in_list(
                current_row,
                found_sets
            )) {
                found_sets <- c(found_sets, list(current_row))
            }
        }
    }

    found_sets
}

#' Shorthand function to check if a set part of any of the list of sets
is_in_list <- function(set, list_of_sets) {
    any(vapply(
        list_of_sets,
        function(check_set) {
            all(set %in% check_set)
        },
        FUN.VALUE = logical(1)
    ))
}

#' Build the community graph using the resulting community structure of `nested_node_comms`.
#' @param results The community structure obtained using `nested_node_comms`
#' @param sameneigh_original_g (Optional) The original graph `nested_node_comms` was performed on. If given, bidirectional edges will be added between vertices with the same neighborhood.
#' @return An `igraph` graph of the community graph
#' @export
nested_node_community_graph <- function(results, sameneigh_original_g = NULL) {
    vertices <- results |>
        unlist() |>
        unique.default()

    n <- length(vertices)

    graph <- igraph::make_empty_graph(n = n, directed = TRUE)
    # graph <- igraph::set_vertex_attr(graph, name = "name", value = vertices)

    edgelist <- c()
    neighs <- if (!is.null(sameneigh_original_g))
        lapply(igraph::adjacent_vertices(sameneigh_original_g, v = igraph::V(sameneigh_original_g), mode = "out"), function(neighs) {
            sort(as.numeric(neighs))
        }) else NULL

    for (i in seq_along(results)) {
        item <- results[[i]]
        if (length(item) >= 2) {
            # repeat all items in the path, except the first and last ones
            # e.g. 1 2 3 4 -> 1 2 2 3 3 4
            repeatvector <- rep_len(2L, length(item))
            repeatvector[c(1, length(repeatvector))] <- 1L
            new_edges <- unname(unlist(mapply(rep, item, repeatvector)))
            edgelist <- c(edgelist, new_edges)
            graph <- igraph::add_edges(graph, new_edges, community = i)
        }

        if (length(item) >= 2 && !is.null(sameneigh_original_g)) {
            for (j in seq_len(length(item) - 1)) {
                if (identical(neighs[[item[j]]], neighs[[item[j + 1]]])) {
                    graph <- igraph::add_edges(graph, c(item[j + 1], item[j]), community = i)
                }
            }
        }
    }

    # graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)

    return(graph)
}

#' Build nested communities in a hierarchical fashion
#' @param g The input graph (an `igraph` graph object)
#' @return A list of `data.table` objects where list item `k` corresponds to nested sets of length `k`
#' @import foreach
#' @export
nested_node_comm <- function(g, use.names = FALSE) {
    neighs <- lapply(
        igraph::adjacent_vertices(g, v = igraph::V(g), mode = "out"),
        as.numeric
    )

    results <- get_lk_all_topsort(neighs)

    # Add back vertex names if needed
    if (use.names) {
        vertex_names <- igraph::as_ids(igraph::V(g))
        results <- lapply(results, function(community) {
            vertex_names[community]
        })
    }

    return(results)
}