# .datatable.aware <- TRUE

#' Edge-based nestedness community detection
#' @param g The input `igraph` graph object
#' @param mode Controsl how many checks the algorithm will perform before
#' returning the communities. Might result in better edge community detection
#' but will also require more iterations to complete. Might find some
#' false communities that are subsets of others, these can be removed using
#' `post_process = TRUE`. The three accepted values are `fast` (n1 > 0),
#' `precise` (n1 >= 0) and `strict` (n2 > 0), the last one comparing all
#' accepted node pairings.
#' @param post_process Whether to remove communities that are subsets
#' of other communities. Useful when `mode = "precise"` or `mode = "strict"`.
#' @param n1_threshold The threshold parameter to use in the iteration
#' condition. `1` is the "fast" mode, "0" is the "precise" mode and `-Inf` is
#' the "strict" mode.
#' @import data.table
#' @import igraph
#' @import logger
#' @export
nested_edge_comm <- function(
    g, mode = c("fast", "precise", "strict"), n1_threshold = 1L,
    post_process = TRUE,
    debug_steps = FALSE
) {
    mode <- match.arg(mode)

    debug_log <- NULL
    if (debug_steps) {
        debug_log <- list()
        logger::log_info("Saving iterations result into `debug_steps` variable")
    }

    # 0. Set edge IDs as attributes for easier referencing

    g <- igraph::set_edge_attr(
        graph = g, name = "id", value = seq_len(igraph::ecount(g))
    )

    # Which edge belongs to which communities
    edge_communities <- vector(mode = "list", length = igraph::ecount(g))
    n <- igraph::vcount(g)

    # 1-neighbors of all vertices
    neighs1 <- igraph::neighborhood(
        graph = g,
        nodes = seq_len(n),
        order = 1,
        mode = "out",
        mindist = 1
    )

    # 2-neighbors (nodes with distance 2) = the N++ set
    neighs2 <- igraph::neighborhood(
        graph = g,
        nodes = seq_len(n),
        order = 2,
        mode = "out",
        mindist = 2
    )

    # Neighborhood count data frame: `n1` is the number of adjacent vertices,
    # `n2` is the number of N++ vertices belonging to each node.
    neighs_df <- data.table::data.table(
        v = seq_len(n),
        n1 = vapply(
            X = neighs1,
            FUN = length,
            FUN.VALUE = numeric(1)
        ),
        n2 = vapply(
            X = neighs2,
            FUN = length,
            FUN.VALUE = numeric(1)
        )
    )

    # Current largest assigned color
    current_max_color <- 0

    # which node is incident with which edges
    node_incident_edges <- lapply(
        igraph::incident_edges(
            g, seq_len(n),
            mode = "out"
        ),
        FUN = function(el) {
            # return only the edge IDs
            el$id
        }
    )
    excluded_pairs <- list()

    # TODO implement good node pair selection heuristic

    # Iterate:
    # 1. Take a node
    #   a) That has at least one incident edge that is not marked (n1 > 0)
    #   b) That has the highest 2-neighbor count
    # 2. Take its 2-neighbour with the most unclassified edges (max(n1))
    # 3. Mark their edges with either the same or a different colour
    # based on whether they are nested or not
    #   This includes additional steps to exclude colours that might merge
    #   two non-nested nodes into a single community.
    # 4. Keep going until this node has unclassified edges

    # Number of iterations performed
    iter <- 0

    if (debug_steps) {
        debug_log[[1]] <- list(
            neighs_df = data.table::copy(neighs_df)
        )
    }

    # n1_threshold <- 1
    # iteration_condition <- quote(any(neighs_df$n1 >= 1))
    # if (mode == "precise") {
    #     iteration_condition <- quote(any(neighs_df$n1 >= 0))
    #     n1_threshold <- 0
    # } else if (mode == "strict") {
    #     iteration_condition <- quote(any(neighs_df$n2 >= 1))
    #     n1_threshold <- 1
    # }

    iteration_condition <- quote(
        nrow(neighs_df[n2 >= 1 & n1 >= n1_threshold]) > 0
    )

    while (eval(iteration_condition)) {
        # first node with available edges and maximal n2 count
        # node_a <- ifelse(
        #     mode == "strict",
        #     neighs_df[n2 == max(n2), v][1],
        #     neighs_df[n1 >= n1_threshold][n2 == max(n2), v][1]
        # )
        node_a <- neighs_df[n2 >= 1 & n1 >= n1_threshold][n2 == max(n2), v][1]
        edges_a <- igraph::E(g)[.from(node_a)]$id # nolint
        colour_a <- NULL

        node_b <- NULL
        edges_b <- NULL
        colour_b <- NULL

        logger::log_debug("# {iter + 1}: {node_a}")

        # if node has no 2-neighbors, just give its edges a new colour
        if (neighs_df[node_a, n2] == 0) {
            logger::log_debug("-- NO 2-NEIGH")
            # if all edges are already coloured, don't assign another colour
            if (!all(vapply(
                edge_communities[edges_a],
                function (comm) length(comm) > 0,
                logical(1)
            ))) {
                # new colour for `node_a`'s edges
                colour_a <- c(current_max_color + 1)
                current_max_color <- current_max_color + 1

                logger::log_debug("-> NEW COLOUR")
            } else {
                # All edges had colours, skip assigning
                logger::log_debug("-> SKIP")
            }
        } else {
            # node that is 2-neighbour of `node_a` and has maximal
            # 1-neighbour count
            node_b <- neighs_df[
                setdiff(neighs2[[node_a]], excluded_pairs[[as.character(node_a)]])
            ][n1 == max(n1), v][1]
            edges_b <- igraph::E(g)[.from(node_b)]$id # nolint

            # exclude the (a, b) pair
            excluded_pairs[[as.character(node_a)]] <-
                c(excluded_pairs[[as.character(node_a)]], node_b)
            excluded_pairs[[as.character(node_b)]] <-
                c(excluded_pairs[[as.character(node_b)]], node_a)

            logger::log_debug("-- &: {node_b}")

            # Whether nodes `a` and `b` are nested
            is_nested <- !is_2k2_neigh(
                neighs1[[node_a]],
                neighs1[[node_b]]
            )

            # TODO debug message
            logger::log_debug("--  = {is_nested}")

            if (is_nested) {
                # this is not enough to assign same colours
                # find 2-neighbors with same colours and check nestedness
                # against them, too!
                # this allows us to avoid cases where two nodes would be given
                # a same colour and a third one in the same community would
                # not be nested with `node_a` or `node_b` as a result.

                excluded_colours <- find_excluded_colors(
                    node_a, node_b, edges_a, edges_b,
                    node_communities, edge_communities, 
                    neighs1, neighs2
                )

                if (length(excluded_colours) > 0) {
                    logger::log_debug("-- EXCL ({paste0(
                            excluded_colours,
                            collapse = \" \"
                    )})")
                }

                # assign the same color to both of them
                newcolor <- find_same_color(
                    edge_communities[edges_a], edge_communities[edges_b],
                    current_max_color,
                    excluded_colours = excluded_colours
                )

                # Update current_max_color if a higher color was assigned
                if (max(newcolor) > current_max_color) {
                    current_max_color <- max(newcolor)
                }

                colour_a <- colour_b <- newcolor

                # TODO debug message
                logger::log_debug("-- UPDATE E {paste0(c(edges_a, edges_b), collapse = \" \")} -> {newcolor}")
            } else {
                # Nodes are not nested, assign their edges to different
                # communities


                # TODO do this smarter (like in the case of nested edges):
                # find a colour that can be assigned to the edges by excluding
                # the other node's colours + 2K2 causing colours
                # check if all "a" edges have a colour already
                if (!all(vapply(
                    edge_communities[edges_a],
                    function(x) length(x) > 0,
                    logical(1)
                ))) {
                    colour_a <- c(current_max_color + 1)
                    current_max_color <- current_max_color + 1
                }

                # and the same for all the "b" edges
                if (!all(vapply(
                    edge_communities[edges_b],
                    function(x) length(x) > 0,
                    logical(1)
                ))) {
                    colour_b <- c(current_max_color + 1)
                    current_max_color <- current_max_color + 1
                }

                # TODO debug message
                logger::log_debug("-- UPDATE E {paste0(edges_a, collapse = \" \")} -> {colour_a}")

                # TODO debug message
                logger::log_debug("-- UPDATE E {paste0(edges_b, collapse = \" \")} -> {colour_b}")
            }
        }

        # Update communities on the edges of `a`
        if (!is.null(colour_a)) {
            edge_communities[edges_a] <- edge_communities[edges_a] %>% lapply(
                FUN = function(comm) {
                    union(comm, colour_a)
                }
            )
        }

        # Update communities on the edges of `b` (if needed)
        if (!is.null(edges_b) && !is.null(colour_b)) {
            edge_communities[edges_b] <- edge_communities[edges_b] %>% lapply(
                FUN = function(comm) {
                    union(comm, colour_b)
                }
            )
        }

        # Get number of examined edges incident on each node
        # plus add node_a and node_b so their n1 becomes even lower
        # (so we don't attempt to choose them again next time)

        ends_of_edges <- as.vector(igraph::ends(g, c(edges_a, edges_b), names = FALSE))
        # nodes_edges <- table(c(ends_of_edges, node_b))
        nodes_edges <- table(ends_of_edges)

        # Subtract number of classified incident edges from `n1`
        neighs_df[as.numeric(names(nodes_edges)), n1 := n1 - nodes_edges]
        neighs_df[c(node_a, node_b), n2 := n2 - 1]

        # update node communities
        node_communities <- get_node_communities(
            n, edge_communities, node_incident_edges
        )

        iter <- iter + 1

        if (debug_steps) {
            debug_log[[iter + 1]] <- list(
                iter = iter,
                node_a = node_a,
                node_b = node_b,
                colour_a = colour_a,
                colour_b = colour_b,
                edge_communities = lapply(edge_communities, sort),
                node_communities = lapply(node_communities, sort),
                neighs_df = data.table::copy(neighs_df)
            )
        }
    }

    # post-process: delete communities that are subsets of other communities
    if (post_process) {
        edge_communities <- edge_comm_postprocess(edge_communities)
        # update node communities, too
        node_communities <- get_node_communities(
            n, edge_communities, node_incident_edges
        )
    }

    return(list(
        edge_communities = lapply(edge_communities, sort),
        node_communities = lapply(node_communities, sort),
        iterations = iter,
        debug_steps = debug_log
    ))
}

#' Find colours that cannot be assigned to both `node_a` and `node_b` because
#' that would create a 2K_2.
#' @return A `vector` of exluded colours
#' @export
find_excluded_colors <- function(
    node_a, node_b, edges_a, edges_b,
    node_communities, edge_communities,
    neighs1, neighs2
) {
    # find colours used so far
    possible_colours <- unique.default(unlist(
        edge_communities[c(edges_a, edges_b)]
    ))

    # mask of excluded colours
    excluded_colours <- vector(
        mode = "logical",
        length = length(possible_colours)
    )

    # find 2-neighbours for each of the colours
    local_2neighs <- unlist(
        neighs2[c(node_a, node_b)],
        recursive = FALSE, use.names = FALSE
    )
    # remove `node_a` and `node_b` from 2-neighbours
    local_2neighs <- setdiff(local_2neighs, c(node_a, node_b))

    # check if each colour should be excluded
    for (i in seq_len(length(possible_colours))) {
        possible_col <- possible_colours[i]
        # get 2-neighbours with colour `possible_col`
        neigh2s_with_col <- local_2neighs[
            vapply(
                local_2neighs,
                FUN = function(node) {
                    any(possible_col %in% node_communities[[node]])
                },
                FUN.VALUE = logical(1)
            )
        ]

        logger::log_debug(
            "-- CHK COL {possible_col} ON ({node_a}, {node_b}) VS ({paste0(neigh2s_with_col, collapse = \", \")})"
        )

        # check each 2-neigbhbour with this colour against 2K2
        for (neigh2 in neigh2s_with_col) {
            if (
                is_2k2_neigh(neighs1[[node_a]], neighs1[[neigh2]])
                ||
                is_2k2_neigh(neighs1[[node_b]], neighs1[[neigh2]])
            ) {
                # exclude this colour
                logger::log_debug(
                    "-- EXCL {possible_col} ({node_a} or {node_b}, {neigh2})"
                )
                excluded_colours[i] <- TRUE
                break
            }
        }
    }

    return(possible_colours[excluded_colours])
}

find_same_color <- function(col1, col2, max_used_color, excluded_colours = c()) {
    # TODO find a good algorithm that doesn't create separate color classes
    # for multiple but nested pairs
    col1 <- unlist(col1)
    col2 <- unlist(col2)

    colors_a <- setdiff(unique.default(col1), excluded_colours)
    colors_b <- setdiff(unique.default(col2), excluded_colours)

    final_colours <- NULL

    # there are no colors in A
    if (length(colors_a) == 0 && length(colors_b) > 0) {
        final_colours <- colors_b
    } else if (length(colors_b) == 0 && length(colors_a) > 0) {
        final_colours <- colors_a
    }

    if (!is.null(final_colours)) {
        if (length(final_colours) > 0) {
            return(final_colours)
        } else {
            # there is no available colour
            return(c(max_used_color + 1))
        }
    }

    # get colors found in both groups
    matching_colors <- union(colors_a, colors_b)

    if (length(matching_colors) == 0) {
        return(c(max_used_color + 1))
    } else {
        final_colours <- matching_colors
        if (length(final_colours) > 0) {
            return(final_colours)
        } else {
            return(c(max_used_color + 1))
        }
    }
}

#' turn edge-communities into node communities:
#' a node belongs to all communities its incident edges belong to
#' IDs of incident edges.
#' @param n The number of vertices the graph has
#' @param edge_communities The edge_communities list of the algorithm
#' @param node_incident_edges List of which node is incident with which edges
#' (can be obtained using `igraph::incident_edges()`)
get_node_communities <- function(n, edge_communities, node_incident_edges) {
    lapply(
        X = seq_len(n),
        FUN = function(node) {
            # unwrap the list of communities into a vector and take its
            # unique values
            unique.default(
                unlist(
                    edge_communities[node_incident_edges[[node]]]
                )
            )
        }
    )
}

#' Perform edge community post-processing by removing communities that are
#' just subsets of other communities.
#' @param edge_communities The edge communities list of the algorithm
#' @param check_only Optionally a vector with colours to specifically check.
#' If `NULL`, all colours will be checked, otherwise only the specified ones.
#' @param remap Whether to remap communities after cleaning to the range `1:n`
#' @import RcppAlgos
edge_comm_postprocess <- function(
    edge_communities,
    check_only = NULL,
    remap = TRUE
    ) {
    colours <- unique.default(unlist(edge_communities))

    if (length(colours) < 2) {
        return(edge_communities)
    }

    if (!is.null(check_only)) {
        colours <- intersect(colours, check_only)
        if (length(colours) == 0) {
            # all colours were excluded, don't do anything
            return(edge_communities)
        }
    }
    colours <- sort(colours)
    excluded_colours <- logical(length(colours))

    colour_community_list <- vector(mode = "list", length = length(colours))
    col_pairs <- RcppAlgos::comboGeneral(v = colours, m = 2)

    # to which communities belong which edges
    for (i in seq_len(length(colours))) {
        colour_community_list[[i]] <- which(vapply(
            edge_communities,
            function(comm) i %in% comm,
            logical(1)
        ))
    }

    for (i in seq_len(nrow(col_pairs))) {
        if (all(
            colour_community_list[[col_pairs[i, 1]]] %in%
                colour_community_list[[col_pairs[i, 2]]]
        )) {
            # if the edges of [i, 1] are all among the edges of [i, 2],
            # exclude it
            excluded_colours[col_pairs[i, 1]] <- TRUE
            logger::log_debug("EXCL {col_pairs[i, 1]} %in% {col_pairs[i, 2]}")
        } else if (all(
            colour_community_list[[col_pairs[i, 2]]] %in%
                colour_community_list[[col_pairs[i, 1]]]
        )) {
            # cat("# EXCLUDING", col_pairs[i, 2], "\n")
            excluded_colours[col_pairs[i, 2]] <- TRUE
            logger::log_debug("EXCL {col_pairs[i, 2]} %in% {col_pairs[i, 1]}")
        }
    }

    excluded_colours <- colours[excluded_colours]

    # cat(paste0(
    #     "# POST-PROCESS: EXCLUDE ",
    #     paste0(excluded_colours, collapse = " "), "\n"
    # ))

    logger::log_debug("POST-PROCESS: EXCL {paste0(excluded_colours, collapse = ' ')}")

    # remove excluded colours
    new_communities <- lapply(
        edge_communities,
        function(comms) {
            setdiff(comms, excluded_colours)
        }
    )

    if (remap) {
        old_colours <- sort.default(unique.default(unlist(new_communities)))

        remap_df <- data.table(
            old = old_colours,
            new = seq_len(length(old_colours))
        )

        new_communities <- lapply(
            new_communities,
            function(comms) {
                remap_df$new[match(comms, remap_df$old)]
            }
        )

        logger::log_debug("## REMAP ({remap_df$old}) -> ({remap_df$new})")
    }

    return(new_communities)
}

#' Returns the list of communities with the edges or nodes they contain
#' @export
edges_per_community <- function(communities) {
    n_comm <- max(vapply(communities, max, numeric(1)))

    lapply(seq_len(n_comm), function(idx) {
        return(
            which(vapply(communities, function(comms) {
                return(any(idx %in% comms))
            }, FUN.VALUE = logical(1)))
        )
    })
}

#' Find communities that are stars
#' @import igraph
#' @return A binary mask of communities that are stars
#' @export
find_star_communities <- function(g, edge_communities) {
    n_comm <- max(vapply(edge_communities, max, numeric(1)))

    logger::log_debug("Checking {n_comm} communities...")

    edges_per_community <- lapply(seq_len(n_comm), function(idx) {
        return(
            which(vapply(edge_communities, function(edge_comms) {
                return(any(idx %in% edge_comms))
            }, FUN.VALUE = logical(1)))
        )
    })

    return(vapply(seq_len(n_comm), function(community) {
        checked_graph <- igraph::subgraph.edges(
            graph = g,
            eids = edges_per_community[[community]],
            delete.vertices = TRUE
        )
        degrees <- igraph::degree(checked_graph, mode = "out")
        vcount <- igraph::vcount(checked_graph)

        # return whether graph is a star
        return(
            igraph::count_components(checked_graph) == 1 &&
                sum(degrees) == 2 * (vcount - 1) &&
                any(degrees == (vcount - 1))
        )
    }, FUN.VALUE = logical(1)))
}

#' Compare how much communities overlap with each other
#' @export
community_overlap <- function(communities) {
    n <- length(communities)

    if (n < 2) {
        return(matrix(1))
    }

    community_pairs <- RcppAlgos::comboGeneral(
        v = seq_len(n), m = 2
    )

    result <- matrix(nrow = n, ncol = n)

    for (k in seq_len(nrow(community_pairs))) {
        i <- community_pairs[k, 1]
        j <- community_pairs[k, 2]

        logger::log_debug("Overlap checking communities: {i} -- {j}")

        comm_a <- communities[[i]]
        comm_b <- communities[[j]]
        overlap_percent <-
            length(intersect(comm_a, comm_b)) / length(union(comm_a, comm_b))

        result[i, j] <- result[j, i] <- overlap_percent
    }

    # set diagonal to 1
    diag(result) <- 1

    return(result)
}
