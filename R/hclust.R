#' Create a distance matrix based on nestedness: distance(i, j) = 1 - nest(i, j)
#' @param g The input `igraph` graph object
#' @return A matrix of size `igraph::vcount(g)` x `igraph::vcount(g)`, where the (i, j)-th entry contains the distance between the i-th and j-th vertices
#' @export
dist_nested <- function(g) {
    neighs <- igraph::neighborhood(g, order = 1, mode = "all", mindist = 1)
    n <- igraph::vcount(g)
    # distance = 1 - nestedness
    distmatrix <- matrix(data = 1, nrow = n, ncol = n)
    for (i in as.numeric(igraph::V(g))) {
        for (j in as.numeric(neighs[[i]])) {
            for (k in as.numeric(neighs[[j]])) {
                if (k < i) {
                    next
                }

                nestedness <- nested.comms::nestedness_value(neighs[[i]], neighs[[k]], i, k)
                distmatrix[k, i] <- 1 - nestedness
            }
        }
    }

    return(distmatrix)
}

#' Transform a graph into a "nestedness graph", a full graph where the edge weights are the nestedness values
#' @param g The input `igraph` graph object
#' @return A graph with edge weights equal to the nestedness values
#' @export
nestedness_graph <- function(g) {
    dist_m <- 1 - dist_nested(g)
    return(igraph::graph_from_adjacency_matrix(dist_m, weighted = TRUE, diag = FALSE, mode = "undirected"))
}

#' Bottom-up hierarchical nestedness clustering using R's built-in `hclust` function
#' @param g The input base graph (an `igraph` graph object)
#' @param method The `method` parameter passed to `hclust` (e.g., `single`, `complete`, `average`)
#' @param level_per_node Boolean indicating whether the resulting list should have a level per each node, instead of a level per each unique merge distance. The latter results in much fewer levels.
#' @return A list of each cluster level
#' @export
bottom_up_clustering <- function(g, method = "complete", level_per_node = TRUE) {
    res <- hclust(as.dist(dist_nested(g)), method = method)
    if (level_per_node) {
        return(lapply(rev(seq_len(igraph::vcount(g))), function(i) cutree(res, k = i)))
    } else {
        return(hclust_to_membershiplist(res))
    }
}

#' Top-down hierarchical nestedness clustering using the Girvan-Newman algorithm
#' @param nested_g The input nestedness `igraph` graph object obtained using `nestedness_graph()`
#' @param method Either `full` (closer to the original Girvan-Newman implementation), or `nested` (faster, does not calculate edge betweenness, only uses nestedness values from the edges)
#' @return A list of vectors of length `igraph::vcount(g)`, where the i-th vector contains the membership of the i-th vertex. Each vector corresponds to a level (step) of hierarchical clustering.
#' @export
girvan_newman_nestedness <- function(nested_g, method = c("nested", "full")) {
    method <- match.arg(method)
    result_list <- list(igraph::components(nested_g)$membership)
    if (!("weight" %in% names(igraph::edge.attributes(nested_g)))) {
        igraph::E(nested_g)$weight <- 1
    }

    while (igraph::ecount(nested_g) > 0) {
        if (method == "nested") {
            nested_g <- girvan_newman_step_nested(nested_g)
        } else {
            nested_g <- girvan_newman_step_full(nested_g)
        }
        result_list <- c(result_list, list(igraph::components(nested_g)$membership))
    }

    return(result_list)
}

#' Performs an iteration of the top-down hierarchical nested clustering algorithm
#' @param g The input nestedness graph object
#' @param method The same as `method` in `girvan_newman_nestedness`
girvan_newman_step_nested <- function(g) {
    original_components <- igraph::count_components(g, mode = "weak")
    new_components <- original_components

    if (original_components == igraph::vcount(g)) {
        logger::log_error("Graph is already disconnected")
        return(g)
    }

    igraph::V(g)$original_membership <- igraph::components(g)$membership

    while (new_components <= original_components) {
        # B) remove the edges with the lowest nestedness.
        # If there are multiple, remove those with maximum edge betweenness
        e <- igraph::E(g)
        e_min_nestedness <- which(e$weight == min(e$weight))
        delete_edges <- e[e_min_nestedness]

        if (length(e_min_nestedness) > 1) {
            e_betweenness <- igraph::edge_betweenness(g, e = delete_edges)
            delete_edges <- delete_edges[e_betweenness == max(e_betweenness)]
        }
        logger::log_debug(sprintf("Deleted %d edge(s) (at %d components)", length(delete_edges), original_components))

        g <- igraph::delete_edges(g, delete_edges)
        new_components <- igraph::count_components(g, mode = "weak")
    }

    return(g)
}

#' Performs an iteration of the top-down hierarchical nested clustering algorithm using the "full" version
#' @param g The input nestedness graph object
girvan_newman_step_full <- function(g) {
    original_components <- igraph::count_components(g, mode = "weak")
    new_components <- original_components

    if (original_components == igraph::vcount(g)) {
        logger::log_error("Graph is already disconnected")
        return(g)
    }

    igraph::V(g)$original_membership <- igraph::components(g)$membership

    while (new_components <= original_components) {
        # remove the edge with the highest edge betweenness, scaled with nestedness (= edge weight)
        e_values <- igraph::edge.betweenness(g) / igraph::E(g)$weight
        delete_edges <- igraph::E(g)[which.max(e_values)]
        logger::log_debug(sprintf("Deleted %d edge(s) (at %d components)", length(delete_edges), original_components))

        g <- igraph::delete_edges(g, delete_edges)
        new_components <- igraph::count_components(g, mode = "weak")
    }
    logger::log_info(sprintf("Graph now has %d components", new_components))

    return(g)
}

#' @export
plot_cluster_steps <- function(g, results_list, steps = NULL, include_nestedness = FALSE, file = NA, layout_fn = NULL, layout = NULL, ...) {
    if (!is.na(file)) {
        cairo_pdf(file, width = 7, height = 7, onefile = TRUE, bg = "transparent")
        # old_par <- par(bg = NA)
    }

    if (is.null(steps)) {
        steps <- seq_along(results_list)
    } else {
        steps <- intersect(steps, seq_along(results_list))
    }

    if (include_nestedness) {
        # par(mfrow = c(1, 2))
        nested_plot_base <- ggplot2::ggplot(reshape2::melt(partition_dataframe(g, results_list), id.vars = c("step"), measure.vars = c("mean_weighted", "mean_simple"))) + ggplot2::geom_line(ggplot2::aes(x = step, y = value, color = variable)) + ggplot2::theme_minimal() + ggplot2::labs(color = "", y = "nest") + ggplot2::scale_color_discrete(labels = c("Weighted", "Simple")) # nolint: object_usage_linter.
    }

    if (is.null(layout)) {
        if (is.null(layout_fn)) {
            layout <- if (igraph::is_bipartite(g)) igraph::layout_as_bipartite(g) else igraph::layout_with_fr(g)
        } else {
            layout <- layout_fn(g)
        }
    }

    stabilized_coloring <- results_list
    colors <- sample(rainbow(length(results_list[[1]])))

    for (i in seq.int(from = 2, to = length(results_list))) {
        stabilized_coloring[[i]] <- relabel_membership(stabilized_coloring[[i - 1]], stabilized_coloring[[i]])
    }

    for (i in steps) {
        logger::log_debug("Plotting step {i}")
        plot(
            g,
            vertex.color = colors[stabilized_coloring[[i]]],
            layout = layout,
            ...
        )

        # if (include_nestedness == TRUE) {
        #     print(nested_plot_base + ggplot2::geom_vline(ggplot2::aes(xintercept = i), color = "red"))
        # }
    }

    if (!is.na(file)) {
        dev.off()
        # par(old_par)
    }

    if (include_nestedness == TRUE) {
        if (!is.na(file)) {
            cairo_pdf(paste0(tools::file_path_sans_ext(file), "_nestedness.pdf"), width = 4, height = 3, onefile = TRUE)
        }

        for (i in seq_along(results_list)) {
            print(nested_plot_base + ggplot2::geom_vline(ggplot2::aes(xintercept = i), color = "red"))
        }

        if (!is.na(file)) {
            dev.off()
        }
    }
}

relabel_membership <- function(old_m, new_m) {
    original_labels <- sort(unique.default(old_m))
    available_labels <- sort(unique.default(new_m))
    remap <- data.frame(from = c(), to = c())

    for (i in original_labels) {
        indices <- which(old_m == i)
        matching_label <- which.max(sapply(available_labels, function(j) {
            sum(new_m[indices] == j) / length(indices)
        }))

        remap <- rbind(remap, data.frame(
            from = matching_label,
            to = i
        ))
    }

    # find left-out label
    remap <- rbind(remap, data.frame(
        from = setdiff(available_labels, remap$from),
        to = setdiff(available_labels, remap$to)
    ))

    return(remap$to[match(new_m, remap$from)])
}

partition_nestedness <- function(g, partition, nested_m = NULL) {
    if (is.null(nested_m)) {
        nested_m <- 1 - dist_nested(g)
    }

    # calculate nestedness for each partition
    nestedness_values <- lapply(seq_len(length(unique.default(partition))), function(i) {
        nestedness <- as.vector(as.dist(nested_m[partition == i, partition == i]))
        if (length(nestedness) == 0) {
            return(c(NA))
        }
        return(nestedness)
    })

    return(nestedness_values)
}

partition_nestedness_all <- function(g, results_list) {
    nested_m <- 1 - dist_nested(g)
    return(lapply(results_list, function(res) partition_nestedness(g, res, nested_m)))
}

#' Calculates the mean nestedness for each partition of each level in a result list
#' @param g The original graph
#' @param results_list The result list obtained from `girvan_newman_nestedness`
#' @return A list of vectors (each of length `k` - the number of partitions) of length `igraph::vcount(g)`, where the i-th vector contains the mean nestedness of the i-th level's partitions.
#' @export
mean_partition_nestedness <- function(g, results_list) {
    return(lapply(partition_nestedness_all(g, results_list), function(res) {
        return(sapply(res, function(x) mean(x, na.rm = TRUE)))
    }))
}

#' Compute performance properties for each clustering level in a result list
#' @param g The original graph
#' @param results_list The result list returned by `girvan_newman_nestedness`
#' @return A data frame with properties of the clustering levels (number of clusters, mean nestedness, etc.)
#' @export
partition_dataframe <- function(g, results_list) {
    n <- length(results_list[[1]])
    weights <- lapply(results_list, function(x) table(x))
    means <- mean_partition_nestedness(g, results_list)
    na_values <- lapply(means, function(x) which(is.na(x)))
    means_zero <- lapply(means, function(x) replace(x, is.na(x), 0))
    means_one <- lapply(means, function(x) replace(x, is.na(x), 1))
    means_omitna <- lapply(means, function(x) {
        res <- as.vector(na.omit(x))
        if (length(res) == 0) {
            c(1)
        } else {
            res
        }
    })
    weights_omitna <- lapply(seq_along(weights), function(i) {
        if (length(na_values[[i]]) == 0) {
            return(weights[[i]])
        } else if (length(na_values[[i]]) == length(weights[[i]])) {
            return(c(1))
        }

        weights[[i]][-na_values[[i]]]
    })

    return(data.frame(
        step = seq_along(results_list),
        clusters = vapply(results_list, function(x) length(unique.default(x)), FUN.VALUE = numeric(1)),
        progress = vapply(results_list, function(x) length(unique.default(x)) / n, FUN.VALUE = numeric(1)),
        mean_weighted = vapply(seq_along(results_list), function(i) weighted.mean(means_zero[[i]], w = weights[[i]]), FUN.VALUE = numeric(1)),
        mean_simple = vapply(means_zero, mean, FUN.VALUE = numeric(1)),
        mean_weighted_one = vapply(seq_along(results_list), function(i) weighted.mean(means_one[[i]], w = weights[[i]]), FUN.VALUE = numeric(1)),
        mean_simple_one = vapply(means_one, mean, FUN.VALUE = numeric(1)),
        mean_weighted_omitna = vapply(seq_along(results_list), function(i) weighted.mean(means_omitna[[i]], w = weights_omitna[[i]]), FUN.VALUE = numeric(1)),
        mean_simple_omitna = vapply(means_omitna, mean, FUN.VALUE = numeric(1)),
        ratio_fullynested_weighted = vapply(seq_along(results_list), function(i) weighted.mean(means_one[[i]] == 1, w = weights[[i]]), FUN.VALUE = numeric(1)),
        ratio_fullynested = vapply(means_one, function(x) sum(x == 1) / length(x), FUN.VALUE = numeric(1)),
        ratio_fullynested_weighted_omitna = vapply(seq_along(results_list), function(i) weighted.mean(means_omitna[[i]] == 1, w = weights_omitna[[i]]), FUN.VALUE = numeric(1)),
        ratio_fullynested_omitna = vapply(means_omitna, function(x) sum(x == 1) / length(x), FUN.VALUE = numeric(1)),
        sd = vapply(means_zero, function(x) sd(x, na.rm = TRUE), FUN.VALUE = numeric(1))
    ))
}

partition_plot <- function(g, results_list) {
    return(ggplot2::ggplot(partition_dataframe(g, results_list)))
}

#' Converts a `hclust` results object into a list of cluster levels (one level per unique merge distance)
hclust_to_membershiplist <- function(hclust_res) {
    # merge levels (distances)
    levels <- sort(unique.default(cophenetic(hclust_res)))

    # get membership for each level
    membership_list <- lapply(levels, function(level) {
        cutree(hclust_res, h = level)
    })

    return(membership_list)
}

#' (Experimental) Perform statistical pre-filtering on a graph, removing edges that are unlikely to be nested
#' @param g The input `igraph` graph object
#' @export
stat_prefilter <- function(g, alpha = 0.05) {
    if (alpha <= 0.00) {
        return(g)
    }

    edges <- igraph::E(g)
    neighs <- igraph::neighborhood(g, order = 1, mode = "all", mindist = 1)
    edges_to_delete <- c()
    n <- igraph::vcount(g)


    for (e in edges) {
        ends <- igraph::ends(g, e, names = FALSE)
        i <- ends[1, 1]
        j <- ends[1, 2]

        k_i <- length(neighs[[i]])
        k_j <- length(neighs[[j]])

        if (min(k_i, k_j) == 0) {
            p <- 0
        } else {
            p <- 1 - phyper(
                q = min(k_i, k_j) - 1,
                m = k_i,
                n = n - k_i,
                k = k_j,
                lower.tail = TRUE
            )
        }

        if (p > alpha) {
            edges_to_delete <- c(edges_to_delete, e)
        }
    }

    logger::log_debug(sprintf("Deleting %d edges (%.1f%%)\n", length(edges_to_delete), length(edges_to_delete) / igraph::ecount(g) * 100))
    return(igraph::delete_edges(g, edges_to_delete))
}

#' @export
experiment_comparefilters <- function(g, alpha = c(0.01, 0.05, 0.10)) {
    nested_g <- nestedness_graph(g)
    dt_list <- list(cbind(
        data.frame(method = "original"),
        partition_dataframe(g, girvan_newman_nestedness(nested_g))
    ))

    for (a in alpha) {
        dt_list <- c(dt_list, list(cbind(
            data.frame(method = sprintf("stat_prefilter_%.2f", a)),
            partition_dataframe(g, girvan_newman_nestedness(stat_prefilter(nested_g, alpha = a)))
        )))
    }

    do.call(rbind, dt_list)
}