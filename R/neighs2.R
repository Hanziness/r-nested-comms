#' Return vertices that are reachable from `nodes` in two steps.
#' As opposed to `igraph::neighborhood(..., order = 2)` this functiondoes not return immediate neighbors if they are not reachable in two steps.
#' @param g An `igraph` graph object
#' @param nodes (Optional) The set of nodes for which the list of 2-neighbors will be returned.
#' @export
neighborhood_2 <- function(g, nodes = igraph::V(g)) {
    if (length(nodes) < igraph::vcount(g)) {
        g <- igraph::subgraph(graph = g, vids = nodes)
    }

    neighs1 <- lapply(igraph::ego(
        graph = g, nodes = igraph::V(g), mindist = 1, order = 1
    ), as.numeric)

    return(
        lapply(seq_len(igraph::vcount(g)), function(i) {
            neighs_next <- unlist(neighs1[neighs1[[i]]])
            sort(unique.default(neighs_next[neighs_next != i]))
        })
    )
}