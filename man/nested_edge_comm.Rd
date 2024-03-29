% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nested_edge_comm.R
\name{nested_edge_comm}
\alias{nested_edge_comm}
\title{Edge-based nestedness community detection}
\usage{
nested_edge_comm(
  g,
  mode = c("fast", "precise", "strict"),
  n1_threshold = 1L,
  post_process = TRUE,
  debug_steps = FALSE
)
}
\arguments{
\item{g}{The input \code{igraph} graph object}

\item{mode}{Controsl how many checks the algorithm will perform before
returning the communities. Might result in better edge community detection
but will also require more iterations to complete. Might find some
false communities that are subsets of others, these can be removed using
\code{post_process = TRUE}. The three accepted values are \code{fast} (n1 > 0),
\code{precise} (n1 >= 0) and \code{strict} (n2 > 0), the last one comparing all
accepted node pairings.}

\item{n1_threshold}{The threshold parameter to use in the iteration
condition. \code{1} is the "fast" mode, "0" is the "precise" mode and \code{-Inf} is
the "strict" mode.}

\item{post_process}{Whether to remove communities that are subsets
of other communities. Useful when \code{mode = "precise"} or \code{mode = "strict"}.}
}
\description{
Edge-based nestedness community detection
}
