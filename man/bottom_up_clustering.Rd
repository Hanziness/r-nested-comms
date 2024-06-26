% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hclust.R
\name{bottom_up_clustering}
\alias{bottom_up_clustering}
\title{Bottom-up hierarchical nestedness clustering using R's built-in \code{hclust} function}
\usage{
bottom_up_clustering(g, method = "complete", level_per_node = TRUE)
}
\arguments{
\item{g}{The input base graph (an \code{igraph} graph object)}

\item{method}{The \code{method} parameter passed to \code{hclust} (e.g., \code{single}, \code{complete}, \code{average})}

\item{level_per_node}{Boolean indicating whether the resulting list should have a level per each node, instead of a level per each unique merge distance. The latter results in much fewer levels.}
}
\value{
A list of each cluster level
}
\description{
Bottom-up hierarchical nestedness clustering using R's built-in \code{hclust} function
}
