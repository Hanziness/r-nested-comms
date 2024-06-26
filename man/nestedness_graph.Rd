% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hclust.R
\name{nestedness_graph}
\alias{nestedness_graph}
\title{Transform a graph into a "nestedness graph", a full graph where the edge weights are the nestedness values}
\usage{
nestedness_graph(g)
}
\arguments{
\item{g}{The input \code{igraph} graph object}
}
\value{
A graph with edge weights equal to the nestedness values
}
\description{
Transform a graph into a "nestedness graph", a full graph where the edge weights are the nestedness values
}
