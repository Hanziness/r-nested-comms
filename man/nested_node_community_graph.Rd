% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nested_node_comm.R
\name{nested_node_community_graph}
\alias{nested_node_community_graph}
\title{Build the community graph using the resulting community structure of \code{nested_node_comms}.}
\usage{
nested_node_community_graph(results, sameneigh_original_g = NULL)
}
\arguments{
\item{results}{The community structure obtained using \code{nested_node_comms}}

\item{sameneigh_original_g}{(Optional) The original graph \code{nested_node_comms} was performed on. If given, bidirectional edges will be added between vertices with the same neighborhood.}
}
\value{
An \code{igraph} graph of the community graph
}
\description{
Build the community graph using the resulting community structure of \code{nested_node_comms}.
}
