% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{is_2k2_neigh_value}
\alias{is_2k2_neigh_value}
\title{Checks if there is a $2K_2$ between two vertices in graph \code{g}
by calculating $|N(v1) cap N(v2)| / min(|N(v1)|, |N(v2)|)$.
This version uses only the neighborhoods of said vertices, optionally
excluding \code{v1} and \code{v2} from each other's neighborhoods if \code{v1} or
\code{v2} are non-negative.}
\usage{
is_2k2_neigh_value(n1, n2, v1 = -1L, v2 = -1L)
}
\arguments{
\item{n1}{The neighborhood of \code{v1}}

\item{n2}{The neighborhood of \code{v2}}

\item{v1}{(Optional) The first vertex of the pair}

\item{v2}{(Optional) The second vertex of the pair}

\item{treshold}{A value (between 0 and 1) that is still accepted as nested.}
}
\value{
A boolean indicating whether there is a 2K_2 between \code{v1} and \code{v2}
}
\description{
Checks if there is a $2K_2$ between two vertices in graph \code{g}
by calculating $|N(v1) cap N(v2)| / min(|N(v1)|, |N(v2)|)$.
This version uses only the neighborhoods of said vertices, optionally
excluding \code{v1} and \code{v2} from each other's neighborhoods if \code{v1} or
\code{v2} are non-negative.
}
