% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edgefinder.R
\name{shortestPathDistance}
\alias{shortestPathDistance}
\title{Return a Matrix with the shortest path distance between nodes (check up to numSteps.)}
\usage{
shortestPathDistance(AdjMat, numSteps = 0)
}
\arguments{
\item{AdjMat}{An adjacency Matrix (0/1).}

\item{numSteps}{The maximum number of edges between pairs of nodes. If numSteps=0, returns the input matrix. numSteps=1 adds neighbors of direct neighbors, etc.}
}
\value{
A Matrix containing the shortset paths between nodes i and j
}
\description{
return the adjacency matrix of expMat connecting neighbors up to numSteps away.
}
\examples{
\donttest{
   data(SIM)
   Sres <- edgefinder(SIM, ttl = "hub network")
   AdjMat1 <- shortestPathDistance(Sres$AdjMat, numSteps=50)
   max(AdjMat1)
   Matrix::image(AdjMat1)
}
}
