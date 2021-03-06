% Generated by roxygen2 (4.1.0.9000): do not edit by hand
% Please edit documentation in R/analysis_functions.R
\name{transformation.matrix}
\alias{transformation.matrix}
\title{Bi-partite network analysis tools}
\usage{
transformation.matrix(network.1, network.2, by.genes = F, standardize = T,
  remove.diagonal = T)
}
\arguments{
\item{net1}{starting network, a genes by transcription factors data.frame with scores for confidence in the existence of edges between}

\item{net2}{final network, a genes by transcription factors data.frame with scores for confidence in the existence of edges between}
}
\description{
This function analyzes a bi-partite network, such as a Transcription factor to gene network derived from the PANDA algorithm.
}
\examples{
data(yeast.panda)
t.matrix <- transformation.matrix(yeast.panda$cell.cycle, yeast.panda$stress.response)
hcl.heatmap.plot(t.matrix, method="pearson")
}
\keyword{keywords}

