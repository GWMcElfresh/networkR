\name{AnalysisFunctions}
\alias{BuildBlacklist}
\alias{LearnBayesianNetwork}
\alias{CompareStructures}
\title{Bayesian network analysis helpers}
\usage{
BuildBlacklist(discretizationResult, configuration)

LearnBayesianNetwork(discretizationResult, configuration)

CompareStructures(firstStructure, secondStructure, firstLabel = "First",
  secondLabel = "Second")
}
\arguments{
  \item{discretizationResult}{A result from \code{DiscretizeMeasurements()}.}
  \item{configuration}{A networkR configuration object.}
  \item{firstStructure}{A result from \code{LearnBayesianNetwork()}.}
  \item{secondStructure}{A result from \code{LearnBayesianNetwork()}.}
  \item{firstLabel}{Label for the first structure.}
  \item{secondLabel}{Label for the second structure.}
}
\value{
\code{BuildBlacklist()} returns a list containing the blacklist and the active
exogenous variables.

\code{LearnBayesianNetwork()} returns a list with learned networks, bootstrap
strength estimates, and metadata.

\code{CompareStructures()} returns a list containing a comparison table.
}
\description{
Functions for constructing blacklist constraints, learning Bayesian networks,
and comparing structure strength across cohorts.
}