\name{AnalysisFunctions}
\alias{BuildBlacklist}
\alias{LearnBayesianNetwork}
\alias{CompareBackendBootstrapStrengths}
\alias{CompareStructures}
\alias{CompareStages}
\title{Bayesian network analysis helpers}
\usage{
BuildBlacklist(discretizationResult, configuration)

LearnBayesianNetwork(discretizationResult, configuration)

CompareBackendBootstrapStrengths(firstStructure, secondStructure,
  firstLabel = "First", secondLabel = "Second", tolerance = 0.05)

CompareStructures(firstStructure, secondStructure, firstLabel = "First",
  secondLabel = "Second")

CompareStages(firstStructure, secondStructure, firstLabel = "First",
  secondLabel = "Second")
}
\arguments{
  \item{discretizationResult}{A result from \code{DiscretizeMeasurements()}.}
  \item{configuration}{A networkR configuration object.}
  \item{firstStructure}{A result from \code{LearnBayesianNetwork()}.}
  \item{secondStructure}{A result from \code{LearnBayesianNetwork()}.}
  \item{firstLabel}{Label for the first structure.}
  \item{secondLabel}{Label for the second structure.}
  \item{tolerance}{Absolute difference threshold used to count exceeded directed-edge bootstrap differences.}
}
\value{
\code{BuildBlacklist()} returns a list containing the blacklist and the active
exogenous variables.

\code{LearnBayesianNetwork()} returns a list with learned networks, bootstrap
strength estimates, and metadata.

\code{CompareBackendBootstrapStrengths()} returns a list containing a directed-edge
comparison table and summary metrics for bootstrap parity.

\code{CompareStructures()} returns a list containing a comparison table.

\code{CompareStages()} returns a list containing a comparison table, Hamming
distance, and context-specific independence distance (CID) for staged event
trees.
}
\description{
Functions for constructing blacklist constraints, learning Bayesian networks,
and comparing structure strength across cohorts. \code{CompareBackendBootstrapStrengths()}
quantifies backend parity on directed bootstrap strengths. \code{CompareStages()}
is the stagedtrees-specific variant of \code{CompareStructures()} and uses
Hamming distance and CID metrics.
}