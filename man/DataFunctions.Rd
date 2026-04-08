\name{DataFunctions}
\alias{ReadPanels}
\alias{AddDerivedVariables}
\alias{ImputeMeasurements}
\alias{DiscretizeMeasurements}
\title{Data preparation helpers}
\usage{
ReadPanels(configuration)

AddDerivedVariables(subjectTable, configuration)

ImputeMeasurements(subjectTable, configuration)

DiscretizeMeasurements(imputationResult, configuration)
}
\arguments{
  \item{configuration}{A networkR configuration object.}
  \item{subjectTable}{A subject-level table.}
  \item{imputationResult}{A result from \code{ImputeMeasurements()}.}
}
\value{
\code{ReadPanels()} returns a subject-level data frame.

\code{AddDerivedVariables()} returns the subject table with derived variables.

\code{ImputeMeasurements()} returns a list containing the imputed data and
analysis metadata.

\code{DiscretizeMeasurements()} returns a list containing continuous and
discrete analysis data.
}
\description{
Functions for reading panel files, deriving study variables, imputing missing
values, and discretizing measurements for Bayesian network structure learning.
}