\name{PlotFunctions}
\alias{PlotBayesianNetwork}
\alias{SavePlot}
\title{Plot Bayesian network results}
\usage{
PlotBayesianNetwork(structureResult, configuration, title = NULL)

SavePlot(plotObject, filePath, configuration = NULL)
}
\arguments{
  \item{structureResult}{A result from \code{LearnBayesianNetwork()}.}
  \item{configuration}{A networkR configuration object.}
  \item{title}{Optional plot title.}
  \item{plotObject}{A ggplot object.}
  \item{filePath}{Output file path.}
}
\value{
\code{PlotBayesianNetwork()} returns a ggplot object.

\code{SavePlot()} returns the output file path, invisibly.
}
\description{
Helpers for drawing averaged Bayesian network graphs with ggraph and saving
plots with package defaults.
}