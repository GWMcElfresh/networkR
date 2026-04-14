\name{PlotFunctions}
\alias{PlotBayesianNetwork}
\alias{PlotEventTree}
\alias{PlotStageSplitHeatmap}
\alias{SavePlot}
\title{Plot Bayesian network and event tree results}
\usage{
PlotBayesianNetwork(structureResult, configuration, title = NULL)

PlotEventTree(structureResult, configuration, title = NULL,
  colorBy = "stage")

PlotStageSplitHeatmap(structureResult, configuration, title = NULL)

SavePlot(plotObject, filePath, configuration = NULL)
}
\arguments{
  \item{structureResult}{A result from \code{LearnBayesianNetwork()}.}
  \item{configuration}{A networkR configuration object.}
  \item{title}{Optional plot title.}
  \item{colorBy}{Node fill colour scheme: \code{"stage"} uses the staged event
    tree's stage colours; \code{"group"} uses measurement-group colours from
    the configuration.  Node borders always encode the measurement group.}
  \item{plotObject}{A ggplot object.}
  \item{filePath}{Output file path.}
}
\value{
\code{PlotBayesianNetwork()} returns a ggplot object.

\code{PlotEventTree()} returns a ggplot object, or \code{NULL} if the
structure result does not contain a staged event tree model.

\code{PlotStageSplitHeatmap()} returns a ggplot object, or \code{NULL} if no
stage probability summary can be constructed.

\code{SavePlot()} returns the output file path, invisibly.
}
\description{
Helpers for drawing averaged Bayesian network graphs and staged event trees
with ggraph, visualizing stage-to-split probabilities for stagedtrees models,
and for saving plots with package defaults.
}