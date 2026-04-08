\name{RunPipeline}
\alias{RunPipeline}
\title{Run the networkR analysis pipeline}
\usage{
RunPipeline(configuration)
}
\arguments{
  \item{configuration}{A networkR configuration object or a YAML file path.}
}
\value{
A pipeline result object containing the subject table, imputation result,
discretization result, full model result, and generated plot objects.
}
\description{
Runs the end-to-end networkR analysis workflow from configuration parsing
through structure learning and plot generation.
}