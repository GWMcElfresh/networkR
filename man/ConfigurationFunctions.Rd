\name{ConfigurationFunctions}
\alias{BuildConfiguration}
\alias{ParseConfiguration}
\alias{ValidateConfiguration}
\title{Configuration helpers}
\usage{
BuildConfiguration(overrides = list())

ParseConfiguration(filePath)

ValidateConfiguration(configuration)
}
\arguments{
  \item{overrides}{A named list of configuration overrides.}
  \item{filePath}{Path to a YAML configuration file.}
  \item{configuration}{A networkR configuration object.}
}
\value{
For \code{BuildConfiguration()} and \code{ParseConfiguration()}, a validated
configuration object. For \code{ValidateConfiguration()}, the validated
configuration object, invisibly.
}
\description{
Helpers for constructing and validating the YAML-backed analysis configuration.
}