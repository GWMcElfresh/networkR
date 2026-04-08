# networkR

`networkR` is an R package scaffold for configurable Bayesian network analysis in
factorial and repeated-measures studies. It extracts the main ideas from the
current [Malaria_U01_BayesianNetwork.Rmd](Malaria_U01_BayesianNetwork.Rmd) and
turns them into reusable package functions.

## Current scope

- YAML or programmatic configuration
- Panel ingestion from `data_sheets/`
- Derived design variables
- MICE predictive mean matching imputation
- Hartemink discretization with fallback handling
- Tier-aware blacklist construction
- Bootstrap Bayesian network fitting with `bnlearn`
- `ggraph(layout = "fr")` network plotting with a `strength.plot()` fallback path
- CI scaffolding, test scaffolding, and pkgdown scaffolding

## Quick start

```r
library(networkR)

configuration <- ParseConfiguration("inst/example_configuration.yml")
analysisResult <- RunPipeline(configuration)

analysisResult$plotObjects$fullModel
```

## Status

This repository now contains the package foundation. The remaining work is to
fill out more of the specialized analyses from the R Markdown document,
including joint paired-tissue modeling, mediation workflows, and correlation
reporting.
