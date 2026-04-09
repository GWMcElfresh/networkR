# Bayesian network analysis helpers

Functions for constructing blacklist constraints, learning Bayesian
networks, and comparing structure strength across cohorts.

## Usage

``` r
BuildBlacklist(discretizationResult, configuration)

LearnBayesianNetwork(discretizationResult, configuration)

CompareStructures(firstStructure, secondStructure, firstLabel = "First",
  secondLabel = "Second")
```

## Arguments

- discretizationResult:

  A result from
  [`DiscretizeMeasurements()`](https://gwmcelfresh.github.io/networkR/reference/DataFunctions.md).

- configuration:

  A networkR configuration object.

- firstStructure:

  A result from `LearnBayesianNetwork()`.

- secondStructure:

  A result from `LearnBayesianNetwork()`.

- firstLabel:

  Label for the first structure.

- secondLabel:

  Label for the second structure.

## Value

`BuildBlacklist()` returns a list containing the blacklist and the
active exogenous variables.

`LearnBayesianNetwork()` returns a list with learned networks, bootstrap
strength estimates, and metadata.

`CompareStructures()` returns a list containing a comparison table.
