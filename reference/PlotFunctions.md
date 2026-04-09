# Plot Bayesian network results

Helpers for drawing averaged Bayesian network graphs with ggraph and
saving plots with package defaults.

## Usage

``` r
PlotBayesianNetwork(structureResult, configuration, title = NULL)

SavePlot(plotObject, filePath, configuration = NULL)
```

## Arguments

- structureResult:

  A result from
  [`LearnBayesianNetwork()`](https://gwmcelfresh.github.io/networkR/reference/AnalysisFunctions.md).

- configuration:

  A networkR configuration object.

- title:

  Optional plot title.

- plotObject:

  A ggplot object.

- filePath:

  Output file path.

## Value

`PlotBayesianNetwork()` returns a ggplot object.

`SavePlot()` returns the output file path, invisibly.
