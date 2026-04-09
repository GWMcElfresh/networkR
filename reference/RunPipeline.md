# Run the networkR analysis pipeline

Runs the end-to-end networkR analysis workflow from configuration
parsing through structure learning and plot generation.

## Usage

``` r
RunPipeline(configuration)
```

## Arguments

- configuration:

  A networkR configuration object or a YAML file path.

## Value

A pipeline result object containing the subject table, imputation
result, discretization result, full model result, and generated plot
objects.
