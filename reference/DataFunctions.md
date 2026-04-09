# Data preparation helpers

Functions for reading panel files, deriving study variables, imputing
missing values, and discretizing measurements for Bayesian network
structure learning.

## Usage

``` r
ReadPanels(configuration)

AddDerivedVariables(subjectTable, configuration)

ImputeMeasurements(subjectTable, configuration)

DiscretizeMeasurements(imputationResult, configuration)
```

## Arguments

- configuration:

  A networkR configuration object.

- subjectTable:

  A subject-level table.

- imputationResult:

  A result from `ImputeMeasurements()`.

## Value

`ReadPanels()` returns a subject-level data frame.

`AddDerivedVariables()` returns the subject table with derived
variables.

`ImputeMeasurements()` returns a list containing the imputed data and
analysis metadata.

`DiscretizeMeasurements()` returns a list containing continuous and
discrete analysis data.
