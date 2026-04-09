# Configuration helpers

Helpers for constructing and validating the YAML-backed analysis
configuration.

## Usage

``` r
BuildConfiguration(overrides = list())

ParseConfiguration(filePath)

ValidateConfiguration(configuration)
```

## Arguments

- overrides:

  A named list of configuration overrides.

- filePath:

  Path to a YAML configuration file.

- configuration:

  A networkR configuration object.

## Value

For `BuildConfiguration()` and `ParseConfiguration()`, a validated
configuration object. For `ValidateConfiguration()`, the validated
configuration object, invisibly.
