# Plot Model Comparison Results

Plot Model Comparison Results

## Usage

``` r
# S3 method for class 'bayesma_comparison'
plot(x, type = c("elpd", "pareto_k"), ...)
```

## Arguments

- x:

  A `bayesma_comparison` object.

- type:

  Character. Type of plot: `"elpd"` (default) for ELPD comparison, or
  `"pareto_k"` for LOO diagnostics.

- ...:

  Additional arguments (unused).

## Value

A `ggplot2` object.
