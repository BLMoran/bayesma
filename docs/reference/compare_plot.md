# Plot Model Comparison Results

Plot Model Comparison Results

## Usage

``` r
compare_plot(x, type = c("loso_crps", "calibration", "elpd", "pareto_k"), ...)
```

## Arguments

- x:

  A `bayesma_comparison` object.

- type:

  Character. Plot type: `"loso_crps"` (default): per-study CRPS from
  LOSO-CV. `"calibration"`: LOSO calibration curves. `"elpd"`:
  within-stage ELPD comparison. `"pareto_k"`: LOO Pareto k diagnostics.

- ...:

  Additional arguments (unused).

## Value

A `ggplot2` object.
