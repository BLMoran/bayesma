# Attach RoBMA Sensitivity Fits to a bayesma Object

Attaches pre-computed RoBMA sensitivity fits to a bayesma object. These
fits will be used by
[`sensitivity_plot()`](https://blmoran.github.io/bayesma/reference/sensitivity_plot.md)
when `incl_robma = TRUE`.

## Usage

``` r
attach_robma_sensitivity(model, robma_sensitivity)
```

## Arguments

- model:

  A `bayesma` object.

- robma_sensitivity:

  A `bayesma_robma_sensitivity` object created by
  [`robma_sensitivity()`](https://blmoran.github.io/bayesma/reference/robma_sensitivity.md).

## Value

The modified `bayesma` object with `$robma_sensitivity` attached.
