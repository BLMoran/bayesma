# Default RoBMA prior set

Returns the default list of priors for one of the three RoBMA dimensions
(effect, heterogeneity, bias). Used internally by
[`robma()`](https://blmoran.github.io/bayesma/reference/robma.md) when
the user does not pass explicit priors.

## Usage

``` r
robma_default_priors(
  dimension = c("effect", "heterogeneity", "bias"),
  null = FALSE,
  rescale = 1
)
```

## Arguments

- dimension:

  One of `"effect"`, `"heterogeneity"`, `"bias"`.

- null:

  Logical. If `TRUE`, return the null-model priors (point at zero for
  effect / heterogeneity, no-bias for bias).

- rescale:

  Numeric. Multiplier on default scales.

## Value

A list of `bayesma_prior` or `robma_bias_prior` objects.
