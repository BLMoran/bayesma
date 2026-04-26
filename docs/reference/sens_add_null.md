# Add or update a null range (ROPE) on the sensitivity plot

Adds the null range shading to the density plot and recalculates the
null range probability columns (Pr(Benefit\>δ) and Pr(Harm\>δ)). If
[`sens_add_probs()`](https://blmoran.github.io/bayesma/reference/sens_add_probs.md)
has been called, the δ columns will automatically appear in the table.

## Usage

``` r
sens_add_null(x, range = NULL, color_null_range = "#77bb41")
```

## Arguments

- x:

  A `bayesma_sensitivity_plot` object.

- range:

  Numeric vector of length 2 giving the null range bounds (on the effect
  scale, e.g. `c(0.9, 1.1)` for OR). If `NULL` (default), a
  measure-appropriate default is used: OR/RR/HR/IRR = `c(0.9, 1.1)`, SMD
  = `c(-0.1, 0.1)`. A single value is interpreted as symmetric around
  the null (e.g. `0.1` becomes `c(null - 0.1, null + 0.1)`).

- color_null_range:

  Fill colour for the null range band. Default `"#77bb41"`.

## Value

Modified `bayesma_sensitivity_plot`.
