# Add or remove probability columns from sensitivity plot

Always shows Pr(Benefit) and Pr(Harm). Only shows the null range
probability columns (Pr(Benefit\>δ) and Pr(Harm\>δ)) when a null range
has been specified — either in the original
[`sensitivity_plot()`](https://blmoran.github.io/bayesma/reference/sensitivity_plot.md)
call or via
[`sens_add_null()`](https://blmoran.github.io/bayesma/reference/sens_add_null.md).

## Usage

``` r
sens_add_probs(x, show = TRUE)
```

## Arguments

- x:

  A `bayesma_sensitivity_plot` object.

- show:

  Logical. If TRUE (default), show probability columns. If FALSE, hide
  them.

## Value

Modified `bayesma_sensitivity_plot`.
