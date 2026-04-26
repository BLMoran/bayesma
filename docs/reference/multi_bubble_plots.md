# Multi-panel Bubble Plots

Creates bubble plots for all continuous moderators in a meta-regression.

## Usage

``` r
multi_bubble_plots(object, ncol = 2, ...)
```

## Arguments

- object:

  A `bayesma_reg` object.

- ncol:

  Integer. Number of columns in the plot grid (default: 2).

- ...:

  Additional arguments passed to
  [`bubble_plot()`](https://blmoran.github.io/bayesma/reference/bubble_plot.md).

## Value

A combined ggplot object (using patchwork if available).
