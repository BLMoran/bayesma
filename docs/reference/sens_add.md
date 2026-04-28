# Post-render modifications to a sensitivity plot

Helpers for adjusting an existing `bayesma_sensitivity_plot` object
without re-running the full pipeline.

## Usage

``` r
sens_add_probs(x, show = TRUE)

sens_add_null(x, range = NULL, color_null_range = "#77bb41")

sens_add_titles(x, title = NULL, subtitle = NULL, align = NULL)

sens_add_x_lim(x, xlim, x_breaks = NULL)

sens_add_plot_width(x, plot_width)
```

## Arguments

- x:

  A `bayesma_sensitivity_plot` object.

- show:

  Logical. Whether to display posterior probability columns.

- range:

  Numeric vector of length 2 giving the null/ROPE range, or `NULL` to
  use the default for the model's estimand.

- color_null_range:

  Colour used to shade the null range.

- title, subtitle:

  Character. Plot titles.

- align:

  Title alignment: `"left"`, `"center"`, or `"right"`.

- xlim:

  Numeric vector of length 2. New x-axis limits.

- x_breaks:

  Optional numeric vector of x-axis break locations.

- plot_width:

  Positive numeric. Relative width of the density panel.

## Value

The updated `bayesma_sensitivity_plot` object.
