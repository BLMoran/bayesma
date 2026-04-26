# Plot Method for bayesma_coef_evidence Objects

Creates a forest-style plot for meta-regression coefficients using
posterior density slabs, consistent with the `bayes_forest` visual
style.

## Usage

``` r
metareg_mod_plot(
  x,
  model = NULL,
  include_intercept = FALSE,
  show_null_range = TRUE,
  null_value = 0,
  output = c("density", "pointinterval"),
  split_color_by_null = FALSE,
  color_posterior = "dodgerblue",
  color_posterior_outline = "blue",
  color_favours_positive = "dodgerblue",
  color_favours_negative = "firebrick",
  color_null_range = "grey50",
  color_pointinterval = "blue",
  xlim = NULL,
  x_breaks = NULL,
  xlab = "Coefficient Estimate",
  title = "Meta-Regression Coefficients",
  subtitle = NULL,
  add_table = TRUE,
  table_width = 0.4,
  font = NULL,
  ...
)
```

## Arguments

- x:

  A `bayesma_coef_evidence` object (created by
  [`coefficient_evidence()`](https://blmoran.github.io/bayesma/reference/coefficient_evidence.md)).

- model:

  The `bayesma_metareg` model. Required to extract posterior draws for
  density slabs.

- include_intercept:

  Logical. Include intercept in plot (default: FALSE).

- show_null_range:

  Logical. Show null range region if available (default: TRUE).

- null_value:

  Numeric. Value for the null reference line (default: 0).

- output:

  Character. Either "density" (default) for density slabs or
  "pointinterval" for point + interval display.

- split_color_by_null:

  Logical. If TRUE, posterior densities are split and coloured based on
  whether values fall above or below the null value (default: FALSE).

- color_posterior:

  Color for coefficient posterior densities (default: "dodgerblue").

- color_posterior_outline:

  Color for posterior outlines ( default: "blue").

- color_favours_positive:

  Colour for density regions \> null when `split_color_by_null = TRUE`
  (default: "dodgerblue").

- color_favours_negative:

  Colour for density regions \< null when `split_color_by_null = TRUE`
  (default: "firebrick").

- color_null_range:

  Color for null range shading (default: "grey50").

- color_pointinterval:

  Color for point intervals when `output = "pointinterval"` (default:
  "blue").

- xlim:

  Numeric vector of length 2 specifying x-axis limits. Default is NULL
  (auto-scaled).

- x_breaks:

  Numeric vector specifying custom x-axis break points. Default is NULL
  (auto).

- xlab:

  Character. X-axis label (default: "Coefficient Estimate").

- title:

  Character. Plot title (default: "Meta-Regression Coefficients").

- subtitle:

  Character. Plot subtitle. Default shows null range if provided.

- add_table:

  Logical. Add a table with coefficient summaries on the right (default:
  TRUE).

- table_width:

  Numeric. Relative width of the table vs plot (default: 0.4).

- font:

  Character. Font family for text elements (default: NULL).

- ...:

  Additional arguments (unused).

## Value

A ggplot object (or patchwork object if `add_table = TRUE`).

## Details

This function creates a forest-style visualisation for meta-regression
coefficients that matches the aesthetic of
[`bayes_forest()`](https://blmoran.github.io/bayesma/reference/bayes_forest.md).
Each coefficient is displayed as a posterior density slab (or point
interval), allowing visualisation of the full posterior distribution
rather than just point estimates and credible intervals

When `split_color_by_null = TRUE`, the density is split at the null
value and coloured to show the proportion of the posterior on each side.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- meta_reg(data, studyvar = "author", yi = "yi", vi = "vi",
                mods = ~ year + quality)

# Get evidence summary
ev <- coefficient_evidence(fit, null_range = c(-0.1, 0.1))

# Basic forest-style plot
metareg_mod_plot(ev, model = fit)

# With split coloring
metareg_mod_plot(ev, model = fit, split_color_by_null = TRUE)

# Point interval style
metareg_mod_plot(ev, model = fit, output = "pointinterval")
} # }
```
