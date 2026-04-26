# Create posterior plots for Bayesian meta-analysis

Produces density plots of the pooled effect (\\\mu\\), heterogeneity
(\\\tau\\), and/or prediction distribution from a fitted bayesma model.

## Usage

``` r
overall_plot(
  model,
  estimand,
  add_tau = FALSE,
  add_mu_prior = FALSE,
  add_tau_prior = FALSE,
  add_pred = FALSE,
  color_pred_posterior = "forestgreen",
  color_pred_posterior_outline = "darkgreen",
  null_value = NULL,
  null_range = NULL,
  add_null_range = FALSE,
  color_null_range = "#77bb41",
  label_control = "Control",
  label_intervention = "Intervention",
  title = NULL,
  subtitle = NULL,
  title_align = "left",
  mu_xlim = NULL,
  tau_xlim = NULL,
  pred_xlim = NULL,
  x_breaks = NULL,
  tau_breaks = NULL,
  pred_breaks = NULL,
  color_overall_posterior = "dodgerblue",
  color_overall_posterior_outline = "blue",
  split_color_by_null = FALSE,
  color_favours_control = "firebrick",
  color_favours_intervention = "dodgerblue",
  tau_posterior_color = "#77bb41",
  tau_posterior_outline = NULL,
  tau_slab_scale = 0.7,
  font = NULL,
  plot_arrangement = "vertical"
)
```

## Arguments

- model:

  A `bayesma` object.

- estimand:

  Character string for the effect measure or marginal estimand.
  Relative-effect: `"OR"`, `"HR"`, `"RR"`, `"IRR"`, `"MD"`, `"SMD"`.
  Marginal: `"RD"`/`"ARR"`, `"ATE"`, `"ATT"`, `"CATE"`. For marginal
  estimands the pooled posterior is drawn from `model$marginal$draws`.

- add_tau:

  Logical, whether to include the heterogeneity (tau) plot.

- add_mu_prior:

  Logical, whether to include the prior distribution for mu.

- add_tau_prior:

  Logical, whether to include the prior distribution for tau.

- add_pred:

  Logical, whether to include the prediction distribution plot.

- color_pred_posterior:

  Color for the prediction posterior distribution.

- color_pred_posterior_outline:

  Outline color for the prediction posterior distribution.

- null_value:

  Numeric, the null value for the effect (default from
  get_measure_properties).

- null_range:

  Numeric vector of length 2, range for null region (e.g., c(0.9, 1.1)).

- add_null_range:

  Logical, whether to add shaded null range region.

- color_null_range:

  Color for the null range shading.

- label_control:

  Character string for control group label.

- label_intervention:

  Character string for intervention group label.

- title:

  Plot title.

- subtitle:

  Plot subtitle.

- title_align:

  Alignment for title ("left", "center", "right").

- mu_xlim:

  Numeric vector of length 2 for mu plot x-axis limits.

- tau_xlim:

  Numeric vector of length 2 for tau plot x-axis limits.

- pred_xlim:

  Numeric vector of length 2 for prediction plot x-axis limits.

- x_breaks:

  Numeric vector for x-axis breaks on mu plot.

- tau_breaks:

  Numeric vector for x-axis breaks on tau plot.

- pred_breaks:

  Numeric vector for x-axis breaks on prediction plot.

- color_overall_posterior:

  Color for the overall posterior distribution.

- color_overall_posterior_outline:

  Outline color for the posterior distribution.

- split_color_by_null:

  Logical, whether to split posterior color by null value.

- color_favours_control:

  Color for region favoring control.

- color_favours_intervention:

  Color for region favoring intervention.

- tau_posterior_color:

  Base color for tau posterior.

- tau_posterior_outline:

  Outline color for tau posterior.

- tau_slab_scale:

  Numeric, scaling factor for tau slab height relative to mu (default
  0.7).

- font:

  Font family for text elements.

- plot_arrangement:

  Character, "vertical" (stacked) or "horizontal" (side by side).

## Value

A ggplot object (or combined plot if add_tau = TRUE or add_pred = TRUE).
