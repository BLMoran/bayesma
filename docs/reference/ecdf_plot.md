# Create ECDF Plot for Sensitivity Analysis

Creates an empirical cumulative distribution function (ECDF) plot
showing posterior distributions across different sensitivity analyses.
Accepts both `brmsfit` and `bayesma` objects.

For bayesma objects the function refits the model across the requested
`model_type` variants using `.build_bayesma_sensitivity_draws()`,
exactly as in
[`sensitivity_plot`](https://blmoran.github.io/bayesma/reference/sensitivity_plot.md).

## Usage

``` r
ecdf_plot(
  model,
  data,
  priors,
  measure,
  study_var = NULL,
  rob_var = NULL,
  exclude_high_rob = FALSE,
  incl_common_effect = FALSE,
  incl_random_effect = TRUE,
  incl_bias_corrected = FALSE,
  incl_selection_copas = FALSE,
  incl_selection_weight = FALSE,
  incl_pet_peese = FALSE,
  pet_peese_direction = "negative",
  pet_peese_threshold = 0.1,
  incl_robust = FALSE,
  incl_mixture = FALSE,
  incl_bma = FALSE,
  model_bma = NULL,
  prob_reference = "null",
  null_value = NULL,
  null_range = NULL,
  add_null_range = FALSE,
  color_null_range = "#77bb41",
  label_control = "Control",
  label_intervention = "Intervention",
  title = NULL,
  subtitle = NULL,
  xlim = NULL,
  x_breaks = NULL,
  color_palette = NULL,
  show_density = TRUE,
  prior_to_plot = "weakreg",
  font = NULL
)
```

## Arguments

- model:

  A fitted `bayesma` object. For RoBMA results, must have
  `$robma_sensitivity` attached via
  [`attach_robma_sensitivity()`](https://blmoran.github.io/bayesma/reference/attach_robma_sensitivity.md).

- data:

  A data frame containing the study data used to fit the model.

- priors:

  A named list of prior specifications. Each element must be a list with
  at least `mu_prior` and (optionally) `tau_prior`, and may include
  `name` used for display.

- measure:

  Effect measure string (e.g., "OR", "RR", "HR", "IRR", "MD", "SMD").

- rob_var:

  Optional. Name of the risk-of-bias variable (unquoted).

- exclude_high_rob:

  Logical. If TRUE, runs an "Excluding High RoB" section.

- incl_common_effect, incl_random_effect, incl_bias_corrected,
  incl_selection_copas, incl_selection_weight, incl_pet_peese,
  incl_robust:

  Logical. Which model strategy sections to include.

- prob_reference:

  Character string. What reference to use for probability calculations.
  Either `"null"` (compare to `null_value`) or `"null_range"` (compare
  to `null_range` boundaries). Default is `"null"`.

- null_value, null_range, add_null_range, color_null_range:

  ROPE settings.

- label_control, label_intervention:

  Group labels for the plot.

- xlim, x_breaks:

  Density axis settings.

- show_density:

  Logical. If `TRUE`, includes a density plot below the ECDF. Default is
  `TRUE`.

- prior_to_plot:

  Character string. Which prior specification to plot in the ECDF. One
  of `"vague"`, `"weakreg"`, or `"informative"`. Default is `"weakreg"`.

- font:

  Optional font family.

## Value

A ggplot object (or patchwork object if `show_density = TRUE`).

## Details

The output ECDF plot includes a left y-axis showing P(effect \< x) and a
right y-axis showing P(effect \> x), with optional ROPE shading and an
optional density panel below.
