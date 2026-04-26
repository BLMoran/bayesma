# ECDF Plot Comparing Prior Sensitivity

Creates an empirical cumulative distribution function (ECDF) plot
comparing posterior distributions across different prior specifications
for a single model type (or small set of model types).

Use this plot to assess how sensitive conclusions are to the choice of
prior distributions.

## Usage

``` r
ecdf_prior_plot(
  model,
  data,
  priors,
  measure,
  model_types = "random_effect",
  prior_order = NULL,
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
  linetype_by_model = TRUE,
  show_density = TRUE,
  font = NULL
)
```

## Arguments

- model:

  A fitted `bayesma` object.

- data:

  A data frame containing the study data used to fit the model.

- priors:

  A named list of prior specifications. Each element must be a list with
  at least `mu_prior` and optionally `tau_prior`, and may include `name`
  for display labels.

- measure:

  Effect measure string (e.g., "OR", "RR", "HR", "IRR", "MD", "SMD").

- model_types:

  Character vector specifying which model type(s) to include. Maximum of
  2 model types to avoid clutter. Valid values: `"common_effect"`,
  `"random_effect"`, `"bias_corrected"`, `"selection_copas"`,
  `"selection_weight"`, `"pet_peese"`, `"robust"`, `"robma"`. Default is
  `"random_effect"`.

- prior_order:

  Optional character vector specifying the display order of priors.
  Should contain the names (IDs) of the priors in the desired order. If
  NULL, priors are displayed in the order they appear in `priors`.

- prob_reference:

  Character string. What reference to use for probability axis labels.
  Either `"null"` (compare to `null_value`) or `"null_range"` (compare
  to `null_range` boundaries). Default is `"null"`.

- null_value:

  Null hypothesis value. If NULL, uses measure default.

- null_range:

  Numeric vector of length 2 giving null range bounds.

- add_null_range:

  Logical. If TRUE and `null_range` is NULL, uses measure-appropriate
  defaults.

- color_null_range:

  Fill colour for the null range band. Default `"#77bb41"`.

- label_control:

  Label for control group. Default `"Control"`.

- label_intervention:

  Label for intervention group. Default `"Intervention"`.

- title:

  Optional plot title.

- subtitle:

  Optional plot subtitle. If NULL and a single model type is selected,
  displays the model type.

- xlim:

  Optional x-axis limits.

- x_breaks:

  Optional x-axis break points.

- color_palette:

  Optional named vector of colours for priors.

- linetype_by_model:

  Logical. If TRUE and multiple model types are selected, uses different
  linetypes for each model type. Default TRUE.

- show_density:

  Logical. If TRUE, includes a density plot below the ECDF. Default is
  TRUE.

- font:

  Optional font family.

## Value

A ggplot object (or patchwork object if `show_density = TRUE`).

## Details

The plot shows one ECDF line per prior specification (and per model type
if multiple are selected). This allows direct comparison of how
different prior choices affect the posterior distribution.

When two model types are selected, colour represents prior and linetype
represents model type, allowing comparison of both dimensions
simultaneously.

The left y-axis shows P(effect \< x) and the right y-axis shows P(effect
\> x).

## See also

[`ecdf_model_plot`](https://blmoran.github.io/bayesma/reference/ecdf_model_plot.md)
for comparing models within a prior.
