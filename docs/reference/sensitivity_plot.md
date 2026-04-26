# Generate Sensitivity Analysis Plot for Bayesian Meta-Analysis

Creates a sensitivity analysis visualisation showing how meta-analytic
estimates vary across different model strategies and user-specified
priors.

For RoBMA results, models must be pre-computed using
`run_robma_sensitivity()` and attached to the bayesma object before
calling this function.

## Usage

``` r
sensitivity_plot(
  model,
  data,
  priors,
  measure,
  prior_order = NULL,
  model_order = NULL,
  rob_var = NULL,
  exclude_high_rob = FALSE,
  incl_common_effect = FALSE,
  incl_random_effect = TRUE,
  incl_bias_corrected = FALSE,
  incl_selection_copas = FALSE,
  incl_selection_weight = FALSE,
  incl_pet_peese = FALSE,
  incl_robust = FALSE,
  incl_robma = FALSE,
  parallel = FALSE,
  workers = NULL,
  seed = TRUE,
  add_probs = FALSE,
  null_value = NULL,
  null_range = NULL,
  add_null_range = FALSE,
  color_null_range = "#77bb41",
  label_control = "Control",
  label_intervention = "Intervention",
  title = NULL,
  subtitle = NULL,
  title_align = "left",
  xlim = NULL,
  x_breaks = NULL,
  color_palette = NULL,
  color_overall_posterior = "dodgerblue",
  color_overall_posterior_outline = "blue",
  split_color_by_null = FALSE,
  color_favours_control = "firebrick",
  color_favours_intervention = "dodgerblue",
  plot_width = 4,
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

- prior_order:

  Optional character vector specifying the display order of priors.
  Should contain the names (IDs) of the priors in the desired order,
  e.g., `c("vague", "weak_reg", "informative")`. If NULL (default),
  priors are displayed in the order they appear in the `priors` list.

- model_order:

  Optional. Specifies the display order of model sections. Can be
  provided as unquoted names (using rlang) or as a character vector.

  Valid values match the `incl_*` argument names (without the `incl_`
  prefix): `common_effect`, `random_effect`, `bias_corrected`,
  `selection_copas`, `selection_weight`, `pet_peese`, `robust`, `robma`.

  For RoBMA, use `robma` to include both conditional and model-averaged
  estimates (conditional first), or `robma_conditional` /
  `robma_model_averaged` separately.

  Example:
  `model_order = c(common_effect, random_effect, pet_peese, robma)`

  If NULL (default), sections are displayed in the default order.

- rob_var:

  Optional. Name of the risk-of-bias variable (unquoted).

- exclude_high_rob:

  Logical. If TRUE, runs an "Excluding High RoB" section.

- incl_common_effect, incl_random_effect, incl_bias_corrected,
  incl_selection_copas, incl_selection_weight, incl_pet_peese,
  incl_robust:

  Logical. Which model strategy sections to include.

- incl_robma:

  Logical. If TRUE, include RoBMA section. Requires
  `model$robma_sensitivity` to be present (see Details).

- parallel:

  Logical. If TRUE, uses parallel processing for non-RoBMA refits.

- workers:

  Optional integer. Number of parallel workers.

- seed:

  Logical. If TRUE (default), uses parallel-safe seeding.

- add_probs:

  Logical. Add probability columns to the results table.

- null_value, null_range, add_null_range, color_null_range:

  ROPE settings.

- label_control, label_intervention:

  Group labels for the plot.

- title, subtitle, title_align:

  Title settings.

- xlim, x_breaks:

  Density axis settings.

- color_palette, color_overall_posterior,
  color_overall_posterior_outline:

  Colour settings.

- split_color_by_null:

  Logical. If TRUE, colour the posterior split at the null value using
  `color_favours_control` and `color_favours_intervention`.

- color_favours_control:

  Colour for the side of the posterior favouring control.

- color_favours_intervention:

  Colour for the side favouring intervention.

- plot_width:

  Width ratio for the density plot section.

- font:

  Optional font family.

## Value

A `bayesma_sensitivity_plot` object (patchwork combining tables and
plots).

## Details

### RoBMA Results

To include RoBMA in the sensitivity plot, you must pre-compute the RoBMA
fits:

    # Step 1: Run RoBMA sensitivity analysis
    robma_sens <- run_robma_sensitivity(
      data = my_data,
      priors = my_priors,
      robma_template = my_robma_fit,
      parallel = TRUE
    )

    # Step 2: Attach to bayesma model
    model <- attach_robma_sensitivity(model, robma_sens)

    # Step 3: Create sensitivity plot with RoBMA
    sensitivity_plot(model, data, priors, measure = "OR", incl_robma = TRUE)

This separation allows the computationally expensive RoBMA fitting to be
done once and cached, rather than being re-run every time the plot is
generated.
