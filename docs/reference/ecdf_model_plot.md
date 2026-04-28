# ECDF Plot Comparing Model Strategies

Creates an empirical cumulative distribution function (ECDF) plot
comparing posterior distributions across different model strategies
(e.g., random effects, selection models, PET-PEESE) under a single prior
specification.

Use this plot to assess how robust conclusions are across different
modelling assumptions.

## Usage

``` r
ecdf_model_plot(
  model,
  data,
  prior,
  estimand,
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
  font = NULL
)
```

## Arguments

- model:

  A fitted `bayesma` object.

- data:

  A data frame containing the study data used to fit the model.

- prior:

  A list containing the prior specification with elements:

  `mu_prior`

  :   (Required) Prior specification for mu.

  `tau_prior`

  :   (Optional) Prior specification for tau.

  `name`

  :   Display label for the prior. Required when `incl_robma = TRUE` as
      it is used to match against priors in `model$robma_sensitivity`.
      If not provided and `incl_robma = FALSE`, defaults to
      "User-specified prior".

- estimand:

  Effect estimand string (e.g., "OR", "RR", "HR", "IRR", "MD", "SMD").

- rob_var:

  Optional. Name of the risk-of-bias variable (unquoted).

- exclude_high_rob:

  Logical. If TRUE, includes an "Excluding High RoB" section.

- incl_common_effect:

  Logical. Include common effect model. Default FALSE.

- incl_random_effect:

  Logical. Include random effects model. Default TRUE.

- incl_bias_corrected:

  Logical. Include bias-corrected model. Default FALSE.

- incl_selection_copas:

  Logical. Include Copas selection model. Default FALSE.

- incl_selection_weight:

  Logical. Include weight-function selection model. Default FALSE.

- incl_pet_peese:

  Logical. Include PET-PEESE model. Default FALSE.

- incl_robust:

  Logical. Include robust mixture model. Default FALSE.

- incl_robma:

  Logical. Include RoBMA model. Requires `model$robma_sensitivity` to be
  present and `prior$name` to match one of the priors used when fitting.
  Default FALSE.

- prob_reference:

  Character string. What reference to use for probability axis labels.
  Either `"null"` (compare to `null_value`) or `"null_range"` (compare
  to `null_range` boundaries). Default is `"null"`.

- null_value:

  Null hypothesis value. If NULL, uses estimand default.

- null_range:

  Numeric vector of length 2 giving null range bounds.

- add_null_range:

  Logical. If TRUE and `null_range` is NULL, uses estimand-appropriate
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

  Optional plot subtitle. If NULL and a prior name is available,
  displays the prior name.

- xlim:

  Optional x-axis limits.

- x_breaks:

  Optional x-axis break points.

- color_palette:

  Optional named vector of colours for model sections.

- show_density:

  Logical. If TRUE, includes a density plot below the ECDF. Default is
  TRUE.

- font:

  Optional font family.

## Value

A ggplot object (or patchwork object if `show_density = TRUE`).

## Details

The plot shows one ECDF line per model strategy, all using the same
prior. This allows direct comparison of how different modelling
assumptions (random effects, selection models, bias correction, etc.)
affect the posterior distribution.

The left y-axis shows P(effect \< x) and the right y-axis shows P(effect
\> x).

## See also

[`ecdf_prior_plot`](https://blmoran.github.io/bayesma/reference/ecdf_prior_plot.md)
for comparing priors within a model type.

## Examples

``` r
if (FALSE) { # \dontrun{
priors <- list(
  vague = list(
    name = "Vague",
    mu_prior = normal(0, 10),
    tau_prior = half_cauchy(0, 1)
  ),
  weak_reg = list(
    name = "Weakly Regularising",
    mu_prior = normal(0, 1),
    tau_prior = half_cauchy(0, 0.5)
  )
)

ecdf_model_plot(
  model = model,
  data = dat,
  estimand = "OR",
  prior = priors$weak_reg,
  incl_random_effect = TRUE,
  incl_pet_peese = TRUE
)

ecdf_model_plot(
  model = model,
  data = dat,
  estimand = "OR",
  prior = priors$weak_reg,
  incl_robma = TRUE
)
} # }
```
