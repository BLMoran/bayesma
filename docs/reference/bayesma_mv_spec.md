# Build a bivariate meta-analysis specification object

Build a bivariate meta-analysis specification object

## Usage

``` r
bayesma_mv_spec(
  data,
  studyvar,
  mean_ctrl_1,
  mean_int_1,
  sd_ctrl_1,
  sd_int_1,
  n_ctrl_1,
  n_int_1,
  mean_ctrl_2,
  mean_int_2,
  sd_ctrl_2,
  sd_int_2,
  n_ctrl_2,
  n_int_2,
  outcome_labels = c("outcome_1", "outcome_2"),
  likelihood = c("gaussian"),
  stage = c("two_stage", "one_stage"),
  rho_within = 0.5,
  mu_prior = NULL,
  tau_prior = NULL,
  rho_between_prior = NULL,
  custom_model = NULL,
  custom_data = NULL
)
```

## Arguments

- data:

  A data frame with one row per study containing arm-level data for both
  outcomes.

- studyvar:

  Column for study identifiers (unquoted).

- mean_ctrl_1, mean_int_1:

  Columns for control and intervention means for outcome 1 (unquoted).

- sd_ctrl_1, sd_int_1:

  Columns for control and intervention SDs for outcome 1 (unquoted).

- n_ctrl_1, n_int_1:

  Columns for control and intervention sample sizes for outcome 1
  (unquoted).

- mean_ctrl_2, mean_int_2:

  Columns for control and intervention means for outcome 2 (unquoted).

- sd_ctrl_2, sd_int_2:

  Columns for control and intervention SDs for outcome 2 (unquoted).

- n_ctrl_2, n_int_2:

  Columns for control and intervention sample sizes for outcome 2
  (unquoted).

- outcome_labels:

  Character vector of length 2 with labels for the two outcomes.
  Default: `c("outcome_1", "outcome_2")`.

- likelihood:

  Character. Currently only `"gaussian"`.

- stage:

  Character. `"two_stage"` (effect sizes computed then modelled) or
  `"one_stage"` (marginalised model).

- rho_within:

  Numeric scalar in `[-1, 1]`. Within-study correlation between the two
  outcomes, assumed known. Default: `0.5`.

- mu_prior:

  Prior on pooled effects. Either a single prior (applied to both
  outcomes) or a named list with elements matching `outcome_labels`.

- tau_prior:

  Prior on between-study SDs. Either a single prior or a named list.

- rho_between_prior:

  Prior on the between-study correlation. Default: `uniform(-1, 1)`.

- custom_model:

  Optional character scalar containing complete Stan code to override
  the generated program.

- custom_data:

  Optional named list merged into the Stan data list.

## Value

An object of class `"bayesma_mv_spec"`.
