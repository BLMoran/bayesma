# Build an Egger's test specification object

Build an Egger's test specification object

## Usage

``` r
egger_spec(
  data,
  studyvar,
  n_ctrl,
  n_int,
  event_ctrl = NULL,
  event_int = NULL,
  mean_ctrl = NULL,
  mean_int = NULL,
  sd_ctrl = NULL,
  sd_int = NULL,
  likelihood = c("binomial", "gaussian", "poisson"),
  heterogeneity = c("multiplicative", "additive"),
  alpha_prior = NULL,
  beta_prior = NULL,
  kappa_prior = NULL,
  gamma_prior = NULL,
  d_prior = NULL,
  tau_prior = NULL,
  credible_level = 0.9,
  custom_model = NULL,
  custom_data = NULL
)
```

## Arguments

- data:

  A data frame with one row per study.

- studyvar:

  Character. Column name of the study identifier.

- n_ctrl, n_int:

  Character. Column names of control and intervention sample sizes.

- event_ctrl, event_int:

  Character. Column names of event counts (binomial / Poisson
  likelihoods).

- mean_ctrl, mean_int, sd_ctrl, sd_int:

  Character. Column names of arm means and SDs (Gaussian likelihood).

- likelihood:

  Character. One of `"binomial"`, `"gaussian"`, `"poisson"`.

- heterogeneity:

  Character. `"multiplicative"` (default) or `"additive"`.

- alpha_prior:

  Prior on the intercept.

- beta_prior:

  Prior on the slope (the Egger coefficient).

- kappa_prior:

  Prior on the multiplicative heterogeneity coefficient.

- gamma_prior:

  Prior on the dispersion parameter.

- d_prior:

  Prior on the overdispersion parameter.

- tau_prior:

  Prior on the between-study SD for the additive heterogeneity
  parameterisation.

- credible_level:

  Numeric in `(0, 1)`. Credible-interval level for the summary. Default
  `0.90`.

- custom_model:

  Optional character scalar of Stan code overriding the generated
  program.

- custom_data:

  Optional named list merged into the Stan data list.

## Value

An object of class `"egger_spec"`.
