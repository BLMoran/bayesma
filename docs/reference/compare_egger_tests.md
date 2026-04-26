# Compare Bayesian and Frequentist Egger Tests

Runs both the Bayesian Egger test (Shi et al.) and the standard
frequentist Egger regression for comparison.

## Usage

``` r
compare_egger_tests(
  data,
  study,
  event_ctrl = NULL,
  event_int = NULL,
  n_ctrl,
  n_int,
  mean_ctrl = NULL,
  mean_int = NULL,
  sd_ctrl = NULL,
  sd_int = NULL,
  likelihood = c("binomial", "gaussian", "poisson"),
  sig_level = 0.1,
  ...
)
```

## Arguments

- data:

  A data frame containing the study data.

- study:

  Character. Column name for study identifiers.

- event_ctrl, event_int:

  Character. Column names for event counts in control and intervention
  groups. Required for binary outcomes.

- n_ctrl, n_int:

  Character. Column names for sample sizes in control and intervention
  groups.

- mean_ctrl, mean_int:

  Character. Column names for means in control and intervention groups.
  Required for continuous outcomes.

- sd_ctrl, sd_int:

  Character. Column names for standard deviations. Required for
  continuous outcomes.

- likelihood:

  Character. `"binomial"` (default), `"gaussian"`, or `"poisson"`.

- sig_level:

  Numeric. Significance level for frequentist test (default: 0.10).

- ...:

  Additional arguments passed to
  [`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html).

## Value

A list with Bayesian and frequentist results.
