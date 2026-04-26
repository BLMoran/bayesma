# Run RoBMA Models Across Multiple Prior Specifications

Pre-computes RoBMA models for each prior specification provided. These
fits are stored and can later be attached to a bayesma object for use in
[`sensitivity_plot()`](https://blmoran.github.io/bayesma/reference/sensitivity_plot.md).
This separation allows the computationally expensive RoBMA fitting to be
done once and reused.

## Usage

``` r
robma_sensitivity(
  data,
  priors,
  robma_template = NULL,
  robma_args = NULL,
  parallel = FALSE,
  workers = NULL,
  seed = TRUE,
  .progress = TRUE
)
```

## Arguments

- data:

  A data frame containing the study data.

- priors:

  A named list of prior specifications. Each element must be a list with
  at least `mu_prior` and (optionally) `tau_prior`, and may include
  `name` used for display. Example:
  `list( prior_1 = list(name = "Vague", mu_prior = normal(0, 10), tau_prior = half_cauchy(0, 1)), prior_2 = list(name = "Regularising", mu_prior = normal(0, 1), tau_prior = half_cauchy(0, 0.5)) )`

- robma_template:

  Either:

  - A fitted RoBMA object (class `bayesma_robma`) with stored
    `$meta$call_args` to use as a template for refitting across priors,
    OR

  - `NULL` (default), in which case `robma_args` must be provided

- robma_args:

  A list of arguments to pass to
  [`robma()`](https://blmoran.github.io/bayesma/reference/robma.md) when
  `robma_template` is `NULL`. Must include all required arguments except
  `priors_effect` and `priors_heterogeneity`, which will be set from the
  `priors` argument.

- parallel:

  Logical. If TRUE, uses parallel processing for refits.

- workers:

  Optional integer. Number of parallel workers.

- seed:

  Logical. If TRUE (default), uses parallel-safe seeding.

- .progress:

  Logical. If TRUE, displays a progress bar.

## Value

A `bayesma_robma_sensitivity` object: a named list of RoBMA fits keyed
by prior IDs, with additional metadata in `$meta`.

## Details

The returned object can be attached to a bayesma model using
[`attach_robma_sensitivity()`](https://blmoran.github.io/bayesma/reference/attach_robma_sensitivity.md)
or passed directly to
[`sensitivity_plot()`](https://blmoran.github.io/bayesma/reference/sensitivity_plot.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Define priors for sensitivity analysis
priors <- list(
  vague = list(
    name = "Vague",
    mu_prior = normal(0, 10),
    tau_prior = half_cauchy(0, 1)
  ),
  informative = list(
    name = "Informative",
    mu_prior = normal(0, 0.5),
    tau_prior = half_cauchy(0, 0.3)
  )
)

# Option 1: Using a template RoBMA fit
robma_sens <- robma_sensitivity(
  data = my_data,
  priors = priors,
  robma_template = my_robma_fit,
  parallel = TRUE
)

# Option 2: Specifying robma arguments directly
robma_sens <- robma_sensitivity(
  data = my_data,
  priors = priors,
  robma_args = list(
    studyvar = "study_id",
    event_ctrl = "events_c",
    event_int = "events_i",
    n_ctrl = "n_c",
    n_int = "n_i",
    likelihood = "binomial"
  ),
  parallel = TRUE
)

# Attach to bayesma model
model <- attach_robma_sensitivity(model, robma_sens)

# Now sensitivity_plot can use the pre-computed RoBMA fits
sensitivity_plot(model, data, priors, measure = "OR", incl_robma = TRUE)
} # }
```
