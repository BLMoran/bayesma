# Posterior and Prior Predictive Checks for bayesma Objects

Compares observed data to replicated datasets drawn from either the
posterior predictive distribution (standard use) or the prior predictive
distribution (when the model was fitted with `sample_prior = TRUE`).

**Posterior predictive check** (`sample_prior = FALSE`, the default):
draws `y_rep` from `p(y_rep | y)`. Agreement between `y` and `y_rep`
indicates the fitted model can reproduce the observed data.

**Prior predictive check** (`sample_prior = TRUE`): draws `y_rep` from
`p(y_rep)`, integrating over the prior without conditioning on the data.
This is a tool for prior *elicitation and sanity-checking* — verifying
that the prior does not place substantial mass on data values that are
impossible or implausible given domain knowledge, before the data have
been seen.

## Usage

``` r
pp_check(object, type = "dens_overlay", ndraws = 100, ...)
```

## Arguments

- object:

  A `bayesma` object. If fitted with `sample_prior = TRUE`, a prior
  predictive check is produced; otherwise a posterior predictive check.

- type:

  Character. Plot type passed to `bayesplot::ppc_*`. One of
  `"dens_overlay"` (default), `"hist"`, `"stat"`, `"scatter"`, `"bars"`
  (discrete only), `"ecdf_overlay"`, `"ribbon"`.

- ndraws:

  Integer. Number of draws to display. Default `100`.

- ...:

  Additional arguments passed to the underlying `bayesplot::ppc_*`
  function.

## Value

A `ggplot` object.

## Prior predictive checks and double-dipping

It is legitimate to use a prior predictive check to verify that a prior
is coherent on the observable scale (e.g. that it does not imply
impossible event counts). It is **not** legitimate to iteratively adjust
priors until `y_rep` matches the observed `y`: doing so smuggles the
data into the prior, inflating posterior confidence. Priors should be
specified from external knowledge, expert elicitation, or independent
reference data — not from the analysis dataset itself.

## Examples

``` r
if (FALSE) { # \dontrun{
# Posterior predictive check
fit <- bayesma(data, likelihood = "binomial", ...)
pp_check(fit)
pp_check(fit, type = "stat", stat = "mean")

# Prior predictive check — use to verify priors are sensible,
# not to calibrate them against the analysis data
prior_fit <- bayesma(data, likelihood = "binomial", ..., sample_prior = TRUE)
pp_check(prior_fit)
} # }
```
