# Single-Page Model Diagnostics for bayesma Objects

Produces a single-page diagnostic summary for a Bayesian meta-analysis
model fitted via `cmdstanr`. Combines visual diagnostics (trace plots,
posterior predictive check, Rhat, ESS) with a tabular summary of MCMC
health indicators. All output is arranged on a single page using
`patchwork`.

## Usage

``` r
diagnostics(object, ...)

# S3 method for class 'bayesma'
diagnostics(object, pars = NULL, ndraws = 100, ...)
```

## Arguments

- object:

  A `bayesma` object containing a `CmdStanMCMC` fit in `object$fit`.

- ...:

  Additional arguments (currently unused).

- pars:

  Character vector of parameter names for trace/ACF plots. If `NULL`,
  sensible defaults are chosen based on the model stage.

- ndraws:

  Number of posterior draws for the posterior predictive check. Default
  is 100.

## Value

Invisibly returns a list of diagnostic values. Prints a composite
`patchwork` plot as a side effect.

## Details

The diagnostic page contains six panels:

- Top-left:

  Trace plots for key parameters.

- Top-right:

  Posterior predictive check.

- Middle-left:

  Rhat values for all parameters.

- Middle-right:

  Effective sample size ratios for all parameters.

- Bottom-left:

  Autocorrelation function for key parameters.

- Bottom-right:

  MCMC diagnostics summary table.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- bayesma(data, likelihood = "binomial", measure = "OR", ...)
bayesma_diagnostics(fit)
} # }
```
