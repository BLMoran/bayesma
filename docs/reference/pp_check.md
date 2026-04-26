# Posterior Predictive Check for bayesma Objects

Produces a posterior predictive check plot comparing observed data to
replicated datasets drawn from the posterior predictive distribution.
Mimics the interface of
[`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html).

## Usage

``` r
pp_check(object, type = "dens_overlay", ndraws = 100, ...)
```

## Arguments

- object:

  A `bayesma` object.

- type:

  Character. Plot type passed to `bayesplot::ppc_*`. One of
  `"dens_overlay"` (default), `"hist"`, `"stat"`, `"scatter"`, `"bars"`
  (discrete only), `"ecdf_overlay"`, `"ribbon"`.

- ndraws:

  Number of posterior draws to use. Default 100.

- ...:

  Additional arguments passed to the underlying `bayesplot::ppc_*`
  function.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- bayesma(data, likelihood = "binomial", ...)
pp_check(fit)
pp_check(fit, type = "stat", stat = "mean")
pp_check(fit, type = "scatter")
} # }
```
