# Prior predictive checks

## Introduction

A prior predictive check simulates data from the model using only the
prior distributions — before any observed data are incorporated. The
resulting prior predictive distribution tells you what data the model
considers plausible before seeing the evidence.

Prior predictive checks serve two purposes:

1.  **Prior calibration**: Does the prior predictive distribution cover
    the range of effects that would be scientifically meaningful? Are
    the implied effects absurdly large or implausibly small?
2.  **Model plausibility**: Before fitting, does the model structure
    make sense? Can it generate data that looks anything like what you
    expect to observe?

## Running prior predictive checks

``` r
fit_prior <- bayesma(
  data,
  model_type    = "random_effect",
  sample_prior  = TRUE
)

ecdf_prior_plot(fit_prior)
```

`sample_prior = TRUE` draws from the prior predictive distribution. No
data are used in the likelihood; the `data` argument only provides the
standard errors $`s_i`$ for the simulation scale.

## What to look for

**Effect scale.** On the log-OR scale, effects above 2 or below −2 are
unusual in most clinical literatures. If the prior predictive
distribution assigns substantial probability to $`|\mu| > 3`$, the prior
on $`\mu`$ is too vague.

**Heterogeneity scale.** A prior that generates $`\tau > 2`$ on the
log-OR scale implies studies differ by more than 4 log-ORs — implausible
in most applications. Tighten the $`\tau`$ prior if the prior predictive
places more than 10% of mass above $`\tau = 1`$.

**Study-level effects.** The prior predictive for individual study
effects ($`\theta_i`$) should span a plausible range without being
dominated by the prior tails.

## Interpreting the ECDF prior plot

[`ecdf_prior_plot()`](https://blmoran.github.io/bayesma/reference/ecdf_prior_plot.md)
overlays the empirical CDF of the observed effects with the prior
predictive CDF. If the observed CDF falls far outside the prior
predictive band, the priors are inconsistent with the data — the
posterior will be dominated by the likelihood pulling against the prior,
which can cause slow mixing and poor calibration.

## Adjusting priors

If the prior predictive check reveals that priors are too vague or too
tight:

``` r
fit_adj <- bayesma(
  data,
  model_type   = "random_effect",
  mu_prior     = normal(0, 0.5),
  tau_prior    = half_normal(0, 0.25),
  sample_prior = TRUE
)

ecdf_prior_plot(fit_adj)
```

Iterate until the prior predictive covers the plausible effect range
without being dominated by the tails.

## Informative priors from external data

Empirically-derived priors on $`\tau`$ are available from Turner et
al. (2015) (medical treatments, binary outcomes) and Rhodes et
al. (2015) (psychological interventions, continuous outcomes). Using
domain-calibrated priors avoids the need to rely entirely on prior
predictive simulation:

``` r
bayesma(
  data,
  model_type = "random_effect",
  tau_prior  = half_normal(0, 0.32)  # Turner et al. median for pharmacological treatments
)
```

## Prior predictive checks for publication bias models

For models with selection mechanisms (selection weight, Copas), prior
predictive checks should verify that the implied selection probabilities
are plausible. A prior that generates weight functions concentrated near
zero (no selection) or near one (all studies published regardless of
result) is uninformative.

``` r
ecdf_prior_plot(fit_selection_prior, parameter = "selection_weights")
```
