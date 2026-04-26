# Model fitting and convergence diagnostics

## Introduction

MCMC inference is only valid when the chains have converged to the
target posterior. Convergence failure produces unreliable estimates,
overconfident credible intervals, and incorrect Bayes factors. Every
**bayesma** fit should pass convergence diagnostics before results are
reported.

## Accessing diagnostics

``` r
fit <- bayesma(data, model_type = "random_effect")
diagnostics(fit)
```

[`diagnostics()`](https://blmoran.github.io/bayesma/reference/diagnostics.md)
returns a summary of all convergence statistics. Key statistics are also
printed automatically when `quiet = FALSE` (the default).

## $`\hat{R}`$ (potential scale reduction factor)

$`\hat{R}`$ compares the variance of a parameter within chains to its
variance across chains. A value close to 1 indicates the chains have
mixed. The threshold is $`\hat{R} < 1.01`$ for reliable inference
(Vehtari et al., 2021). Values above 1.05 indicate mixing failure.

``` r
diagnostics(fit) |>
  dplyr::filter_out(rhat < 1.01)
```

## Effective sample size (ESS)

ESS estimates the number of independent samples in the correlated MCMC
chain. Two variants:

- **Bulk ESS**: measures mixing in the centre of the distribution
  (relevant for means and medians).
- **Tail ESS**: measures mixing in the tails (relevant for 95% credible
  intervals).

Thresholds: bulk ESS \> 400 per parameter; tail ESS \> 400 per
parameter. Low tail ESS for $`\tau`$ is common when $`k`$ is small and
the posterior for $`\tau`$ is near-zero.

## Divergent transitions

Divergent transitions occur when the HMC integrator fails to track the
posterior geometry. They indicate a pathological posterior, often caused
by a funnel-shaped geometry near $`\tau = 0`$.

``` r
diagnostics(fit)$divergences
```

**Zero divergences is the target.** Any divergences require
investigation before the fit can be trusted.

Common solutions:

| Problem | Solution |
|----|----|
| Funnel near $`\tau \approx 0`$ | NCP already applied (default); tighten $`\tau`$ prior |
| Complex model | Increase `adapt_delta` to 0.99 |
| Collinear parameters | Reparameterise or centre covariates |
| Model misspecification | Simplify model or check data |

## Trace plots

Visual inspection of trace plots catches mixing failures that summary
statistics miss:

``` r
bayesplot::mcmc_trace(fit$draws(), pars = c("mu", "tau"))
```

Well-mixing chains look like “hairy caterpillars” — rapid, uncorrelated
oscillation across the full posterior range. Signs of poor mixing
include:

- **Stuck chains**: one or more chains staying in one region for many
  iterations.
- **Slow mixing**: autocorrelated, sluggish movement across the
  parameter space.
- **Bimodal chains**: abrupt jumps between two regions.

## Pair plots

Pair plots reveal correlations between parameters and can show the
characteristic funnel geometry:

``` r
bayesplot::mcmc_pairs(fit$draws(), pars = c("mu", "tau"), off_diag_args = list(size = 0.5))
```

A banana-shaped joint distribution for $`(\mu, \tau)`$ near $`\tau = 0`$
is expected and harmless when the NCP is used. A tight funnel that the
chains fail to explore is the warning sign.

## LOO Pareto $`\hat{k}`$ diagnostics

After fitting, LOO-IC computes per-observation Pareto $`\hat{k}`$
statistics. Values above 0.7 flag observations whose LOO contribution is
unreliable:

``` r
fit$loo()$diagnostics
```

Studies with $`\hat{k} > 0.7`$ are influential. Refit without these
studies (LOSO-CV) to verify that the pooled estimate is robust.

## Reporting convergence

In any published meta-analysis using **bayesma**, report:

- Number of chains and iterations (warmup + sampling).
- `adapt_delta` value used.
- Maximum $`\hat{R}`$ across all parameters.
- Minimum bulk and tail ESS.
- Number of divergent transitions.
- Any parameter with $`\hat{R} > 1.01`$ or ESS \< 400, and the
  resolution taken.
