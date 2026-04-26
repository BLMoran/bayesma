# Gaussian random-effects model

## Introduction

The Gaussian random-effects model is the default meta-analysis model in
**bayesma**. It assumes that the true effects across studies are
exchangeable draws from a normal distribution, with between-study
variability governed by a heterogeneity parameter $`\tau`$. This is the
correct starting point when there is no strong reason to assume that
every study estimates an identical effect.

The common-effect model (all studies share one $`\theta`$) is a limiting
case at $`\tau = 0`$. See [Common-Effect
Model](https://blmoran.github.io/bayesma/articles/common-effect-model.md)
for when it is preferable.

## Model specification

Let $`y_i`$ be the effect estimate and $`s_i`$ the known standard error
from study $`i = 1, \ldots, k`$. The two-level hierarchy is:

**Likelihood:**
``` math

y_i \mid \theta_i \sim \mathcal{N}(\theta_i,\; s_i^2)
```

**Random-effects distribution:**
``` math

\theta_i \mid \mu, \tau \sim \mathcal{N}(\mu,\; \tau^2)
```

Marginalising over the study-specific effects $`\theta_i`$ gives the
integrated likelihood used for inference:

``` math

y_i \mid \mu, \tau \sim \mathcal{N}(\mu,\; s_i^2 + \tau^2)
```

Studies with smaller $`s_i`$ still receive more weight, but the weight
differential narrows as $`\tau`$ grows.

## Priors

**bayesma** defaults:

``` math

\mu \sim \mathcal{N}(0,\; 1), \qquad \tau \sim \text{Half-Cauchy}(0,\; 0.5)
```

The Half-Cauchy(0, 0.5) prior on $`\tau`$ is weakly informative: it
assigns most prior mass to small heterogeneity while keeping
moderate-to-large values within the plausible range. For effect sizes on
standardised scales (log-OR, SMD), heterogeneity above $`\tau = 1`$ is
already very large.

Priors can be overridden via the `prior_mu` and `prior_tau` arguments:

``` r
#| eval: false
fit <- bayesma(
  data,
  prior_mu  = normal(0, 0.5),
  prior_tau = half_normal(0, 0.5)
)
```

See [Prior Predictive
Checks](https://blmoran.github.io/bayesma/articles/prior-predictive-checks.md)
and [Sensitivity
Analysis](https://blmoran.github.io/bayesma/articles/sensitivity-analysis.md)
for guidance.

## Fitting the model

``` r
#| eval: false
fit_re <- bayesma(data)  # Gaussian RE is the default

summary(fit_re)
```

The [`summary()`](https://rdrr.io/r/base/summary.html) output reports
the posterior median and 95% credible interval for $`\mu`$, $`\tau`$,
and $`I^2`$, along with convergence diagnostics.

## Key estimands

### Overall effect $`\mu`$

The posterior for $`\mu`$ quantifies the mean of the distribution of
true effects. This is the best estimate of what a new, exchangeable
study would find on average.

### Between-study heterogeneity $`\tau`$

$`\tau`$ is the standard deviation of the true-effect distribution.
Larger $`\tau`$ means more variation in true effects across studies. The
posterior for $`\tau`$ should be examined directly — if it has
substantial mass near zero, the data are consistent with negligible
heterogeneity.

### $`I^2`$

**bayesma** reports the posterior distribution of $`I^2`$, defined as:

``` math

I^2 = \frac{\tau^2}{\tau^2 + \tilde{s}^2}
```

where $`\tilde{s}^2`$ is the typical within-study variance. $`I^2`$ is a
relative measure: the proportion of total observed variance attributable
to between-study heterogeneity. It depends on both $`\tau`$*and* the
precision of the primary studies, so it should be interpreted alongside
$`\tau`$ rather than in isolation.

### Prediction interval

The 95% prediction interval for a new study’s true effect is

``` math

\mu \pm 1.96\,\tau
```

This is wider than the credible interval for $`\mu`$ and is reported by
[`summary()`](https://rdrr.io/r/base/summary.html) when
`prediction_interval = TRUE`.

## Assumptions

The Gaussian random-effects model assumes:

1.  **Exchangeability**: the true effects $`\theta_i`$ are drawn from
    the same distribution. This requires that the studies are
    sufficiently similar in population, intervention, comparator, and
    outcome to justify pooling.
2.  **Normal random-effects distribution**: the distribution of true
    effects is symmetric. When this is implausible — for example, when
    outlier studies are suspected — consider [Alternative RE
    Distributions](https://blmoran.github.io/bayesma/articles/random-effects-distributions.md)
    or [RE Mixture
    Models](https://blmoran.github.io/bayesma/articles/robust-mixture-models.md).
3.  **Known within-study variances**: $`s_i^2`$ is treated as fixed at
    the estimated value. This approximation deteriorates for very small
    studies.

## Heterogeneity assessment

``` r
#| eval: false
bayesma_output(fit_re, type = "heterogeneity")
```

This produces a table of $`\tau`$, $`I^2`$, and $`H^2`$ with posterior
intervals. For a visual assessment see [Funnel
Plots](https://blmoran.github.io/bayesma/articles/funnel-plots.md) and
[Posterior Predictive
Checks](https://blmoran.github.io/bayesma/articles/posterior-predictive-checks.md).

## Choosing between common-effect and random-effects

The random-effects model is the appropriate default for most
meta-analyses. Reserve the common-effect model for planned sensitivity
analyses, when exchangeability is genuinely implausible to question
(e.g., replications of the same study by the same team), or when the
number of studies is very small and $`\tau`$ cannot be estimated
reliably.

Bayesian model comparison via
[`compare_models()`](https://blmoran.github.io/bayesma/reference/compare_models.md)
provides formal evidence for or against heterogeneity:

``` r
#| eval: false
fit_ce <- bayesma(data, model = "common_effect")
compare_models(fit_ce, fit_re, labels = c("Common-effect", "Random-effects"))
```
