# Common-effect model

## Introduction

The common-effect model is the simplest Bayesian meta-analysis model. It
assumes that every study in the synthesis estimates the same underlying
true effect $`\theta`$, and that observed differences between studies
arise entirely from sampling error.

This assumption is strong. It implies that if every study had infinite
sample size, all estimates would converge on the same value. In practice
it is most defensible when:

- studies share the same population, intervention, comparator, and
  outcome
- $`k`$ is small and between-study variation cannot be reliably
  estimated
- a sensitivity check against random-effects estimates is wanted

When heterogeneity is plausible, prefer a random-effects model with
informative priors on $`\tau`$.

## Model specification

Let $`y_i`$ denote the effect estimate from study $`i = 1, \ldots, k`$,
and $`s_i`$ its known standard error. The common-effect likelihood is

``` math

y_i \mid \theta \sim \mathcal{N}(\theta,\, s_i^2)
```

A weakly informative prior is placed on the shared effect:

``` math

\theta \sim \mathcal{N}(0,\, \sigma_\theta^2)
```

where $`\sigma_\theta`$ is set to reflect plausible effect magnitudes on
the analysis scale (log-OR, SMD, etc.). The default in **bayesma** is
$`\mathcal{N}(0, 1)`$.

The posterior is

``` math

\theta \mid \mathbf{y} \propto \mathcal{N}(0, \sigma_\theta^2) \prod_{i=1}^{k} \mathcal{N}(y_i \mid \theta, s_i^2)
```

For a conjugate Gaussian prior the posterior is also Gaussian:

``` math

\theta \mid \mathbf{y} \sim \mathcal{N}\!\left(\hat{\mu},\, \hat{\sigma}^2\right)
```

``` math

\hat{\sigma}^{-2} = \sigma_\theta^{-2} + \sum_{i=1}^{k} s_i^{-2}, \qquad
\hat{\mu} = \hat{\sigma}^2 \sum_{i=1}^{k} \frac{y_i}{s_i^2}
```

This is precision-weighted averaging: studies with smaller standard
errors contribute more to the posterior.

## The two-stage common-effect model

In two-stage meta-analysis the raw outcome data are first summarised to
$`(y_i, s_i)`$ and then combined using the Gaussian likelihood above.
This is the default in **bayesma** for continuous outcomes and log-scale
binary outcomes.

## The one-stage common-effect model

For binary outcomes, two-stage analysis discards within-study
information by treating $`s_i`$ as known. The one-stage common-effect
model replaces the normal approximation with the exact binomial
likelihood:

``` math

r_{ij} \sim \text{Binomial}(n_{ij},\, \pi_{ij}), \quad j \in \{\text{ctrl}, \text{trt}\}
```

``` math

\text{logit}(\pi_{ij}) = \gamma_i + \theta \cdot \mathbb{1}[j = \text{trt}]
```

where $`\gamma_i`$ are study-specific baseline log-odds (nuisance
parameters) and $`\theta`$ is the common log-odds ratio. Because the
$`\gamma_i`$ appear as fixed effects rather than random effects, this is
the one-stage common-effect model — all heterogeneity in baseline risk
is absorbed into the study dummies.

## Estimands

The posterior for $`\theta`$ lives on the analysis scale (log-OR,
log-RR, SMD, etc.). **bayesma** reports the posterior median and
equal-tailed 95% credible interval on this scale and, where applicable,
back-transforms to the natural scale (OR, RR, etc.).

## Prior sensitivity

Because the common-effect model has a single parameter, the posterior is
sensitive to $`\sigma_\theta`$ only when $`k`$ is very small or standard
errors are large. A sensitivity check with wider and narrower priors
(e.g., $`\mathcal{N}(0, 0.5)`$ and $`\mathcal{N}(0, 2)`$) is
recommended.

## Limitations

The common-effect model does not estimate between-study heterogeneity.
If heterogeneity exists, the pooled estimate $`\hat{\mu}`$ is still
consistent but its posterior credible interval is too narrow — it
quantifies uncertainty about the study-average effect, not the
distribution of effects across a target population of studies.

For most applied meta-analyses, a random-effects model is the preferred
starting point. The common-effect model is best viewed as a boundary
case and a diagnostic tool.
