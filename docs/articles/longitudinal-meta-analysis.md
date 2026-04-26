# Longitudinal meta-analysis

## Introduction

Many clinical and psychological interventions produce effects that
evolve over time. When studies report outcomes at multiple follow-up
points, a longitudinal meta-analysis (LMA) models the trajectory of the
treatment effect across time, rather than selecting a single time-point
for pooling.

LMA answers questions that univariate meta-analysis cannot:

- Does the treatment effect grow, decay, or plateau over time?
- At which time-point is the effect largest?
- Is there evidence that long-term and short-term effects differ?

## Model specification

Let $`y_{it}`$ be the effect estimate from study $`i`$ at time
$`t \in \{t_1, \ldots, t_{T_i}\}`$, with known standard error
$`s_{it}`$. The longitudinal model has two levels.

**Level 1 (within-study):**

``` math

y_{it} \mid \theta_{it} \sim \mathcal{N}(\theta_{it},\, s_{it}^2)
```

**Level 2 (between-study trajectory):**

``` math

\theta_{it} = f(t;\, \boldsymbol{\mu}) + u_i + v_{it}
```

where:

- $`f(t;\, \boldsymbol{\mu})`$ is the population-level trajectory
  (parametric or non-parametric)
- $`u_i \sim \mathcal{N}(0, \tau_u^2)`$ is a study-level random
  intercept
- $`v_{it} \sim \mathcal{N}(0, \tau_v^2)`$ is a study-by-time random
  deviation

## Trajectory specifications

### Linear trajectory

``` math

f(t;\, \boldsymbol{\mu}) = \mu_0 + \mu_1 t
```

Appropriate when the effect changes at a roughly constant rate.
$`\mu_0`$ is the effect at $`t = 0`$ (baseline) and $`\mu_1`$ is the
rate of change per unit time.

### Exponential decay

``` math

f(t;\, \boldsymbol{\mu}) = \mu_\infty + (\mu_0 - \mu_\infty) e^{-\lambda t}
```

Appropriate for effects that peak near treatment initiation and decay
toward a long-run asymptote $`\mu_\infty`$. $`\lambda > 0`$ controls the
decay rate.

### Piecewise linear (spline)

``` math

f(t;\, \boldsymbol{\mu}) = \mu_0 + \mu_1 t + \sum_j \mu_{2j} (t - \kappa_j)_+
```

where $`\kappa_j`$ are knot positions and $`(x)_+ = \max(0, x)`$.
Flexible enough to capture non-monotone trajectories.

## Priors

``` math

\mu_0, \mu_1 \sim \mathcal{N}(0, 1), \quad
\tau_u, \tau_v \sim \text{Half-Cauchy}(0, 0.5)
```

Trajectory-specific parameters (e.g., $`\lambda`$, $`\mu_\infty`$)
require domain-informed priors. Exponential decay requires
$`\lambda > 0`$; a $`\text{Half-Normal}(0, 1)`$ prior is a reasonable
weakly informative choice.

## Fitting longitudinal models

``` r
fit_lma <- bayesma(
  data,
  model_type = "random_effect",
  time_col   = "weeks",
  trajectory = "linear"
)
```

The `time_col` argument names the column containing the time variable.
`trajectory` accepts `"linear"`, `"exponential"`, or `"spline"`.

## Estimands

[`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md)
returns the posterior trajectory $`f(t)`$ evaluated on a user-specified
grid via
[`bayesma_marginal()`](https://blmoran.github.io/bayesma/reference/bayesma_marginal.md).
This can be used to:

- plot the trajectory with pointwise credible bands
- identify the time of maximum effect
- compare trajectories across subgroups (via
  [`meta_reg()`](https://blmoran.github.io/bayesma/reference/meta_reg.md))

## Alignment across studies

Studies rarely measure outcomes at identical time-points. The model
handles this naturally — each study contributes $`(y_{it}, s_{it}, t)`$
triplets, and the trajectory is estimated across all studies
simultaneously. No imputation for missing time-points is needed.

## Limitations

- Longitudinal models require that time is measured comparably across
  studies (same origin, same units).
- The trajectory form must be specified. Model comparison between
  trajectory specifications is supported via
  [`compare_models()`](https://blmoran.github.io/bayesma/reference/compare_models.md).
- With $`k < 10`$ and heterogeneous time-points, both $`\tau_u`$ and
  $`\tau_v`$ are poorly identified.
