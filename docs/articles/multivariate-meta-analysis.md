# Multivariate meta-analysis

## Introduction

Many systematic reviews report multiple outcomes per study — for
example, both a primary efficacy endpoint and a safety endpoint, or the
same outcome measured at multiple time-points. Analysing each outcome
independently discards the within-study correlation between outcomes and
can yield inconsistent conclusions (e.g., significant benefit on the
primary endpoint but not the secondary, when the two are strongly
correlated within studies).

Multivariate meta-analysis (MVMA) models all outcomes jointly,
exploiting the within-study correlation to improve estimation efficiency
and providing coherent inference across outcomes.

## Model specification

Let $`\mathbf{y}_i = (y_{i1}, \ldots, y_{iP})^\top`$ be the vector of
$`P`$ effect estimates from study $`i`$, with known within-study
covariance matrix $`\mathbf{S}_i`$. The MVMA model is

``` math

\mathbf{y}_i \mid \boldsymbol{\theta}_i \sim \mathcal{N}_P(\boldsymbol{\theta}_i,\, \mathbf{S}_i)
```

``` math

\boldsymbol{\theta}_i \sim \mathcal{N}_P(\boldsymbol{\mu},\, \boldsymbol{\Sigma})
```

where:

- $`\boldsymbol{\mu} = (\mu_1, \ldots, \mu_P)^\top`$ is the vector of
  pooled effects
- $`\boldsymbol{\Sigma} = \text{diag}(\boldsymbol{\tau}) \cdot \Omega \cdot \text{diag}(\boldsymbol{\tau})`$
  is the between-study covariance matrix
- $`\tau_p`$ is the between-study standard deviation for outcome $`p`$
- $`\Omega`$ is the between-study correlation matrix

Priors:

``` math

\mu_p \sim \mathcal{N}(0, \sigma_p^2), \quad
\tau_p \sim \text{Half-Cauchy}(0, 0.5), \quad
\Omega \sim \text{LKJ}(\eta)
```

## Within-study covariance

The within-study covariance $`\mathbf{S}_i`$ is typically unknown and
must be approximated. **bayesma** supports:

1.  **Known correlations.** If the within-study correlation
    $`\rho_\text{within}`$ is known (e.g., from individual participant
    data or published correlation matrices), $`\mathbf{S}_i`$ is
    constructed from the marginal variances $`s_{ip}^2`$ and
    $`\rho_\text{within}`$.

2.  **Imputed correlations.** When $`\rho_\text{within}`$ is not
    reported, a plausible value (typically 0.5) is imputed. Sensitivity
    to this choice should be reported.

3.  **Riley approximation.** Riley et al. (2008) proposed marginalising
    over $`\rho_\text{within}`$ with a uniform prior, avoiding the need
    for a point estimate.

## Fitting multivariate models

``` r
fit_mv <- bayesma_mv(
  data,
  outcomes  = c("lnOR_primary", "lnOR_safety"),
  se_cols   = c("se_primary", "se_safety"),
  rho_within = 0.5
)
```

The `outcomes` argument names the columns containing the per-outcome
effect estimates; `se_cols` names the standard error columns.
`rho_within` accepts a scalar (applied uniformly) or a list of
study-specific values.

## Borrowing strength

The key benefit of MVMA is borrowing strength across outcomes. When
outcome $`p`$ is observed in only a subset of studies, the model uses
the between-study correlation to impute $`\theta_{ip}`$ for studies that
did not report it. This is only valid when the missing outcomes are
missing at random — i.e., the probability of reporting does not depend
on the unreported outcome value.

## When univariate meta-analysis is adequate

Multivariate meta-analysis adds complexity. It is most useful when:

- outcomes are strongly correlated within studies
  ($`|\rho_\text{within}| > 0.3`$)
- some outcomes are partially missing across studies
- joint inference across outcomes is required (e.g., benefit–risk
  trade-off)

When outcomes are weakly correlated and fully observed, separate
univariate analyses produce nearly identical estimates.

## Limitations

- The between-study correlation $`\Omega`$ is estimated from
  between-study variation. With few studies ($`k < 10`$) or many
  outcomes ($`P > 3`$), $`\Omega`$ is poorly identified and posterior
  estimates depend heavily on the LKJ prior.
- Robust heterogeneity models (Student-$`t`$ RE, skew-normal RE) for the
  multivariate case are computationally demanding and not currently
  supported in all model-type combinations.
