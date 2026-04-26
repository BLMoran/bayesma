# Bayesian Egger’s test

## Introduction

Egger’s regression (Egger et al., 1997) detects funnel plot asymmetry by
regressing standardised effect estimates on precision. In a symmetric
funnel — the expected pattern under no publication bias — the regression
passes through the origin. A non-zero intercept indicates that
small-precision (small-sample) studies yield systematically different
effects than large-precision (large-sample) studies, a pattern
consistent with publication bias or small-study effects.

**bayesma** implements a fully Bayesian version of Egger’s regression
via [`egger()`](https://blmoran.github.io/bayesma/reference/egger.md).
The Bayesian formulation provides:

- posterior distributions for the Egger intercept and slope, not just
  p-values
- flexible modelling of residual heterogeneity (multiplicative or
  additive)
- coherent quantification of evidence for asymmetry via Bayes factors

## Model specification

Following Shi et al. (2020), the Bayesian Egger model is

``` math

\frac{y_i}{s_i} = \alpha \cdot s_i^{-1} + \beta + \varepsilon_i
```

where $`y_i / s_i`$ is the standardised effect (z-score) and $`s_i`$ is
the standard error. Rearranging:

``` math

y_i = \alpha + \beta \cdot s_i + \varepsilon_i
```

Here $`\alpha`$ is the effect at infinite precision (no publication
bias) and $`\beta`$ is the Egger regression coefficient: the rate at
which the effect changes as a function of standard error.

### Residual heterogeneity

Two models for the error term $`\varepsilon_i`$:

**Multiplicative heterogeneity** (Egger’s original formulation):
``` math

\varepsilon_i \sim \mathcal{N}(0, \kappa^2 s_i^2)
```
where $`\kappa`$ is a multiplicative overdispersion factor.
$`\kappa = 1`$ corresponds to the standard fixed-effect Egger model.

**Additive heterogeneity**:
``` math

\varepsilon_i \sim \mathcal{N}(0, s_i^2 + \tau^2)
```
where $`\tau`$ is the between-study heterogeneity, identical to the
standard random-effects parameterisation.

Priors:

``` math

\alpha \sim \mathcal{N}(0, 1), \qquad \beta \sim \mathcal{N}(0, 1), \qquad \kappa \sim \text{Half-Normal}(0, 1), \qquad \tau \sim \text{Half-Cauchy}(0, 0.5)
```

## Fitting Egger’s test

``` r
egger_fit <- egger(
  data,
  heterogeneity = "additive"
)

summary(egger_fit)
egger_plot(egger_fit)
```

`heterogeneity = "multiplicative"` reproduces the original Egger model.
`"additive"` is preferable when the goal is to estimate the
publication-bias adjusted effect $`\alpha`$ while allowing for genuine
heterogeneity.

## Interpreting the output

[`summary()`](https://rdrr.io/r/base/summary.html) returns:

| Parameter | Interpretation                                               |
|-----------|--------------------------------------------------------------|
| `alpha`   | Pooled effect adjusted for publication bias (intercept)      |
| `beta`    | Egger slope: change in effect per unit SE                    |
| `tau`     | Residual heterogeneity (additive model)                      |
| `kappa`   | Overdispersion (multiplicative model)                        |
| `BF_beta` | Bayes factor for $`H_1: \beta \neq 0`$ vs $`H_0: \beta = 0`$ |

A large $`|\beta|`$ with posterior mass away from zero indicates funnel
asymmetry. The Bayes factor quantifies evidence on a continuous scale:

| $`\text{BF}_{10}`$ | Interpretation             |
|--------------------|----------------------------|
| $`< 1`$            | Evidence against asymmetry |
| $`1`$–$`3`$        | Anecdotal                  |
| $`3`$–$`10`$       | Moderate                   |
| $`10`$–$`30`$      | Strong                     |
| $`> 30`$           | Very strong                |

## Egger plot

[`egger_plot()`](https://blmoran.github.io/bayesma/reference/egger_plot.md)
produces a scatter plot of standardised effects vs standard error with
the fitted Egger regression line and its 95% credible band.

## Limitations

Egger’s test has low power when $`k < 10`$. A non-significant result
does not rule out publication bias; it only indicates that the data are
insufficient to detect asymmetry at the measured precision.

Funnel asymmetry can also arise from causes other than publication bias:
genuine small-study effects, heterogeneity correlated with study size,
or artefacts in effect size computation. A positive Egger test should
prompt further investigation (selection models, PET-PEESE) rather than a
reflexive conclusion of publication bias.

For a visual complement to the Egger test, see [Funnel
Plots](https://blmoran.github.io/bayesma/articles/funnel-plots.md).
