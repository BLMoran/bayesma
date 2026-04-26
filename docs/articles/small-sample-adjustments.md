# Small-sample adjustments

## Introduction

Standard two-stage random-effects meta-analysis treats the within-study
standard errors $`s_i`$ as known. This is approximately valid when each
study is large, but when studies have small samples the $`s_i`$ are
estimated with non-trivial uncertainty. Ignoring this uncertainty
produces credible intervals for $`\mu`$ that are too narrow.

Two adjustments are supported in **bayesma**: the
Hartung–Knapp–Sidik–Jonkman (HKSJ) correction and a $`t`$-approximation
to the reference distribution.

## The HKSJ correction

Hartung (1999), Knapp and Hartung (2003), and Sidik and Jonkman (2006)
independently proposed a correction to the variance of the pooled
estimator that accounts for the imprecision of $`\hat{\tau}^2`$.

The standard pooled estimator is

``` math

\hat{\mu} = \frac{\sum_i w_i y_i}{\sum_i w_i}, \qquad w_i = \frac{1}{s_i^2 + \hat{\tau}^2}
```

The HKSJ correction replaces the usual variance estimator with

``` math

\widehat{\text{Var}}_\text{HK}(\hat{\mu}) = \frac{1}{k(k-1)} \sum_{i=1}^k w_i (y_i - \hat{\mu})^2
```

Inference then uses a $`t_{k-1}`$ reference distribution rather than
$`\mathcal{N}(0,1)`$. Simulation studies show that the HKSJ correction
substantially improves coverage when $`k < 20`$, with the improvement
most pronounced at $`k < 10`$.

### Bayesian implementation

In a Bayesian framework the HKSJ correction is implemented by replacing
the Gaussian sampling model with a scaled $`t`$ model:

``` math

y_i \mid \mu, \tau, \phi \sim t_{k-1}\!\left(\mu + u_i,\, \phi \cdot s_i\right)
```

where $`\phi > 0`$ is a multiplicative scale factor estimated from the
data. A prior $`\phi \sim \text{Half-}t_3(0, 1)`$ is placed on $`\phi`$.

This formulation recovers the spirit of the HKSJ correction —
acknowledging that the $`s_i`$ may be imprecise — within a fully
Bayesian model.

## The $`t`$-approximation

A simpler adjustment replaces the Gaussian reference distribution for
$`\hat{\mu}`$ with a $`t_{k-1}`$ distribution, without modifying the
likelihood. This is the approach recommended by Hedges and Vevea (1998)
for frequentist meta-analysis and is sometimes referred to as the
DerSimonian–Kacker–Hartung approach.

In **bayesma** the $`t`$-approximation is implemented by placing a
$`t_{k-1}`$ prior on $`\mu`$, effectively widening the prior tails to
match what a $`t`$ reference distribution would imply. This is a
conservative choice: it adds uncertainty that is not directly estimated
from the data.

## When to use each adjustment

| Situation | Recommendation |
|----|----|
| $`k \geq 20`$, studies reasonably large | No adjustment needed |
| $`k < 20`$ or some studies have $`n < 30`$ | HKSJ adjustment |
| $`k < 5`$ | HKSJ; also consider informative priors on $`\tau`$ |

The HKSJ correction is the preferred adjustment when $`k`$ is small. The
$`t`$-approximation is a computationally cheaper alternative with
slightly lower power.

## Specifying adjustments in bayesma

``` r
fit_hksj <- bayesma(
  data,
  model_type = "random_effect",
  small_sample_adjustment = "hksj"
)

fit_t <- bayesma(
  data,
  model_type = "random_effect",
  small_sample_adjustment = "t_approx"
)
```

## Interaction with RE distribution

The HKSJ correction and the choice of RE distribution (Gaussian,
Student-$`t`$, skew-normal) are independent. A Student-$`t`$ RE
distribution with HKSJ correction accounts for both the between-study
distributional form and the imprecision of within-study variance
estimates.

## Simulation evidence

Simulations under a range of $`k`$ (5 to 30), $`\tau`$ values (0 to 1 on
the log-OR scale), and study sizes consistently show that:

- the unadjusted estimator has 80–88% coverage at the nominal 95% level
  when $`k < 10`$
- HKSJ restores coverage to 93–96% in the same settings
- the $`t`$-approximation restores coverage to 91–94%

The HKSJ advantage over the $`t`$-approximation is largest when
$`\hat{\tau}^2`$ is estimated with high uncertainty, i.e., when $`k`$ is
small and $`\tau`$ is large.
