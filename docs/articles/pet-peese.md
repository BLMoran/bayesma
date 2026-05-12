# PET-PEESE

## Introduction

PET-PEESE (Stanley & Doucouliagos, 2014) is a regression-based method
for detecting and correcting small-study effects — the tendency for
smaller, less precise studies to report larger effects, a pattern
consistent with publication bias or inflated estimates in underpowered
research.

The key insight is that under publication bias, small-sample studies
survive the publication filter only if they produce large effects. This
creates a correlation between standard error and effect size. PET and
PEESE use regression to estimate the effect that would be reported if a
study were infinitely large (i.e., $`SE \to 0`$).

**bayesma** implements a fully Bayesian version of PET-PEESE via
[`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md)
with `model = "pet_peese"`.

## PET — Precision Effect Test

PET regresses the effect estimate on its standard error:

``` math

y_i = \alpha + \beta \cdot SE_i + \varepsilon_i
```

The intercept $`\alpha`$ estimates the true effect at infinite precision
($`SE_i = 0`$). The slope $`\beta`$ captures the association between
standard error and effect size: a positive $`\beta`$ indicates that less
precise studies report larger effects.

The error term allows for residual heterogeneity:

``` math

\varepsilon_i \sim \mathcal{N}(0,\; s_i^2 + \tau^2)
```

**When to use PET**: PET is the correct specification when the true
effect is zero. Its intercept estimate is unbiased in that case. If the
true effect is non-zero, the $`SE_i`$ term in the regression absorbs
some of the real effect, biasing $`\hat{\alpha}`$ downward.

## PEESE — Precision Effect Estimate with Standard Error

PEESE replaces the linear $`SE_i`$ term with $`SE_i^2`$:

``` math

y_i = \alpha + \beta \cdot SE_i^2 + \varepsilon_i
```

The quadratic specification is better calibrated when the true effect is
non-zero: the $`SE_i^2`$ term changes more slowly for imprecise studies,
leaving more of the true effect in the intercept. PEESE is therefore
preferred as the primary bias-correction model when prior evidence
suggests a non-zero effect.

## Priors

``` math

\alpha \sim \mathcal{N}(0,\; 1), \qquad \beta \sim \mathcal{N}(0,\; 1), \qquad \tau \sim \text{Half-Cauchy}(0,\; 0.5)
```

These are the defaults. The prior on $`\beta`$ is weakly informative and
symmetric — it does not enforce the direction of small-study effects.

## Fitting PET-PEESE in bayesma

``` r
#| eval: false
fit_pet   <- bayesma(data, model = "pet")
fit_peese <- bayesma(data, model = "peese")

summary(fit_peese)
```

[`summary()`](https://rdrr.io/r/base/summary.html) reports the posterior
for $`\alpha`$ (bias-corrected effect), $`\beta`$ (small-study effect
slope), and $`\tau`$ (residual heterogeneity).

## The conditional PET-PEESE rule

Stanley & Doucouliagos (2014) proposed a two-step decision rule:

1.  Fit PET. If the PET intercept $`\alpha`$ is clearly different from
    zero (95% CI excludes 0), proceed to PEESE.
2.  Report the PEESE intercept as the bias-corrected effect estimate.

The Bayesian version of this rule uses the posterior probability that
$`\alpha_\text{PET} \neq 0`$:

``` r
#| eval: false
# Fit both
fit_pet   <- bayesma(data, model = "pet")
fit_peese <- bayesma(data, model = "peese")

# Posterior probability of non-zero effect (PET intercept)
interpret(fit_pet, parameter = "alpha")

# Report PEESE if PET provides clear evidence of non-zero effect
summary(fit_peese)
```

In **bayesma**, Bayes factors for $`\alpha`$ are computed using the
Savage-Dickey density ratio, providing a continuous measure of evidence
for vs. against a non-zero effect.

## Interpreting the output

| Parameter  | Interpretation                                                 |
|------------|----------------------------------------------------------------|
| `alpha`    | Bias-corrected pooled effect ($`SE \to 0`$)                    |
| `beta`     | Small-study effect coefficient                                 |
| `tau`      | Residual between-study heterogeneity                           |
| `BF_alpha` | Bayes factor for $`H_1: \alpha \neq 0`$ vs $`H_0: \alpha = 0`$ |

## Limitations

- PET-PEESE assumes that small-study effects are the *only* source of
  funnel asymmetry. Genuine small-study effects (e.g., smaller studies
  targeting different populations) can mimic publication bias.
- The method has low power when $`k < 20`$ and performs poorly when
  heterogeneity is large relative to the precision range of the included
  studies.
- If true heterogeneity and publication bias both exist, PET-PEESE may
  overcorrect.

For a more comprehensive approach that models multiple sources of bias
simultaneously, see [Robust Bayesian Meta-Analysis
(RoBMA)](https://blmoran.github.io/bayesma/articles/robma.md), which
averages over PET-PEESE and selection model assumptions.

## References

Stanley TD, Doucouliagos H (2014). Meta-regression approximations to
reduce publication selection bias. *Research Synthesis Methods*, 5(1),
60–78.

## Stan Code

``` stan
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> tau;
}

model {
  target += normal_lpdf(alpha | 0, 1);
  target += normal_lpdf(beta  | 0, 1);
  target += cauchy_lpdf(tau   | 0, 0.5);

  target += normal_lpdf(y | alpha + beta * square(se),
                           sqrt(square(se) + square(tau)));
}

generated quantities {
  real b_Intercept = alpha;
}
```

For the PET model, replace `square(se)` with `se` in both the mean and
(if needed) a different specification.

## Parameterisation

- `alpha` is the bias-corrected pooled effect (effect at $`s_i = 0`$).
- `beta` is the publication-bias slope: the rate at which effects grow
  as standard error increases.
- `tau` captures residual between-study heterogeneity not explained by
  the precision-effect relationship.
- `b_Intercept = alpha` is the estimand reported by `bayesma_output()`.

## Identifiability

PET-PEESE is identified when there is meaningful variation in $`s_i`$
across studies. When all studies have similar precision, the slope
$`\beta`$ is weakly identified and the intercept $`\alpha`$ is
uncertain. This is a limitation of the method, not a model specification
error.

## Known Sampling Difficulties

PET-PEESE is a linear regression model and samples efficiently. No
divergences are expected. The only potential issue is when $`\tau`$ is
near zero and the posterior becomes very flat; a slightly more
informative prior on $`\tau`$ helps.

## How bayesma calls this model

``` r
bayesma(data, model_type = "pet_peese", pet_peese_form = "peese")
bayesma(data, model_type = "pet_peese", pet_peese_form = "pet")
```

The default is PEESE, which is preferred when the meta-analytic effect
is expected to be non-zero. PET is preferred under the null or when
testing whether there is any effect after adjusting for publication
bias.
