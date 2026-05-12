# Meta-regression

## Introduction

Meta-regression extends the random-effects model by relating
between-study variation in true effects to study-level covariates
(moderators). It is the primary tool for investigating *why* effects
vary across studies.

Common moderators include:

- participant characteristics (mean age, proportion female, disease
  severity)
- intervention characteristics (dose, duration, delivery mode)
- study characteristics (year, risk of bias, country income level)
- outcome characteristics (follow-up duration, outcome instrument)

## Model specification

Let $`x_{ij}`$ denote the $`j`$-th moderator for study $`i`$. The
meta-regression model is

``` math

y_i \mid \theta_i \sim \mathcal{N}(\theta_i,\, s_i^2)
```

``` math

\theta_i = \mu + \sum_{j=1}^{p} \beta_j x_{ij} + u_i, \qquad u_i \sim \mathcal{N}(0, \tau^2)
```

where:

- $`\mu`$ is the intercept (effect when all moderators are at their
  reference values)
- $`\beta_j`$ is the moderating effect of variable $`j`$
- $`\tau`$ is the residual between-study heterogeneity unexplained by
  the moderators

Priors:

``` math

\mu \sim \mathcal{N}(0, 1), \qquad \beta_j \sim \mathcal{N}(0, 0.5), \qquad \tau \sim \text{Half-Cauchy}(0, 0.5)
```

The default $`\beta`$ prior is weakly informative. Domain-specific
priors can be specified via `beta_priors`.

## Fitting meta-regression

``` r
fit_mr <- meta_reg(
  data,
  formula     = ~ intervention_duration + mean_age + risk_of_bias,
  model_type  = "random_effect",
  center      = TRUE,
  scale       = TRUE
)
```

The `formula` argument uses standard R formula syntax. `center = TRUE`
centres continuous moderators at their mean; `scale = TRUE` scales them
to unit standard deviation. Centering is strongly recommended to improve
sampling efficiency and interpretability of the intercept.

## Interpreting coefficients

[`coefficient_evidence()`](https://blmoran.github.io/bayesma/reference/coefficient_evidence.md)
returns the posterior median, 95% credible interval, and a Bayes factor
against the null hypothesis $`\beta_j = 0`$:

``` r
coefficient_evidence(fit_mr)
```

A $`\beta_j`$ credible interval that excludes zero indicates evidence
for moderation. The Bayes factor quantifies the evidence in favour of
$`H_1 : \beta_j \neq 0`$ relative to $`H_0 : \beta_j = 0`$.

## Visualisation

``` r
metareg_mod_plot(fit_mr, moderator = "intervention_duration")
bubble_plot(fit_mr, moderator = "mean_age")
```

[`metareg_mod_plot()`](https://blmoran.github.io/bayesma/reference/metareg_mod_plot.md)
plots the posterior regression line and credible band against each
moderator.
[`bubble_plot()`](https://blmoran.github.io/bayesma/reference/bubble_plot.md)
produces a bubble plot where bubble size represents study precision.

## Categorical moderators

Categorical moderators are dummy-coded automatically. The reference
category is the first level of the factor. Effect sizes for each
category relative to the reference are the $`\beta_j`$ posteriors for
the corresponding dummy variables.

``` r
fit_mr_cat <- meta_reg(
  data,
  formula = ~ intervention_type,
  model_type = "random_effect"
)
```

## Residual heterogeneity

A meta-regression model with moderators should be compared to the
intercept-only random-effects model:

- If $`\tau`$ decreases substantially after adding moderators, the
  moderators explain a meaningful proportion of heterogeneity.
- $`R^2_\tau = 1 - \hat{\tau}^2_\text{adjusted} / \hat{\tau}^2_\text{unadjusted}`$
  quantifies the proportion of heterogeneity explained.

If $`\tau`$ is still large after adding all measured moderators,
substantial unexplained heterogeneity remains.

## Ecological fallacy

Meta-regression estimates the association between study-level moderators
and study-level effects. This is an ecological association: it does not
identify causal patient-level effects. An intervention that works better
in studies with older mean age does not necessarily work better for
older individuals within any study.

## Common pitfalls

- **Overfitting.** With $`k < 10`$ per moderator, regression estimates
  are unreliable. As a rough rule, allow $`k/10`$ moderators.
- **Multiple testing.** Testing many moderators inflates false discovery
  rates. Report all tested moderators, not only significant ones.
- **Confounding.** Study-level moderators may be correlated. Interpret
  individual $`\beta_j`$ estimates cautiously when moderators share
  variance.
- **Missing data.** Studies with missing moderator values are dropped by
  default. Imputation or sensitivity analysis is needed when missingness
  is substantial.

## Stan Code

``` stan
data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> P;
  vector[N] y;
  vector<lower=0>[N] se;
  matrix[N, P] X;
  array[N] int<lower=1> study;
}

parameters {
  real mu;
  vector[P] beta;
  real<lower=0> tau;
  vector[K] z;
}

transformed parameters {
  vector[K] u = tau * z;
}

model {
  target += normal_lpdf(mu   | 0, 1);
  target += normal_lpdf(beta | 0, 0.5);
  target += cauchy_lpdf(tau  | 0, 0.5);
  target += std_normal_lpdf(z);

  target += normal_lpdf(y | mu + X * beta + u[study], se);
}

generated quantities {
  real b_Intercept = mu;
  vector[P] b      = beta;
}
```

## Parameterisation

`mu` is the intercept: the expected effect when all continuous
moderators are at their centred values (typically mean) and categorical
moderators are at their reference level. This interpretation depends on
centring; always use `center = TRUE` for interpretable intercepts.

`tau` is the residual between-study SD unexplained by the moderators.
When $`\tau`$ is substantially smaller than the $`\tau`$ from the
intercept-only model, the moderators account for a meaningful proportion
of heterogeneity.

## Known Sampling Difficulties

Meta-regression posteriors are well-behaved when the moderators are
uncorrelated and the sample size is adequate. Problems arise when:

- Moderators are highly collinear (VIF \> 5): the joint posterior for
  $`\boldsymbol{\beta}`$ is elongated and mixing is slow.
- The design matrix is nearly rank-deficient: Stan will warn of
  numerical issues in the matrix operations.
- $`P`$ is large relative to $`k`$: $`\tau`$ is overparameterised and
  divergences may occur.

## How bayesma calls this model

[`meta_reg()`](https://blmoran.github.io/bayesma/reference/meta_reg.md)
constructs the design matrix `X` from the user’s formula:

``` r
fit_mr <- meta_reg(
  data,
  formula = ~ duration + mean_age,
  center  = TRUE,
  scale   = TRUE
)
```

`X` has columns for each term in the formula after model matrix
expansion (dummy coding for factors, polynomial expansion for
[`poly()`](https://rdrr.io/r/stats/poly.html) terms). Centring and
scaling are applied before passing to Stan.

Per-coefficient priors can be overridden:

``` r
meta_reg(
  data,
  formula      = ~ duration + mean_age,
  beta_prior   = normal(0, 0.5),
  beta_priors  = list(duration = normal(0, 0.2), mean_age = normal(0, 1))
)
```
