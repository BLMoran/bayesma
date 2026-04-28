# Stan code — PET-PEESE

## Model description

PET-PEESE (Stanley & Doucouliagos, 2014) corrects for publication bias
by regressing effects on their standard errors (PET) or squared standard
errors (PEESE). The intercept estimates the effect at infinite precision
— the unbiased estimate.

## Mathematical specification

**PET likelihood:**

``` math

y_i \mid \alpha, \beta \sim \mathcal{N}(\alpha + \beta \cdot s_i,\, s_i^2 + \tau^2)
```

**PEESE likelihood:**

``` math

y_i \mid \alpha, \beta \sim \mathcal{N}(\alpha + \beta \cdot s_i^2,\, s_i^2 + \tau^2)
```

**Priors:**

``` math

\alpha \sim \mathcal{N}(0,\, 1), \qquad \beta \sim \mathcal{N}(0,\, 1), \qquad \tau \sim \text{Half-Cauchy}(0,\, 0.5)
```

## Stan code (PEESE)

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

## How bayesma calls this model

``` r
bayesma(data, model_type = "pet_peese", pet_peese_form = "peese")
bayesma(data, model_type = "pet_peese", pet_peese_form = "pet")
```

The default is PEESE, which is preferred when the meta-analytic effect
is expected to be non-zero. PET is preferred under the null or when
testing whether there is any effect after adjusting for publication
bias.

## Parameterisation notes

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

## Known sampling difficulties

PET-PEESE is a linear regression model and samples efficiently. No
divergences are expected. The only potential issue is when $`\tau`$ is
near zero and the posterior becomes very flat; a slightly more
informative prior on $`\tau`$ helps.
