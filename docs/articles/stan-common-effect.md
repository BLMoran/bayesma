# Stan code — common-effect model

## Model description

The common-effect meta-analysis model assumes all studies estimate the
same true effect $`\mu`$. Observed differences between study estimates
arise solely from sampling error. There is no between-study
heterogeneity parameter.

## Mathematical specification

**Likelihood:**

``` math

y_i \mid \mu \sim \mathcal{N}(\mu,\, s_i^2), \quad i = 1, \ldots, k
```

**Prior:**

``` math

\mu \sim \mathcal{N}(0,\, 1)
```

**Derived quantity:**

``` math

b\_\text{Intercept} = \mu
```

The posterior is analytically tractable (conjugate Gaussian) but
**bayesma** fits it via MCMC for consistency with the rest of the
pipeline.

## Stan code

``` stan
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
}

parameters {
  real mu;
}

model {
  target += normal_lpdf(mu | 0, 1);
  target += normal_lpdf(y  | mu, se);
}

generated quantities {
  real b_Intercept = mu;
}
```

## How bayesma calls this model

[`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md)
selects this model when `model_type = "common_effect"`. The data list
passed to Stan contains:

``` r
list(
  N  = nrow(data),
  y  = data$yi,       # pre-computed effect estimates
  se = data$sei       # pre-computed standard errors
)
```

Internally `bayesma_stan_data()` computes $`y_i`$ and $`s_i`$ from the
raw outcome columns (events, means, SDs, sample sizes) before passing
them to Stan.

## Parameterisation notes

- The prior $`\mu \sim \mathcal{N}(0, 1)`$ is weakly informative on the
  log-OR or SMD scale. For studies on other scales (e.g., raw mean
  differences in physical units), the prior should be widened via
  `mu_prior = normal(0, sigma)` where $`\sigma`$ reflects the plausible
  range of effects.
- The common-effect model has a single parameter. With $`k \geq 5`$, the
  likelihood dominates the prior and the posterior is largely
  data-driven.
- No random-effects structure; no $`\tau`$ parameter; no non-centred
  parameterisation needed.

## Identifiability

The model is fully identified with a single observation ($`k = 1`$). As
$`k`$ increases, the posterior for $`\mu`$ concentrates around the
precision-weighted average of the $`y_i`$.

## Known sampling difficulties

None. The posterior is unimodal and near-Gaussian. MCMC mixes rapidly
with the default settings.
