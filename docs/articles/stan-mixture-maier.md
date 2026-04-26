# Stan code — Mixture model (Maier)

## Model description

The Maier mixture model (Maier et al., 2023) represents observed effects
as arising from two populations: genuinely unbiased studies and
publication-bias-inflated studies. The proportion of biased studies and
the inflated effect mean are estimated jointly with the true effect. See
[Mixture model
(Maier)](https://blmoran.github.io/bayesma/articles/mixture-model-maier.md)
for the full statistical rationale.

## Mathematical specification

**Mixture likelihood:**

``` math

p(y_i \mid \theta, \sigma_\theta, \mu_b, \sigma_b, \pi_b) = (1 - \pi_b) \cdot \mathcal{N}\!\left(y_i;\, \theta,\, \sigma_\theta^2 + s_i^2\right) + \pi_b \cdot \mathcal{N}\!\left(y_i;\, \mu_b,\, \sigma_b^2 + s_i^2\right)
```

with the constraint $`\mu_b > \theta`$ (biased studies overestimate the
effect).

**Priors:**

``` math

\theta \sim \mathcal{N}(0,\, 1), \qquad \sigma_\theta \sim \text{Half-Cauchy}(0,\, 0.5)
```

``` math

\delta_b = \mu_b - \theta \sim \text{Half-Normal}(0,\, 0.5), \qquad \sigma_b \sim \text{Half-Cauchy}(0,\, 0.5)
```

``` math

\pi_b \sim \text{Beta}(1,\, 4)
```

## Stan code

``` stan
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
}

parameters {
  real theta;
  real<lower=0> sigma_theta;
  real<lower=0> delta_b;
  real<lower=0> sigma_b;
  real<lower=0, upper=1> pi_b;
}

transformed parameters {
  real mu_b = theta + delta_b;
}

model {
  target += normal_lpdf(theta       | 0, 1);
  target += cauchy_lpdf(sigma_theta | 0, 0.5);
  target += normal_lpdf(delta_b     | 0, 0.5);
  target += cauchy_lpdf(sigma_b     | 0, 0.5);
  target += beta_lpdf(pi_b          | 1, 4);

  for (i in 1:N) {
    target += log_mix(
      pi_b,
      normal_lpdf(y[i] | mu_b,   sqrt(square(sigma_b)     + square(se[i]))),
      normal_lpdf(y[i] | theta,  sqrt(square(sigma_theta) + square(se[i])))
    );
  }
}

generated quantities {
  real b_Intercept = theta;
  real b_mu_b      = mu_b;
  real b_pi_b      = pi_b;
}
```

## How bayesma calls this model

``` r
bayesma(
  data,
  model_type   = "mixture_model",
  p_bias_prior = beta(1, 4)
)
```

## Parameterisation notes

- `delta_b = mu_b - theta` is constrained positive, enforcing the
  assumption that biased studies overestimate the true effect. This
  eliminates one of the two label-switching modes.
- `log_mix()` marginalises over the latent component assignment $`z_i`$,
  allowing continuous sampling.
- The `Beta(1, 4)` prior on $`\pi_b`$ encodes a prior expectation that
  the majority of published studies are unbiased.

## Known sampling difficulties

The mixture constraint `mu_b > theta` removes one label-switching mode
but the posterior for $`\pi_b`$ and $`\sigma_b`$ can still be multimodal
when the two components are not well-separated. Use
`adapt_delta = 0.99`, inspect per-chain posteriors, and report the
sensitivity to the prior on $`\pi_b`$.
