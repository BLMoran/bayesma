# Stan code — Gaussian random-effects model

## Model description

The Gaussian random-effects model is the standard Bayesian meta-analysis
model. Study-level true effects $`\theta_i`$ are drawn from a normal
distribution with mean $`\mu`$ and standard deviation $`\tau`$. Observed
effects $`y_i`$ are noisy realisations of $`\theta_i`$.

## Mathematical specification

**Likelihood:**

``` math

y_i \mid \theta_i \sim \mathcal{N}(\theta_i,\, s_i^2)
```

**Random effects:**

``` math

\theta_i \sim \mathcal{N}(\mu,\, \tau^2) \quad \Longleftrightarrow \quad \theta_i = \mu + \tau z_i, \quad z_i \sim \mathcal{N}(0, 1)
```

**Priors:**

``` math

\mu \sim \mathcal{N}(0,\, 1), \qquad \tau \sim \text{Half-Cauchy}(0,\, 0.5)
```

**Derived quantities:**

``` math

u_i = \tau z_i \quad (\text{study-level deviation}), \qquad b\_\text{Intercept} = \mu
```

## Stan code

``` stan
data {
  int<lower=1> N;
  int<lower=1> K;
  vector[N] y;
  vector<lower=0>[N] se;
  array[N] int<lower=1> study;
}

parameters {
  real mu;
  real<lower=0> tau;
  vector[K] z;
}

transformed parameters {
  vector[K] u = tau * z;
}

model {
  target += normal_lpdf(mu  | 0, 1);
  target += cauchy_lpdf(tau | 0, 0.5);
  target += std_normal_lpdf(z);

  target += normal_lpdf(y | mu + u[study], se);
}

generated quantities {
  real b_Intercept = mu;
}
```

## How bayesma calls this model

This is the default for `model_type = "random_effect"` with
`re_dist = "normal"`. The Stan data list includes the `study` index
array, mapping each row in the data to a study. For two-stage
meta-analysis, each study contributes one row; $`N = K`$. For one-stage
models with multiple arms per study, $`N > K`$.

The default prior on $`\tau`$ is Half-Cauchy(0, 0.5). This places
substantial prior mass below $`\tau = 1`$ while allowing larger values.
It can be changed with `tau_prior`:

``` r
bayesma(data, model_type = "random_effect", tau_prior = half_normal(0, 0.25))
```

## Non-centred parameterisation

The `z` parameterisation separates $`\tau`$ from the shape of the random
effects. When $`\tau`$ is near zero, the centred parameterisation (using
`u` directly) creates a funnel-shaped posterior where the $`u_i`$ must
all approach zero together as $`\tau \to 0`$. The NCP avoids this by
sampling $`z`$ on a fixed $`\mathcal{N}(0,1)`$ scale.

## Identifiability and small-$`k`$ behaviour

When $`k < 5`$, $`\tau`$ is weakly identified. The posterior for
$`\tau`$ reflects the prior more heavily. In this case:

- Report the sensitivity of conclusions to the $`\tau`$ prior.
- Consider using a more informative $`\tau`$ prior from an external
  source (e.g., Turner et al. 2015).
- The posterior mean of $`\mu`$ is still a valid estimator, but its
  credible interval may be misleading.

## Known sampling difficulties

The funnel-shaped posterior near $`\tau = 0`$ can cause divergent
transitions in the centred parameterisation. The NCP (used by default)
resolves this in most cases. If divergences persist:

- Increase `adapt_delta` to 0.99.
- Place a more informative prior on $`\tau`$ (e.g.,
  `half_normal(0, 0.5)` instead of `half_cauchy(0, 0.5)`).
