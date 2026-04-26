# Stan code — mixture random-effects model

## Model description

The two-component Gaussian mixture random-effects model allows the
distribution of true study effects to be bimodal. Studies are assumed to
arise from one of two populations — for example, short-term and
long-term follow-up populations, or low-dose and high-dose populations.
The mixing weight, component means, and component SDs are all estimated
from the data.

## Mathematical specification

**Likelihood:**

``` math

y_i \mid \theta_i \sim \mathcal{N}(\theta_i,\, s_i^2)
```

**Random effects (mixture):**

``` math

p(\theta_i) = \pi \cdot \mathcal{N}(\theta_i \mid \mu_1,\, \tau_1^2) + (1 - \pi) \cdot \mathcal{N}(\theta_i \mid \mu_2,\, \tau_2^2)
```

**Priors:**

``` math

\pi \sim \text{Beta}(2,\, 2), \qquad \mu_j \sim \mathcal{N}(0,\, 1), \qquad \tau_j \sim \text{Half-Cauchy}(0,\, 0.5)
```

with the constraint $`\mu_1 < \mu_2`$ imposed to resolve label
switching.

## Stan code

``` stan
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
}

parameters {
  real<upper=0> mu1;
  real<lower=0> mu2;
  real<lower=0> tau1;
  real<lower=0> tau2;
  real<lower=0, upper=1> pi_mix;
}

model {
  target += beta_lpdf(pi_mix | 2, 2);
  target += normal_lpdf(mu1  | 0, 1);
  target += normal_lpdf(mu2  | 0, 1);
  target += cauchy_lpdf(tau1 | 0, 0.5);
  target += cauchy_lpdf(tau2 | 0, 0.5);

  for (i in 1:N) {
    target += log_mix(
      pi_mix,
      normal_lpdf(y[i] | mu1, sqrt(square(tau1) + square(se[i]))),
      normal_lpdf(y[i] | mu2, sqrt(square(tau2) + square(se[i])))
    );
  }
}

generated quantities {
  real b_Intercept = pi_mix * mu1 + (1 - pi_mix) * mu2;
}
```

## How bayesma calls this model

Selected by `model_type = "random_effect"` with `re_dist = "mixture"`.

``` r
bayesma(
  data,
  model_type = "random_effect",
  re_dist    = "mixture"
)
```

The constraint $`\mu_1 < 0`$ and $`\mu_2 > 0`$ resolves label switching
by anchoring component 1 to negative effects and component 2 to positive
effects. When both components are expected to have the same sign, use a
custom model with an ordered constraint: `ordered[2] mu`.

`b_Intercept` is the mixture-weighted mean: $`(1-\pi)\mu_1 + \pi\mu_2`$.

## Parameterisation notes

The `log_mix()` function marginalises over the discrete component
assignment $`z_i \in \{1, 2\}`$, which avoids the poor mixing that
occurs when sampling discrete parameters directly.

The constraint on component means ($`\mu_1 < 0 < \mu_2`$) is a
sufficient condition to eliminate label switching but may be
inappropriate if both components are expected to have the same sign. In
that case, use an ordered vector:

``` stan
ordered[2] mu;
```

with a prior that allows both components to range over plausible values.

## Identifiability

Two-component mixture models require $`k \geq 20`$ for reliable
separation of the components. With smaller $`k`$:

- The mixing weight $`\pi`$ is poorly identified.
- The posterior may be multimodal (the two label-switching modes not
  fully eliminated by the constraint).
- Results are strongly prior-dependent.

## Known sampling difficulties

Mixture models are among the most challenging posteriors to sample
efficiently in Stan. Divergences and slow mixing are common.
Mitigations:

- Use `adapt_delta = 0.99`.
- Increase warmup iterations (`iter_warmup = 2000`).
- Inspect all chain traces individually.
- Run prior predictive simulation to verify identifiability before
  fitting.
