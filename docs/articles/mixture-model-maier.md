# Mixture model (Maier)

## Introduction

The mixture model approach to publication bias adjustment (Maier et al.,
2023, *metamix*) models the observed distribution of effects as a
mixture of two components: studies reporting genuine effects and studies
reflecting publication-bias-inflated estimates.

Unlike selection weight models (which reweight existing studies) or
PET-PEESE (which adjusts via regression), the mixture model directly
represents the two-population structure hypothesised under publication
bias: a component of unbiased studies and a component of overestimates
that survived the selection filter.

## Model specification

Each study $`i`$ is assigned latent membership $`z_i \in \{0, 1\}`$:

- $`z_i = 0`$: study draws from the unbiased distribution
- $`z_i = 1`$: study draws from the biased distribution

The likelihood is

``` math

p(y_i \mid \theta, \sigma_\theta, \mu_b, \sigma_b, \pi_b) = (1 - \pi_b) \cdot \mathcal{N}(y_i \mid \theta, \sigma_\theta^2 + s_i^2) + \pi_b \cdot \mathcal{N}(y_i \mid \mu_b, \sigma_b^2 + s_i^2)
```

where:

- $`\theta`$ is the true pooled effect (unbiased component mean)
- $`\sigma_\theta`$ is the true between-study heterogeneity (unbiased
  component SD)
- $`\mu_b > \theta`$ is the mean of the biased component (constrained to
  exceed the true effect)
- $`\sigma_b`$ is the SD of the biased component
- $`\pi_b`$ is the probability that a published study belongs to the
  biased component

The constraint $`\mu_b > \theta`$ encodes the assumption that
publication bias inflates effects upward (positive direction). For
negative effects, the constraint is $`\mu_b < \theta`$.

## Priors

``` math

\theta \sim \mathcal{N}(0, 1), \qquad \sigma_\theta \sim \text{Half-Cauchy}(0, 0.5)
```

``` math

\mu_b \sim \mathcal{N}(\theta + \delta_b, \sigma_b^2), \quad \delta_b \sim \text{Half-Normal}(0, 0.5)
```

``` math

\sigma_b \sim \text{Half-Cauchy}(0, 0.5), \qquad \pi_b \sim \text{Beta}(1, 4)
```

The $`\text{Beta}(1, 4)`$ prior on $`\pi_b`$ assigns prior mass
primarily to small proportions of biased studies, encoding the
assumption that most published studies are genuine.

## Key estimands

- **$`\theta`$**: The bias-corrected pooled effect, estimated from the
  unbiased component only.
- **$`\pi_b`$**: The posterior probability that a typical published
  study is drawn from the biased component.
- **Mixture-averaged effect**: $`(1 - \pi_b)\theta + \pi_b \mu_b`$ — the
  effect that would be estimated by a naive meta-analysis ignoring the
  mixture structure.

## Fitting the mixture model

``` r
fit_mix <- bayesma(
  data,
  model_type   = "mixture_model",
  p_bias_prior = beta(1, 4)
)

summary(fit_mix)
```

## Interpreting results

The posterior for $`\pi_b`$ quantifies the evidence for a biased
subpopulation. A 95% credible interval for $`\pi_b`$ that excludes zero
provides evidence for the existence of biased studies. The posterior for
$`\theta`$ provides the bias-corrected estimate.

| $`\pi_b`$ posterior median | Interpretation |
|----|----|
| $`< 0.05`$ | Little evidence of a biased component |
| $`0.05`$–$`0.20`$ | Moderate evidence; some inflation likely |
| $`> 0.20`$ | Strong evidence; substantial proportion of biased studies |

## Comparison with selection models

The mixture model and selection weight models (Vevea-Hedges, Copas)
approach the same problem from different angles:

| Aspect | Mixture model | Selection model |
|----|----|----|
| Mechanism | Latent two-population structure | Reweighting by $`p`$-value or precision |
| Estimand | Effect in unbiased component | Effect corrected for reweighting |
| Assumption | Biased studies inflate effect magnitude | Studies selected with prob $`w(p_i)`$ |
| $`k`$ requirement | $`k \geq 10`$ | $`k \geq 5`$ |

When selection is primarily $`p`$-value based, selection models are
better theoretically motivated. When the selection mechanism is unknown
or involves factors beyond significance, the mixture model is more
agnostic.

## Prior sensitivity

The mixture model is sensitive to the prior on $`\pi_b`$. A sensitivity
analysis comparing $`\text{Beta}(1, 1)`$ (uniform) and
$`\text{Beta}(1, 9)`$ (strong prior towards no bias) is recommended.

## Stan Code

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

## Parameterisation

- `delta_b = mu_b - theta` is constrained positive, enforcing the
  assumption that biased studies overestimate the true effect. This
  eliminates one of the two label-switching modes.
- `log_mix()` marginalises over the latent component assignment $`z_i`$,
  allowing continuous sampling.
- The `Beta(1, 4)` prior on $`\pi_b`$ encodes a prior expectation that
  the majority of published studies are unbiased.

## Known Sampling Difficulties

The mixture constraint `mu_b > theta` removes one label-switching mode
but the posterior for $`\pi_b`$ and $`\sigma_b`$ can still be multimodal
when the two components are not well-separated. Use
`adapt_delta = 0.99`, inspect per-chain posteriors, and report the
sensitivity to the prior on $`\pi_b`$.

## How bayesma calls this model

``` r
bayesma(
  data,
  model_type   = "mixture_model",
  p_bias_prior = beta(1, 4)
)
```
