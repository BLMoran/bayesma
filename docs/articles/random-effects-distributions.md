# Alternative random-effects distributions

## Introduction

The standard random-effects model places a Gaussian distribution on the
study-level true effects:

``` math

\theta_i \sim \mathcal{N}(\mu, \tau^2)
```

This assumption is convenient but not always appropriate. When the
distribution of true effects is heavy-tailed, asymmetric, or bimodal —
for example, because a minority of studies target fundamentally
different populations or use qualitatively different implementations —
the Gaussian assumption can distort the pooled estimate and understate
predictive uncertainty.

**bayesma** supports three alternatives: Student-$`t`$, skew-normal, and
two-component mixture.

## Gaussian (default)

``` math

\theta_i \sim \mathcal{N}(\mu, \tau^2)
```

The Gaussian is the natural starting point. It implies symmetric
heterogeneity and exponentially light tails: very large deviations from
$`\mu`$ are treated as extremely unlikely. Use it as the default; switch
to an alternative only when there is substantive or data-driven reason.

## Student-$`t`$

``` math

\theta_i \sim t_\nu(\mu, \tau^2)
```

The Student-$`t`$ distribution adds a degrees-of-freedom parameter
$`\nu`$ that controls tail heaviness. As $`\nu \to \infty`$ it converges
to the Gaussian; at $`\nu = 3`$ the tails are substantially heavier.
Heavy tails mean that large study-level deviations from $`\mu`$ are more
plausible, reducing the influence of outlier studies on the pooled
estimate.

A prior on $`\nu`$ must be specified. **bayesma** uses

``` math

\nu \sim \text{Gamma}(2, 0.1)
```

by default, which assigns most mass to $`\nu \in (3, 30)`$ while
allowing both very heavy and near-Gaussian tails. This can be overridden
via `nu_prior`.

The Student-$`t`$ random-effects model is appropriate when:

- a small number of studies have effects far from the bulk
- the source of those deviations is unknown and possibly real (not
  artefact)
- the analyst wants to downweight outliers without excluding them

## Skew-normal

``` math

\theta_i \sim \text{SN}(\mu, \tau, \alpha)
```

The skew-normal distribution generalises the Gaussian with a shape
parameter $`\alpha`$ that controls the direction and degree of
asymmetry. When $`\alpha = 0`$ it reduces to
$`\mathcal{N}(\mu, \tau^2)`$. Positive $`\alpha`$ yields a right-skewed
distribution of effects; negative $`\alpha`$ a left-skewed distribution.

The skew-normal is appropriate when:

- there is a directional floor or ceiling effect (e.g., effects cannot
  be negative)
- a literature has a mixture of small and large positive effects but few
  negative effects
- theoretical reasons exist to expect asymmetric heterogeneity

A prior on the shape parameter is required. **bayesma** uses
$`\alpha \sim \mathcal{N}(0, 1)`$ by default.

## Two-component mixture

``` math

\theta_i \sim \pi \cdot \mathcal{N}(\mu_1, \tau_1^2) + (1 - \pi) \cdot \mathcal{N}(\mu_2, \tau_2^2)
```

The mixture model allows the distribution of true effects to be bimodal.
The mixing weight $`\pi`$ estimates the proportion of studies belonging
to each component. This is the most flexible option and the hardest to
fit reliably.

The mixture is appropriate when:

- substantive theory predicts two distinct populations of studies (e.g.,
  short vs long follow-up; low vs high dose)
- the forest plot shows a clear gap or bimodal distribution
- the analyst suspects a qualitative moderator that was not measured

**Caution**: the two-component mixture requires adequate $`k`$
(typically $`k \geq 20`$) and informative priors on $`\pi`$ and the
component parameters. With small $`k`$, the components are poorly
identified and the posterior is strongly prior-dependent. A prior
sensitivity analysis is essential.

Priors:

``` math

\pi \sim \text{Beta}(2, 2), \qquad
\mu_j \sim \mathcal{N}(0, 1), \qquad
\tau_j \sim \text{Half-Cauchy}(0, 0.5)
```

## Choosing a distribution

See [Assessment of RE
Distributions](https://blmoran.github.io/bayesma/articles/re-distribution-assessment.md)
for a model comparison workflow. As a heuristic:

| Situation                                                 | Distribution  |
|-----------------------------------------------------------|---------------|
| No strong reason to depart from Gaussian                  | Normal        |
| One or two studies with extreme effects                   | Student-$`t`$ |
| Effects are concentrated near zero with a long right tail | Skew-normal   |
| Forest plot shows two clusters of effects                 | Mixture       |

## Effect on prediction

The choice of RE distribution affects not only the pooled estimate but
also the predictive distribution for a new study. Heavier-tailed
distributions produce wider prediction intervals. This is the primary
reason to prefer heavier tails when outliers are present: the Gaussian
understates how variable a new study’s result might be.

## Stan Code

### Student-t random-effects model

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
  real<lower=2> nu;
  vector[K] z;
  vector<lower=0>[K] v;
}

transformed parameters {
  vector[K] u = tau * z ./ sqrt(v / nu);
}

model {
  target += normal_lpdf(mu   | 0, 1);
  target += cauchy_lpdf(tau  | 0, 0.5);
  target += gamma_lpdf(nu    | 2, 0.1);
  target += std_normal_lpdf(z);
  target += chi_square_lpdf(v | nu);

  target += normal_lpdf(y | mu + u[study], se);
}

generated quantities {
  real b_Intercept = mu;
}
```

### Skew-normal random-effects model

``` stan
data {
  int<lower=1> N;
  int<lower=1> K;
  vector[N] y;
  vector<lower=0>[N] se;
  array[N] int<lower=1> study;
}

parameters {
  real xi;
  real<lower=0> omega;
  real alpha_sk;
  vector[K] theta_raw;
}

transformed parameters {
  vector[K] theta;
  {
    real delta  = alpha_sk / sqrt(1 + square(alpha_sk));
    real sigma1 = omega * sqrt(1 - square(delta));
    vector[K] mu_sn = xi + omega * delta * abs(theta_raw);
    theta = mu_sn + sigma1 * theta_raw;
  }
}

model {
  target += normal_lpdf(xi       | 0, 1);
  target += cauchy_lpdf(omega    | 0, 0.5);
  target += normal_lpdf(alpha_sk | 0, 1);
  target += std_normal_lpdf(theta_raw);

  target += skew_normal_lpdf(theta | xi, omega, alpha_sk);
  target += normal_lpdf(y | theta[study], se);
}

generated quantities {
  real b_Intercept = xi + omega * (alpha_sk / sqrt(1 + square(alpha_sk))) * sqrt(2.0 / pi());
}
```

### Two-component mixture random-effects model

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

## Parameterisation

**Student-t**: The scale-mixture representation samples $`\nu`$ and the
auxiliary $`v_i`$ jointly. This is preferable to directly sampling from
a $`t`$ distribution in Stan because it avoids the non-standard $`t`$
log-density computation and produces more efficient sampling. The
constraint `real<lower=2> nu` ensures the variance of the $`t`$
distribution is finite.

**Skew-normal**: Stan’s built-in `skew_normal_lpdf` is used for
efficiency. When $`|\alpha_\text{sk}|`$ is large (above 5), the
skew-normal is highly asymmetric and MCMC can be slow. In practice,
meaningful skewness is captured by $`|\alpha_\text{sk}| \in [0.5, 3]`$.
`b_Intercept` is set to the mean of the skew-normal distribution (not
the location $`\xi`$), so that the reported pooled effect is comparable
across models.

**Mixture**: The `log_mix()` function marginalises over the discrete
component assignment $`z_i \in \{1, 2\}`$, which avoids the poor mixing
that occurs when sampling discrete parameters directly. The constraint
on component means ($`\mu_1 < 0 < \mu_2`$) eliminates label switching
but may be inappropriate if both components are expected to have the
same sign. In that case, use an ordered vector: `ordered[2] mu`.

## Identifiability

**Student-t**: The degrees-of-freedom parameter $`\nu`$ is poorly
identified when $`k < 15`$. With few studies, the data are consistent
with a wide range of $`\nu`$ values, and the posterior for $`\nu`$ is
largely prior-driven. Consider setting $`\nu`$ to a fixed value (e.g.,
$`\nu = 5`$) via a tight prior: `gamma(50, 10)`.

**Skew-normal**: The shape parameter $`\alpha_\text{sk}`$ is identified
only when $`k`$ is sufficient to observe the distributional tail.
Simulation studies suggest $`k \geq 20`$ for reliable estimation.

**Mixture**: Two-component mixture models require $`k \geq 20`$ for
reliable separation of the components. With smaller $`k`$, the mixing
weight $`\pi`$ is poorly identified and results are strongly
prior-dependent.

## Known Sampling Difficulties

**Student-t**: The joint sampling of $`\nu`$ and $`v_i`$ can be slow
when $`\nu`$ is large (near-Gaussian tail). Increasing the number of
chains or iterations helps. Persistent divergences near $`\tau = 0`$ are
handled by the non-centred parameterisation as in the Gaussian model.

**Skew-normal**: The posterior for $`\alpha_\text{sk}`$ and $`\omega`$
can be multimodal when $`k`$ is small, because the data are consistent
with both $`(\omega \text{ large}, \alpha_\text{sk} \approx 0)`$ (wide
symmetric) and
$`(\omega \text{ moderate}, |\alpha_\text{sk}| \text{ large})`$ (narrow
skewed). Multiple chains and trace plot inspection are essential.

**Mixture**: Mixture models are among the most challenging posteriors to
sample efficiently in Stan. Divergences and slow mixing are common.
Mitigations: use `adapt_delta = 0.99`, increase warmup iterations
(`iter_warmup = 2000`), inspect all chain traces individually, and run
prior predictive simulation to verify identifiability before fitting.

## How bayesma calls this model

``` r
# Student-t
bayesma(
  data,
  model_type = "random_effect",
  re_dist    = "t",
  nu_prior   = gamma(2, 0.1)
)

# Skew-normal
bayesma(
  data,
  model_type = "random_effect",
  re_dist    = "skew_normal",
  alpha_prior = normal(0, 1)
)

# Mixture
bayesma(
  data,
  model_type = "random_effect",
  re_dist    = "mixture"
)
```

The `Gamma(2, 0.1)` prior on $`\nu`$ places most mass on
$`\nu \in (5, 40)`$, allowing substantial flexibility between
near-Gaussian ($`\nu \approx 30`$) and heavy-tailed ($`\nu \approx 5`$)
behaviour.
