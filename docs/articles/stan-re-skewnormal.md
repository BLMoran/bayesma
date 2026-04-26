# Stan code — skew-normal random-effects model

## Model description

The skew-normal random-effects model allows the distribution of true
study effects to be asymmetric. This is appropriate when effects are
bounded below (or above) by a natural floor or ceiling, or when the
literature is characterised by many small effects and a long tail of
large effects in one direction.

The skew-normal distribution generalises the Gaussian with a shape
parameter $`\alpha_\text{sk}`$. Positive $`\alpha_\text{sk}`$ produces
right skew; negative produces left skew; $`\alpha_\text{sk} = 0`$
recovers the Gaussian.

## Mathematical specification

**Likelihood:**

``` math

y_i \mid \theta_i \sim \mathcal{N}(\theta_i,\, s_i^2)
```

**Random effects:**

``` math

\theta_i \sim \text{SN}(\xi,\, \omega,\, \alpha_\text{sk})
```

where $`\xi`$ is the location, $`\omega > 0`$ is the scale, and
$`\alpha_\text{sk}`$ is the shape. The mean and variance of the
skew-normal are:

``` math

\mathbb{E}[\theta_i] = \xi + \omega \cdot \delta \cdot \sqrt{2/\pi}, \quad \delta = \frac{\alpha_\text{sk}}{\sqrt{1 + \alpha_\text{sk}^2}}
```

``` math

\text{Var}(\theta_i) = \omega^2 \left(1 - \frac{2\delta^2}{\pi}\right)
```

**Priors:**

``` math

\xi \sim \mathcal{N}(0,\, 1), \qquad \omega \sim \text{Half-Cauchy}(0,\, 0.5), \qquad \alpha_\text{sk} \sim \mathcal{N}(0,\, 1)
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

## How bayesma calls this model

Selected by `model_type = "random_effect"` with
`re_dist = "skew_normal"`.

``` r
bayesma(
  data,
  model_type = "random_effect",
  re_dist    = "skew_normal",
  alpha_prior = normal(0, 1)
)
```

`b_Intercept` is set to the mean of the skew-normal distribution (not
the location $`\xi`$), so that the reported pooled effect is comparable
to the Gaussian and Student-$`t`$ models.

## Parameterisation notes

Stan’s built-in `skew_normal_lpdf` is used for efficiency. The
non-centred parameterisation for the skew-normal requires generating
truncated normal auxiliary variables; the current implementation uses
the direct skew-normal log-density instead, which is well-calibrated for
moderate $`|\alpha_\text{sk}|`$.

When $`|\alpha_\text{sk}|`$ is large (above 5), the skew-normal is
highly asymmetric and MCMC can be slow. In practice, meaningful skewness
is captured by $`|\alpha_\text{sk}| \in [0.5, 3]`$.

## Identifiability

The shape parameter $`\alpha_\text{sk}`$ is identified only when $`k`$
is sufficient to observe the distributional tail. Simulation studies
suggest $`k \geq 20`$ for reliable estimation. With smaller $`k`$, set
$`\alpha_\text{sk}`$ to a fixed value or use the Gaussian model.

## Known sampling difficulties

The posterior for $`\alpha_\text{sk}`$ and $`\omega`$ can be multimodal
when $`k`$ is small, because the data are consistent with both
$`(\omega \text{ large}, \alpha_\text{sk} \approx 0)`$ (wide symmetric)
and $`(\omega \text{ moderate}, |\alpha_\text{sk}| \text{ large})`$
(narrow skewed). Multiple chains and trace plot inspection are
essential.
