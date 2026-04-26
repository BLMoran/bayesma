# Stan code — Bias-adjusted model (Jung & Aloe 2026)

## Model description

The Jung and Aloe (2026) bias-adjusted model extends Verde (2021) by
replacing the composite risk-of-bias score with domain-specific binary
RoB indicators. Each RoB domain contributes its own bias coefficient
$`\xi_d`$, allowing the model to estimate which domains are most
strongly associated with effect inflation.

See [Bias-adjusted models (Jung & Aloe
2026)](https://blmoran.github.io/bayesma/articles/bias-corrected-jung-verde.md)
for the full statistical rationale.

## Mathematical specification

**Decomposition:**

``` math

y_i = \theta_i + b_i + \varepsilon_i, \quad \varepsilon_i \sim \mathcal{N}(0, s_i^2)
```

**True effects:**

``` math

\theta_i \sim \mathcal{N}(\mu, \tau^2)
```

**Bias component (domain-specific):**

``` math

b_i \sim \mathcal{N}\!\left(\sum_{d=1}^{D} R_{id} \cdot \xi_d,\, \sigma_b^2\right)
```

where $`R_{id} \in \{0, 1\}`$ is the RoB indicator for study $`i`$ in
domain $`d`$, and $`\xi_d`$ is the mean bias contribution from domain
$`d`$.

**Priors:**

``` math

\mu \sim \mathcal{N}(0,\, 1), \qquad \tau \sim \text{Half-Cauchy}(0,\, 0.5)
```

``` math

\xi_d \sim \mathcal{N}(0,\, 0.5), \quad d = 1, \ldots, D, \qquad \sigma_b \sim \text{Half-Normal}(0,\, 0.5)
```

## Stan code

``` stan
data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> D;
  vector[N] y;
  vector<lower=0>[N] se;
  matrix[N, D] R;              // domain-specific RoB indicators
  array[N] int<lower=1> study;
}

parameters {
  real mu;
  real<lower=0> tau;
  vector[D] xi;
  real<lower=0> sigma_b;
  vector[K] z_theta;
  vector[N] z_b;
}

transformed parameters {
  vector[K] u_theta  = tau * z_theta;
  vector[N] b_mean   = R * xi;
  vector[N] b        = b_mean + sigma_b * z_b;
}

model {
  target += normal_lpdf(mu      | 0, 1);
  target += cauchy_lpdf(tau     | 0, 0.5);
  target += normal_lpdf(xi      | 0, 0.5);
  target += normal_lpdf(sigma_b | 0, 0.5);
  target += std_normal_lpdf(z_theta);
  target += std_normal_lpdf(z_b);

  target += normal_lpdf(y | mu + u_theta[study] + b, se);
}

generated quantities {
  real b_Intercept  = mu;
  vector[D] b_xi    = xi;
}
```

## How bayesma calls this model

``` r
bayesma(
  data,
  model_type      = "bias_corrected",
  bias_source     = "jung_aloe_2026",
  rob_domain_cols = c("rob_randomisation", "rob_deviations",
                      "rob_missing", "rob_measurement", "rob_selection")
)
```

`rob_domain_cols` names the binary RoB domain columns. Each should be
coded 0 (low risk) or 1 (high risk).

## Parameterisation notes

- `R * xi` is a matrix-vector product:
  `b_mean[i] = sum_d R[i,d] * xi[d]`.
- A positive $`\xi_d`$ means high risk on domain $`d`$ is associated
  with upward bias.
- `b_xi` in `generated quantities` provides the posterior for each
  domain coefficient.

## Model comparison with Verde (2021)

The Jung & Aloe (2026) model has $`D`$ additional parameters compared to
Verde (2021). Use
[`compare_models()`](https://blmoran.github.io/bayesma/reference/compare_models.md)
to assess whether domain-specific coefficients improve predictive
accuracy:

``` r
compare_models(
  verde = fit_verde,
  jung  = fit_jung
)
```

In small meta-analyses ($`k < 20`$), the additional parameters may not
be estimable, and Verde (2021) is preferable.

## Known sampling difficulties

Same as the Verde (2021) model. With many RoB domains ($`D > 5`$)
relative to studies ($`k`$), the $`\xi_d`$ are weakly identified.
Consider a group-level prior (horseshoe or regularised horseshoe) on
$`\xi_d`$ for shrinkage towards zero.
