# Stan code — small-sample adjustments

## Model description

Small-sample adjustments correct the pooled estimate’s credible interval
when the within-study standard errors $`s_i`$ are estimated rather than
known. Two adjustments are implemented: the Hartung–Knapp–Sidik–Jonkman
(HKSJ) multiplicative correction and the $`t`$-approximation. See
[Small-sample
adjustments](https://blmoran.github.io/bayesma/articles/small-sample-adjustments.md)
for the statistical rationale.

## HKSJ adjustment

### Mathematical specification

**Likelihood (HKSJ):**

``` math

y_i \mid \mu, \tau, \phi \sim t_{k-1}\!\left(\mu + u_i,\, \phi \cdot s_i\right)
```

``` math

u_i \sim \mathcal{N}(0,\, \tau^2)
```

**Priors:**

``` math

\mu \sim \mathcal{N}(0,\, 1), \qquad \tau \sim \text{Half-Cauchy}(0,\, 0.5), \qquad \phi \sim \text{Half-}t_3(0,\, 1)
```

### Stan code (HKSJ)

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
  real<lower=0> phi;
  vector[K] z;
}

transformed parameters {
  vector[K] u = tau * z;
}

model {
  target += normal_lpdf(mu  | 0, 1);
  target += cauchy_lpdf(tau | 0, 0.5);
  target += student_t_lpdf(phi | 3, 0, 1);
  target += std_normal_lpdf(z);

  for (i in 1:N) {
    target += student_t_lpdf(y[i] | K - 1, mu + u[study[i]], phi * se[i]);
  }
}

generated quantities {
  real b_Intercept = mu;
}
```

## $`t`$-approximation

### Mathematical specification

The $`t`$-approximation widens the prior on $`\mu`$ to a $`t_{k-1}`$
distribution, effectively matching the tail behaviour of the frequentist
$`t`$-based confidence interval:

``` math

\mu \sim t_{k-1}(0,\, \sigma_\mu)
```

``` math

y_i \mid \theta_i \sim \mathcal{N}(\theta_i,\, s_i^2), \quad \theta_i \sim \mathcal{N}(\mu,\, \tau^2)
```

### Stan code ($`t`$-approximation)

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
  target += student_t_lpdf(mu | K - 1, 0, 1);
  target += cauchy_lpdf(tau   | 0, 0.5);
  target += std_normal_lpdf(z);

  target += normal_lpdf(y | mu + u[study], se);
}

generated quantities {
  real b_Intercept = mu;
}
```

## How bayesma calls these models

``` r
bayesma(data, model_type = "random_effect", small_sample_adjustment = "hksj")
bayesma(data, model_type = "random_effect", small_sample_adjustment = "t_approx")
```

## Parameterisation notes

The HKSJ model estimates $`\phi`$ as an additional parameter. Values
$`\phi < 1`$ indicate less variability than expected from the $`s_i`$
alone (uncommon); $`\phi > 1`$ indicates overdispersion. The $`\phi`$
posterior serves as a diagnostic: if $`\phi`$ is concentrated well above
1, the standard errors are systematically underestimated.

## Known sampling difficulties

The HKSJ model involves a loop over observations with per-observation
$`t`$ log-densities. This is slower than the vectorised Gaussian but
typically converges without difficulty. The $`t`$-approximation has no
additional computational cost over the standard RE model.
