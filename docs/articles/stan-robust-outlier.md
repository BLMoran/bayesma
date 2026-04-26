# Stan code — Robust outlier mixture (RoBMA)

## Model description

The robust Bayesian model averaging (RoBMA) approach models the evidence
across multiple competing meta-analytic models: models that include or
exclude heterogeneity, and models that include or exclude publication
bias. Model weights are estimated via Bayes factors (bridge sampling)
and used to produce a model-averaged posterior for the true effect.

The component models in RoBMA are:

1.  Common effect, no publication bias
2.  Random effects, no publication bias
3.  Common effect + selection bias
4.  Random effects + selection bias

Each component is fitted separately; Bayes factors assess the evidence
for heterogeneity and for publication bias independently.

## Mathematical specification

**Model-averaged posterior:**

``` math

p(\theta \mid \mathbf{y}) = \sum_{k=1}^{K} p(M_k \mid \mathbf{y}) \cdot p(\theta \mid \mathbf{y}, M_k)
```

**Posterior model probability:**

``` math

p(M_k \mid \mathbf{y}) \propto p(\mathbf{y} \mid M_k) \cdot p(M_k)
```

where $`p(\mathbf{y} \mid M_k)`$ is the marginal likelihood of model
$`k`$ (computed via bridge sampling) and $`p(M_k)`$ is the prior model
probability.

## Stan code (component RE model with spike-and-slab heterogeneity)

``` stan
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  real<lower=0, upper=1> p_het;  // prior probability of heterogeneity
}

parameters {
  real mu;
  real<lower=0, upper=1> include_tau;  // inclusion indicator
  real<lower=0> tau_raw;
}

transformed parameters {
  real tau = include_tau > p_het ? tau_raw : 0.0;
}

model {
  target += normal_lpdf(mu      | 0, 1);
  target += bernoulli_lpdf(include_tau > p_het | p_het);
  target += cauchy_lpdf(tau_raw | 0, 0.5);

  target += normal_lpdf(y | mu, sqrt(square(se) + square(tau)));
}

generated quantities {
  real b_Intercept = mu;
}
```

Note: **bayesma** uses bridge sampling rather than spike-and-slab for
RoBMA, computing the marginal likelihood for each component model
separately and then combining via posterior model probabilities. The
above is a simplified illustration.

## How bayesma calls RoBMA

``` r
fit_robma <- robma(
  data,
  method         = "bridge",
  bias_models    = c("none", "pet_peese", "selection_weight"),
  priors_effect  = normal(0, 1),
  null_range     = c(-0.1, 0.1)
)
```

`method = "bridge"` fits each component model separately and computes
Bayes factors via bridge sampling. `method = "ss"` uses a spike-and-slab
formulation.

## Bridge sampling for Bayes factors

Bridge sampling (Meng & Wong, 1996) estimates the marginal likelihood:

``` math

p(\mathbf{y} \mid M_k) = \int p(\mathbf{y} \mid \boldsymbol{\theta}, M_k) \cdot p(\boldsymbol{\theta} \mid M_k) \cdot d\boldsymbol{\theta}
```

by fitting both the model and a reference distribution and computing a
ratio estimator. Stan’s output is passed to
[`bridgesampling::bridge_sampler()`](https://rdrr.io/pkg/bridgesampling/man/bridge_sampler.html).

## Interpreting RoBMA output

[`robma()`](https://blmoran.github.io/bayesma/reference/robma.md)
returns:

- **Component model weights**: posterior probability of each model.
- **Bayes factor for heterogeneity**:
  $`\text{BF}_\tau = P(\text{heterogeneity}) / P(\text{no heterogeneity})`$.
- **Bayes factor for bias**:
  $`\text{BF}_b = P(\text{bias}) / P(\text{no bias})`$.
- **Model-averaged posterior for $`\mu`$**: weighted average across all
  models.

## Known sampling difficulties

Bridge sampling requires that the posterior is well-explored by MCMC.
Use `iter_sampling = 5000` for reliable bridge sampling estimates. The
variance of the bridge sampling estimate can be assessed via
[`robma_sensitivity()`](https://blmoran.github.io/bayesma/reference/robma_sensitivity.md).
