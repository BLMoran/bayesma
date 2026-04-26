# Stan code — Robust outlier mixture model (Cruz et al.)

## Model description

The robust outlier mixture model (Cruz et al.) extends the standard
Gaussian random-effects model by allowing a small proportion of studies
to be outliers — studies whose true effects are so discrepant from the
bulk of the literature that they would distort estimates of $`\mu`$ and
$`\tau`$ if analysed under a single-component Gaussian.

Unlike the [two-component RE mixture
model](https://blmoran.github.io/bayesma/articles/stan-re-mixture-model.md),
which assumes two substantive subpopulations with different means, the
outlier model treats the second component as a nuisance: its role is to
absorb anomalous studies and protect inference on $`\mu`$ for the main
component.

The model is fit via
[`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md)
with `model = "robust_outlier"`.

## Mathematical specification

**Likelihood:**

``` math

y_i \mid \theta_i \sim \mathcal{N}(\theta_i,\; s_i^2)
```

**Outlier mixture prior on true effects:**

``` math

p(\theta_i) = (1 - \pi) \cdot \mathcal{N}(\theta_i \mid \mu,\; \tau^2) + \pi \cdot \mathcal{N}(\theta_i \mid \mu,\; \tau_{\text{out}}^2)
```

The two components share the mean $`\mu`$ but have different scales:
$`\tau`$ for typical studies and $`\tau_{\text{out}} \gg \tau`$ for
outlier studies. The outlier component is much wider — it accommodates
extreme effect sizes without pulling $`\mu`$ toward them.

**Marginalised likelihood:**

``` math

p(y_i \mid \mu, \tau, \tau_{\text{out}}, \pi) = (1-\pi) \cdot \mathcal{N}(y_i \mid \mu,\; \tau^2 + s_i^2) + \pi \cdot \mathcal{N}(y_i \mid \mu,\; \tau_{\text{out}}^2 + s_i^2)
```

**Outlier scale parameterisation:**

``` math

\tau_{\text{out}} = C \cdot \tau, \quad C > 1
```

The scale multiplier $`C`$ is typically fixed (default $`C = 10`$) or
assigned a prior. This ensures that the outlier component is always
wider than the main component.

**Priors:**

``` math

\mu \sim \mathcal{N}(0,\; 1), \qquad \tau \sim \text{Half-Cauchy}(0,\; 0.5)
```

``` math

\pi \sim \text{Beta}(1,\; 9)
```

The Beta(1, 9) prior places prior expectation on $`\pi`$ at 0.10,
reflecting the assumption that outliers are uncommon.

## Stan code

``` stan
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  real<lower=1> C;  // outlier scale multiplier (default 10)
}

parameters {
  real mu;
  real<lower=0> tau;
  real<lower=0, upper=1> pi_out;
}

transformed parameters {
  real tau_out = C * tau;
}

model {
  target += normal_lpdf(mu     | 0, 1);
  target += cauchy_lpdf(tau    | 0, 0.5);
  target += beta_lpdf(pi_out   | 1, 9);

  for (i in 1:N) {
    target += log_mix(
      pi_out,
      normal_lpdf(y[i] | mu, sqrt(square(tau_out) + square(se[i]))),
      normal_lpdf(y[i] | mu, sqrt(square(tau)     + square(se[i])))
    );
  }
}

generated quantities {
  real b_Intercept = mu;
  real b_tau       = tau;
  real b_pi_out    = pi_out;

  // Posterior outlier probability for each study
  vector[N] p_outlier;
  for (i in 1:N) {
    real lp_out = log(pi_out)      + normal_lpdf(y[i] | mu, sqrt(square(tau_out) + square(se[i])));
    real lp_reg = log1m(pi_out)    + normal_lpdf(y[i] | mu, sqrt(square(tau)     + square(se[i])));
    p_outlier[i] = exp(lp_out - log_sum_exp(lp_out, lp_reg));
  }
}
```

## How bayesma calls this model

``` r
#| eval: false
fit_outlier <- bayesma(
  data,
  model   = "robust_outlier",
  C       = 10,
  prior_pi_out = beta(1, 9)
)

summary(fit_outlier)
```

The `generated quantities` block computes `p_outlier[i]` — the posterior
probability that study $`i`$ belongs to the outlier component. These are
extracted and reported by
[`bayesma_output()`](https://blmoran.github.io/bayesma/reference/bayesma_output.md).

## Identifying outlier studies

``` r
#| eval: false
bayesma_output(fit_outlier, type = "outlier_probabilities")
```

Studies with `p_outlier > 0.5` are flagged as probable outliers in the
summary table. These should be investigated for data quality issues,
coding errors, or genuine moderators that explain the discrepancy.

## Key output parameters

| Parameter      | Interpretation                                     |
|----------------|----------------------------------------------------|
| `mu`           | Pooled effect for the main (non-outlier) component |
| `tau`          | Between-study heterogeneity in main component      |
| `pi_out`       | Proportion of studies in the outlier component     |
| `p_outlier[i]` | Per-study posterior outlier probability            |

## Relation to the Student-t random-effects model

A Student-t random-effects distribution (see [Alternative RE
Distributions](https://blmoran.github.io/bayesma/articles/random-effects-distributions.md))
achieves similar robustness through heavier tails rather than an
explicit mixture. The mixture parameterisation is more interpretable —
it assigns each study an outlier probability — but both approaches are
defensible when outliers are a concern.

## References

Cruz KS, et al. A robust outlier mixture model for Bayesian
meta-analysis. *Manuscript in preparation*.
