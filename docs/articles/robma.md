# Robust Bayesian Meta-Analysis (RoBMA)

## Introduction

Any single publication-bias correction method embeds strong assumptions
about the nature of the selection mechanism. PET-PEESE assumes selection
is correlated with standard error. Weight-function models assume
selection operates through p-values. Choosing between them requires
subjective judgment, and reporting only the best-fitting model ignores
model uncertainty.

Robust Bayesian Meta-Analysis (RoBMA; Maier, Bartoš & Wagenmakers, 2023)
addresses this by averaging over a family of meta-analytic models rather
than committing to one. RoBMA produces a single model-averaged posterior
for the true effect that reflects uncertainty about both the
heterogeneity structure and the publication-bias mechanism.

**bayesma** implements RoBMA via
[`robma()`](https://blmoran.github.io/bayesma/reference/robma.md).

## Model family

RoBMA fits a grid of models formed by crossing:

- **Heterogeneity**: common-effect vs. random-effects
- **Publication bias**: none vs. PET-PEESE vs. selection weight function

This yields a default grid of $`2 \times 3 = 6`$ component models. Each
model $`M_k`$ is fitted separately.

| Model   | Heterogeneity  | Bias correction  |
|---------|----------------|------------------|
| $`M_1`$ | Common effect  | None             |
| $`M_2`$ | Random effects | None             |
| $`M_3`$ | Common effect  | PET-PEESE        |
| $`M_4`$ | Random effects | PET-PEESE        |
| $`M_5`$ | Common effect  | Selection weight |
| $`M_6`$ | Random effects | Selection weight |

Custom model grids can be specified via the `bias_models` and
`heterogeneity` arguments.

## Model averaging

Each component model is assigned a prior probability $`p(M_k)`$ (equal
by default: $`1/6`$). After fitting, Bayes factors are used to compute
posterior model probabilities:

``` math

p(M_k \mid \mathbf{y}) \propto p(\mathbf{y} \mid M_k) \cdot p(M_k)
```

where $`p(\mathbf{y} \mid M_k)`$ is the marginal likelihood of model
$`k`$, estimated via bridge sampling.

The model-averaged posterior for any quantity $`\psi`$ (e.g., $`\mu`$)
is:

``` math

p(\psi \mid \mathbf{y}) = \sum_{k=1}^{K} p(M_k \mid \mathbf{y}) \cdot p(\psi \mid \mathbf{y}, M_k)
```

## Bayes factors for heterogeneity and bias

RoBMA computes two summary Bayes factors by collapsing over the model
grid:

**Bayes factor for heterogeneity** (random-effects vs. common-effect):

``` math

\text{BF}_\tau = \frac{P(\tau > 0 \mid \mathbf{y})}{P(\tau = 0 \mid \mathbf{y})}
```

**Bayes factor for publication bias**:

``` math

\text{BF}_b = \frac{P(\text{bias} \mid \mathbf{y})}{P(\text{no bias} \mid \mathbf{y})}
```

These are the primary inferential summaries: $`\text{BF}_\tau`$
quantifies the evidence that heterogeneity exists; $`\text{BF}_b`$
quantifies the evidence that publication bias is present.

## Priors

Component model priors:

``` math

\mu \sim \mathcal{N}(0,\; 1), \qquad \tau \sim \text{Half-Cauchy}(0,\; 0.5)
```

PET-PEESE coefficient:

``` math

\beta_\text{PET} \sim \mathcal{N}(0,\; 1)
```

Weight-function selection weights:

``` math

\omega_j \sim \text{Dirichlet}(\mathbf{1})
```

Priors can be modified via
[`robma_default_priors()`](https://blmoran.github.io/bayesma/reference/robma_default_priors.md)
and passed to
[`robma()`](https://blmoran.github.io/bayesma/reference/robma.md).

## Fitting RoBMA

``` r
#| eval: false
fit_robma <- robma(data)

robma_output(fit_robma)
robma_table(fit_robma)
```

`robma_output()` prints the model-averaged summary.
[`robma_table()`](https://blmoran.github.io/bayesma/reference/robma_table.md)
returns a `gt` table with posterior model probabilities, the Bayes
factors for heterogeneity and bias, and the model-averaged posterior for
$`\mu`$.

## Interpreting RoBMA output

### Model-averaged posterior for $`\mu`$

The model-averaged posterior for $`\mu`$ integrates over all component
models, weighted by their posterior model probabilities. This is the
primary bias-corrected estimate.

### Posterior model probabilities

[`robma_table()`](https://blmoran.github.io/bayesma/reference/robma_table.md)
reports the posterior model probability for each component model
$`M_k`$. Large posterior probability on models without a bias component
($`M_1`$, $`M_2`$) indicates little evidence for publication bias.

### Bayes factors

| $`\text{BF}`$ | Interpretation                                        |
|---------------|-------------------------------------------------------|
| $`< 1`$       | Evidence against the presence of heterogeneity / bias |
| $`1`$–$`3`$   | Anecdotal                                             |
| $`3`$–$`10`$  | Moderate                                              |
| $`10`$–$`30`$ | Strong                                                |
| $`> 30`$      | Very strong                                           |

## Sensitivity analysis

Bridge sampling estimates have Monte Carlo error.
[`robma_sensitivity()`](https://blmoran.github.io/bayesma/reference/robma_sensitivity.md)
reruns the bridge sampler with different random seeds and reports the
variance of the log-marginal-likelihood estimates:

``` r
#| eval: false
fit_robma <- robma(data, iter_sampling = 5000)  # more samples → more stable BF
robma_sensitivity(fit_robma)
```

Posterior model probabilities that change substantially across seeds
indicate that the MCMC chains need more samples before bridge sampling
is reliable.

## Custom model grids

``` r
#| eval: false
fit_custom <- robma(
  data,
  bias_models     = c("none", "selection_weight"),  # exclude PET-PEESE
  heterogeneity   = c("random"),                    # random-effects only
  priors_effect   = normal(0, 0.5),
  null_range      = c(-0.1, 0.1)
)
```

The `null_range` argument specifies a region of practical equivalence
around zero. Model probabilities are then computed for models that
predict the effect is in vs. outside this range.

## Limitations

- Bridge sampling requires well-mixed MCMC chains. Use
  `iter_sampling ≥ 4000` and check convergence diagnostics before
  interpreting Bayes factors.
- RoBMA is computationally intensive — fitting six component models with
  bridge sampling takes several minutes per dataset. Use
  `robma(data, parallel = TRUE)` to fit component models in parallel.
- The default model grid assumes that the bias direction is positive
  (published effects are inflated upward). For outcomes where bias
  deflates effects, specify appropriate priors on the weight functions.

## References

Maier M, Bartoš F, Wagenmakers EJ (2023). Robust Bayesian meta-analysis:
Addressing publication bias with model-averaging. *Psychological
Methods*, 28(1), 107–122.

Bartoš F, Maier M, Wagenmakers EJ, Doucouliagos H, Stanley TD (2023).
Robust Bayesian meta-analysis: Model-averaging across formulations of
publication bias. *Psychological Methods*.

## Stan Code

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

## Parameterisation

The two components share the mean $`\mu`$ but have different scales:
$`\tau`$ for typical studies and $`\tau_{\text{out}} = C \cdot \tau`$
for outlier studies. The outlier component is much wider — it
accommodates extreme effect sizes without pulling $`\mu`$ toward them.

The Beta(1, 9) prior places prior expectation on $`\pi`$ at 0.10,
reflecting the assumption that outliers are uncommon.

The `generated quantities` block computes `p_outlier[i]` — the posterior
probability that study $`i`$ belongs to the outlier component. These are
extracted and reported by `bayesma_output()`.

## Known Sampling Difficulties

The mixture parameterisation is generally well-behaved because both
components share the same mean $`\mu`$, reducing label-switching risk.
Use `adapt_delta = 0.95` and inspect per-chain trace plots for `pi_out`
and `tau`.

## How bayesma calls this model

``` r
fit_outlier <- bayesma(
  data,
  model   = "robust_outlier",
  C       = 10,
  prior_pi_out = beta(1, 9)
)

summary(fit_outlier)
```

Studies with `p_outlier > 0.5` are flagged as probable outliers in the
summary table via
`bayesma_output(fit_outlier, type = "outlier_probabilities")`.
