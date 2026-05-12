# Selection models for publication bias

## Introduction

Selection models address publication bias by explicitly modelling the
mechanism by which studies enter (or fail to enter) the published
literature. Rather than adjusting the effect estimate via regression,
selection models embed the publication process into the likelihood,
correcting inference about $`\mu`$ by down-weighting observations that
are over-represented due to selection.

**bayesma** implements two complementary families:

1.  The **Copas selection model** — a latent-variable model in which
    selection depends on study precision.
2.  **Vevea-Hedges weight-function models** — a step-function model in
    which selection depends on the reported p-value.

Both are available through
[`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md).

## Copas selection model

### Mechanism

Copas & Shi (2000) model each study as having been published if a latent
selection variable exceeds zero:

``` math

Z_i = \gamma_0 + \frac{\gamma_1}{SE_i} + \delta_i, \qquad \delta_i \sim \mathcal{N}(0, 1)
```

Study $`i`$ is observed if $`Z_i > 0`$. The selection parameters
$`(\gamma_0, \gamma_1)`$ control the shape of the selection mechanism:

- $`\gamma_0`$ is a baseline selection intercept (more negative → more
  severe selection overall).
- $`\gamma_1`$ governs how much easier large-precision (small-$`SE`$)
  studies are to publish. Positive $`\gamma_1`$ means precise studies
  are more likely to be published regardless of their result.

The selection variable $`Z_i`$ is correlated with the observed effect
$`y_i`$ through a joint bivariate normal:

``` math

\begin{pmatrix} y_i \\ Z_i \end{pmatrix} \sim \mathcal{N}\!\left(\begin{pmatrix} \theta_i \\ \gamma_0 + \gamma_1 / SE_i \end{pmatrix},\; \begin{pmatrix} s_i^2 & \rho\,s_i \\ \rho\,s_i & 1 \end{pmatrix}\right)
```

The correlation $`\rho`$ captures the tendency for studies with large
effects to be published (effect-dependent selection). Under the null
$`\rho = 0`$, selection depends only on precision, not on the effect
magnitude.

### Priors

``` math

\mu \sim \mathcal{N}(0, 1), \qquad \tau \sim \text{Half-Cauchy}(0, 0.5)
```

``` math

\gamma_0, \gamma_1 \sim \mathcal{N}(0, 1), \qquad \rho \sim \text{Uniform}(-1, 1)
```

Informative priors can be supplied via
[`prior_copas()`](https://blmoran.github.io/bayesma/reference/prior_bias.md):

``` r
#| eval: false
fit_copas <- bayesma(
  data,
  model        = "copas",
  prior_copas  = prior_copas(gamma0 = normal(0, 1), rho = uniform(-1, 0))
)
```

Restricting $`\rho \leq 0`$ encodes the expectation that publication
bias inflates effects upward.

## Vevea-Hedges weight-function models

### Mechanism

Vevea & Hedges (1995) model selection as a step function of the
two-tailed p-value. Studies are retained with probability $`\omega_j`$
if their p-value falls in interval $`j`$:

``` math

p(\text{study observed} \mid p_i) = \omega_j, \quad p_i \in (a_{j-1},\; a_j]
```

The likelihood is reweighted so that over-represented p-value regions
(e.g., $`p < .05`$) contribute proportionally less to the pooled
estimate.

**bayesma** implements the one-sided step function by default, with
cut-points at $`.025,\; .05,\; .10,\; .25,\; .50,\; 1.0`$:

``` math

\omega_1 \geq \omega_2 \geq \cdots \geq \omega_6 = 1
```

The final weight $`\omega_6 = 1`$ is fixed as the reference
(non-significant studies), and each $`\omega_j`$ gives the relative
probability that a study in p-value category $`j`$ is published compared
to a study with $`p > .50`$.

### Priors

The weight parameters are assigned a Dirichlet prior to ensure proper
identification:

``` math

(\omega_1, \ldots, \omega_5) \sim \text{Dirichlet}(\mathbf{1})
```

Custom weights can be specified via
[`prior_weight_function()`](https://blmoran.github.io/bayesma/reference/prior_bias.md).

### One-sided vs. two-sided step function

``` r
#| eval: false
fit_wf1 <- bayesma(
  data,
  model       = "weight_function",
  weight_type = "one_sided"   # default
)

fit_wf2 <- bayesma(
  data,
  model       = "weight_function",
  weight_type = "two_sided"
)
```

The two-sided version uses the absolute p-value and is more appropriate
when effects can be in either direction.

## Fitting and comparing selection models

``` r
#| eval: false
fit_re    <- bayesma(data)
fit_copas <- bayesma(data, model = "copas")
fit_wf    <- bayesma(data, model = "weight_function")

compare_models(fit_re, fit_copas, fit_wf,
               labels = c("Random-effects", "Copas", "Weight function"))
```

[`compare_models()`](https://blmoran.github.io/bayesma/reference/compare_models.md)
overlays the posteriors for $`\mu`$ under each model, making it easy to
see how much the publication-bias correction shifts the estimate and
whether the two selection models agree.

## Interpreting the output

For the Copas model, [`summary()`](https://rdrr.io/r/base/summary.html)
returns:

| Parameter | Interpretation                         |
|-----------|----------------------------------------|
| `mu`      | Bias-corrected pooled effect           |
| `tau`     | Between-study heterogeneity            |
| `gamma0`  | Selection intercept                    |
| `gamma1`  | Precision-dependent selection          |
| `rho`     | Effect-dependent selection correlation |

For the weight-function model:

| Parameter | Interpretation                                              |
|-----------|-------------------------------------------------------------|
| `mu`      | Bias-corrected pooled effect                                |
| `tau`     | Between-study heterogeneity                                 |
| `omega_j` | Relative publication probability for p-value category $`j`$ |

## Limitations

- Both models require sufficient studies ($`k \geq 20`$ is a practical
  minimum for the weight-function model; Copas can work with fewer but
  is sensitive to the prior on $`\rho`$).
- The Copas model assumes a specific functional form for the selection
  mechanism. If the true mechanism differs, the correction can be
  imprecise.
- Weight-function models assume that selection operates exclusively
  through p-values. Selection based on effect size direction or
  magnitude (without reference to p-values) is not captured.
- Neither model identifies publication bias if the filed studies are
  completely missing (i.e., the selection mechanism is so strong that no
  unpublished studies share observable characteristics with the
  published set).

For a model-averaged approach that includes both PET-PEESE and weight
functions as bias components, see [Robust Bayesian Meta-Analysis
(RoBMA)](https://blmoran.github.io/bayesma/articles/robma.md).

## References

Copas J, Shi JQ (2000). Meta-analysis, funnel plots and sensitivity
analysis. *Biostatistics*, 1(3), 247–262.

Vevea JL, Hedges LV (1995). A general linear model for estimating effect
size in the presence of publication bias. *Psychometrika*, 60(3),
419–435.

## Stan Code

### Vevea-Hedges weight-function model

``` stan
data {
  int<lower=1> N;
  int<lower=1> J;
  vector[N] y;
  vector<lower=0>[N] se;
  array[N] int<lower=1> bin;        // p-value bin for each study
  matrix[N, J] bin_probs;           // P(p in bin j | mu, tau) for each study
}

parameters {
  real mu;
  real<lower=0> tau;
  simplex[J] w_raw;
}

transformed parameters {
  vector[J] w = w_raw;
  w[1] = 1.0;
}

model {
  target += normal_lpdf(mu  | 0, 1);
  target += cauchy_lpdf(tau | 0, 0.5);
  for (j in 2:J) {
    target += beta_lpdf(w[j] | 1, 1);
  }

  for (i in 1:N) {
    real denom = dot_product(w, bin_probs[i]);
    target += log(w[bin[i]]) + normal_lpdf(y[i] | mu, sqrt(square(se[i]) + square(tau)))
            - log(denom);
  }
}

generated quantities {
  real b_Intercept = mu;
}
```

### Copas selection model

``` stan
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
}

parameters {
  real mu;
  real<lower=0> tau;
  real gamma0;
  real<lower=0> gamma1;
  real<lower=-1, upper=1> rho;
}

model {
  target += normal_lpdf(mu     | 0, 1);
  target += cauchy_lpdf(tau    | 0, 0.5);
  target += normal_lpdf(gamma0 | 0, 1);
  target += normal_lpdf(gamma1 | 0, 1);
  target += uniform_lpdf(rho   | -1, 1);

  for (i in 1:N) {
    real sigma_i = sqrt(square(se[i]) + square(tau));
    real z_sel   = gamma0 + gamma1 / se[i];
    real z_adj   = (z_sel + rho * (y[i] - mu) / sigma_i) / sqrt(1 - square(rho));

    target += normal_lpdf(y[i] | mu, sigma_i) + normal_lcdf(z_adj | 0, 1)
            - normal_lcdf(z_sel | 0, 1);
  }
}

generated quantities {
  real b_Intercept = mu;
}
```

## Parameterisation

**Vevea-Hedges**: The normalising denominator in the weighted likelihood
ensures the model correctly accounts for the fact that only published
studies are observed. The $`w_j`$ parameterisation treats the first bin
(most significant) as the reference with $`w_1 = 1`$. Weaker
significance bins have $`w_j \leq 1`$ by construction.

**Copas**: $`\gamma_1 > 0`$ enforces the direction constraint: larger
studies (smaller $`s_i`$) are more likely to be published. $`\rho`$
captures the relationship between effect magnitude and publication.
Positive $`\rho`$ means larger effects are more likely to be published
(a common publication bias mechanism). The normalising denominator
`normal_lcdf(z_sel | 0, 1)` corrects for the fact that only published
studies are observed.

## Known Sampling Difficulties

**Vevea-Hedges**: The normalising denominator $`\sum_j w_j P(\ldots)`$
depends on both $`\mu`$ and $`\tau`$, creating a complex likelihood
surface. With many bins and small $`k`$, the $`w_j`$ posterior is
multimodal. Increasing `adapt_delta` to 0.99 and using
`iter_warmup = 2000` is recommended.

**Copas**: The Copas likelihood is non-convex and can have multiple
local modes in the $`(\mu, \rho)`$ space. Use at least 4 chains, inspect
all traces, and run with `adapt_delta = 0.99`. The Copas model is not
identified without prior information on either $`\gamma_0`$ or $`\rho`$;
**bayesma** provides a sensitivity plot over a grid of
$`(\gamma_0, \gamma_1)`$ values via
`sensitivity_plot(fit_copas, type = "copas_grid")`.

## How bayesma calls these models

``` r
# Vevea-Hedges
bayesma(
  data,
  model_type       = "selection_weight",
  p_cutoffs        = c(0.025, 0.05, 0.10, 0.25, 0.50, 1.0),
  selection_priors = list(w2 = beta(1, 1), w3 = beta(1, 1))
)

# Copas
bayesma(
  data,
  model_type = "selection_copas",
  gamma_prior = list(gamma0 = normal(0, 1), gamma1 = half_normal(0, 1))
)
```

`bin_probs` for the Vevea-Hedges model is computed internally by
`bayesma_stan_data()` using the normal CDF evaluated at the bin
boundaries.
