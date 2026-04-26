# Bias-adjusted models (Jung & Aloe)

## Introduction

Publication bias and selection models address systematic inflation
driven by the publication process. A different but equally important
problem is **within-study bias**: some studies may report inflated
effects due to methodological limitations (low risk-of-bias items,
selective reporting, inadequate blinding) rather than because they were
selectively published.

Jung & Aloe (2026) extend the bias-adjusted framework of Verde (2021) to
handle this problem. The model assigns each study a latent bias
indicator and estimates how much the bias inflates the reported effect,
while still recovering a valid estimate of the true underlying effect
from the unbiased component.

**bayesma** implements this model via
[`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md)
with `model = "jung_aloe"`.

## Model specification

### Observation model

Each study $`i`$ has a latent bias indicator $`I_i \in \{0, 1\}`$, where
$`I_i = 1`$ means the study is biased. The observed effect is:

``` math

y_i \mid I_i \sim \mathcal{N}\!\left(\theta_i^B,\; s_i^2\right)
```

where the apparent effect under bias is

``` math

\theta_i^B = \theta_i + I_i \cdot \beta_i
```

- $`\theta_i`$ is the true study effect (what the study would have
  reported under ideal conditions).
- $`\beta_i`$ is the bias magnitude for study $`i`$ (zero for unbiased
  studies, positive for upwardly biased ones).
- $`I_i`$ is the binary bias indicator.

### Random-effects structure

True effects are exchangeable across studies:

``` math

\theta_i \mid \mu, \tau \sim \mathcal{N}(\mu,\; \tau^2)
```

Bias magnitudes are also treated as random:

``` math

\beta_i \mid \mu_\beta, \tau_\beta \sim \mathcal{N}(\mu_\beta,\; \tau_\beta^2), \quad \mu_\beta \geq 0
```

The constraint $`\mu_\beta \geq 0`$ encodes the assumption that bias
inflates rather than deflates effects.

Bias prevalence:

``` math

I_i \mid \phi \sim \text{Bernoulli}(\phi), \qquad \phi \sim \text{Beta}(1,\; 4)
```

The Beta(1, 4) prior places most prior mass on low bias prevalence,
reflecting the assumption that most studies in a well-conducted
systematic review are not severely biased.

### Priors

``` math

\mu \sim \mathcal{N}(0,\; 1), \qquad \tau \sim \text{Half-Cauchy}(0,\; 0.5)
```

``` math

\mu_\beta \sim \text{Half-Normal}(0,\; 0.5), \qquad \tau_\beta \sim \text{Half-Cauchy}(0,\; 0.5)
```

## Incorporating risk-of-bias information

The key extension in Jung & Aloe (2026) is that the bias indicator
probability $`\phi_i`$ can be informed by risk-of-bias assessments
(e.g., from the Cochrane RoB 2 tool or INSPECT-SR). High-risk-of-bias
studies receive a higher prior probability of being biased:

``` math

\text{logit}(\phi_i) = \alpha + \gamma \cdot \text{RoB}_i
```

where $`\text{RoB}_i`$ is a numeric risk-of-bias score for study $`i`$.

``` r
#| eval: false
fit_jung <- bayesma(
  data,
  model      = "jung_aloe",
  rob_scores = data$rob_score   # numeric RoB score per study
)
```

Without `rob_scores`, all studies share the same $`\phi`$ (Verde’s
original formulation).

## Fitting the model

``` r
#| eval: false
# Without RoB scores (Verde 2021 version)
fit_verde <- bayesma(data, model = "jung_aloe")

# With risk-of-bias scores (Jung & Aloe 2026)
fit_jung <- bayesma(
  data,
  model      = "jung_aloe",
  rob_scores = data$rob_score
)

summary(fit_jung)
```

## Key estimands

| Parameter | Interpretation |
|----|----|
| `mu` | True pooled effect (bias-corrected) |
| `tau` | Between-study heterogeneity in true effects |
| `phi` | Posterior probability of bias prevalence |
| `mu_beta` | Mean bias magnitude across biased studies |
| `tau_beta` | Between-study variation in bias magnitudes |
| `gamma` | Log-odds increase in bias probability per unit RoB (if RoB scores supplied) |

The posterior for `mu` is the primary quantity of interest: it estimates
the pooled true effect after removing the contribution of study-level
bias.

## Comparing with the uncorrected model

``` r
#| eval: false
fit_re   <- bayesma(data)
fit_jung <- bayesma(data, model = "jung_aloe", rob_scores = data$rob_score)

compare_models(fit_re, fit_jung, labels = c("Random-effects", "Jung & Aloe"))
```

A downward shift in `mu` from the standard RE model to the bias-adjusted
model indicates that within-study bias is inflating the naive pooled
estimate.

## Integration with INSPECT-SR

**bayesma** integrates with INSPECT-SR assessments. Risk-of-bias scores
from
[`inspect_sr()`](https://blmoran.github.io/bayesma/reference/inspect_sr.md)
can be passed directly:

``` r
#| eval: false
rob <- inspect_sr(data)
fit_jung <- bayesma(data, model = "jung_aloe", rob_scores = rob$rob_score)
```

See [INSPECT-SR
Assessment](https://blmoran.github.io/bayesma/articles/inspect_sr.md)
for details on computing `rob_score`.

## Limitations

- The model requires informative priors on $`\phi`$ (or per-study RoB
  scores) to identify the bias component. Without any prior information
  about bias prevalence, the bias parameters are weakly identified.
- The assumption $`\beta_i \geq 0`$ (upward bias) should be checked
  against subject-matter knowledge. Downward bias (e.g., measurement
  attenuation) requires setting $`\mu_\beta \leq 0`$.
- The model assumes that the bias process is independent of the true
  effect. This may not hold if more biased studies are conducted
  precisely when true effects are small.

## References

Verde PE (2021). A bias-corrected meta-analysis model for combining
studies of different types and quality. *Biometrical Journal*, 63(2),
406–422.

Jung Y, Aloe AM (2026). Bayesian bias-adjusted meta-analysis with
risk-of-bias information. *Research Synthesis Methods*.
