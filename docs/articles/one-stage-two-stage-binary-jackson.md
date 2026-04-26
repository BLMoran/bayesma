# One-stage and two-stage models for binary outcomes

## Introduction

Meta-analyses of binary outcomes (events vs. non-events) can be
conducted using either a **two-stage** or a **one-stage** approach. The
distinction matters when studies are small, baseline event rates vary
widely, or the outcome is rare.

**Two-stage**: each study is first summarised to a log-scale effect
estimate $`y_i`$ and standard error $`s_i`$, then the meta-analysis uses
a Gaussian likelihood on those summaries. This is the default in
**bayesma** and adequate for most applications.

**One-stage**: the original binomial counts are modelled directly,
avoiding the normal approximation and handling studies with zero events
without continuity corrections. **bayesma** implements the one-stage
model following Jackson et al. (2018).

## Two-stage binary meta-analysis

### Effect size computation

For study $`i`$ with treatment arm $`(r_{1i}, n_{1i})`$ and control arm
$`(r_{0i}, n_{0i})`$, the log odds ratio is

``` math

y_i = \log\!\left(\frac{r_{1i}/(n_{1i}-r_{1i})}{r_{0i}/(n_{0i}-r_{0i})}\right)
```

with approximate variance

``` math

s_i^2 = \frac{1}{r_{1i}} + \frac{1}{n_{1i}-r_{1i}} + \frac{1}{r_{0i}} + \frac{1}{n_{0i}-r_{0i}}
```

Log-risk ratios and risk differences are supported via the
`effect_measure` argument.

### Meta-analysis model

The two-stage random-effects model applies the standard Gaussian
hierarchy to the log-scale summaries:

``` math

y_i \mid \theta_i \sim \mathcal{N}(\theta_i,\; s_i^2), \qquad \theta_i \mid \mu, \tau \sim \mathcal{N}(\mu,\; \tau^2)
```

See [Gaussian Random-Effects
Model](https://blmoran.github.io/bayesma/articles/common-effect-random-effect.md)
for full details.

### Limitations

- The normal approximation to the log-OR is poor when cell counts are
  small ($`r_i < 5`$).
- Studies with zero events require a continuity correction, which shifts
  the estimate in ways that depend on the correction chosen.
- Within-study information about baseline event rates is discarded.

## One-stage model (Jackson et al.)

The one-stage model replaces the Gaussian approximation with an exact
binomial likelihood. The formulation follows Jackson et al. (2018).

### Model specification

For study $`i \in \{1, \ldots, k\}`$ and arm
$`j \in \{0\;\text{(ctrl)},\; 1\;\text{(trt)}\}`$:

``` math

r_{ij} \sim \text{Binomial}(n_{ij},\; \pi_{ij})
```

``` math

\text{logit}(\pi_{ij}) = \gamma_i + \theta_i \cdot j
```

where $`\gamma_i`$ is the study-specific log-odds on control (baseline
log-odds) and $`\theta_i`$ is the study-specific log-odds ratio.

Treatment effects are exchangeable across studies:

``` math

\theta_i \mid \mu, \tau \sim \mathcal{N}(\mu,\; \tau^2)
```

Baseline log-odds receive weakly informative priors:

``` math

\gamma_i \sim \mathcal{N}(0,\; \sigma_\gamma^2), \quad \sigma_\gamma = 2
```

covering a wide range of baseline event rates on the logit scale.

### Priors

``` math

\mu \sim \mathcal{N}(0,\; 1), \qquad \tau \sim \text{Half-Cauchy}(0,\; 0.5)
```

These match the defaults for the two-stage model. See [Prior Predictive
Checks](https://blmoran.github.io/bayesma/articles/prior-predictive-checks.md)
for outcome-specific guidance.

### Why the one-stage model handles zero events

A study with $`r_{ij} = 0`$ contributes a valid binomial log-likelihood
term — the probability $`\text{Binomial}(0 \mid n, \pi)`$ is defined and
informative. The nuisance parameters $`\gamma_i`$ are integrated out by
MCMC, so uncertainty about baseline risk propagates into the posteriors
for $`\mu`$ and $`\tau`$.

## Fitting in bayesma

``` r
#| eval: false

# Two-stage (default for binary outcomes with pre-computed summaries)
fit_2s <- bayesma(
  data,
  outcome        = "binary",
  effect_measure = "log_or",
  model_stage    = "two_stage"
)

# One-stage (Jackson et al.) — expects raw count columns r1, n1, r0, n0
fit_1s <- bayesma(
  data,
  outcome        = "binary",
  effect_measure = "log_or",
  model_stage    = "one_stage"
)

summary(fit_1s)
```

## When to prefer the one-stage model

| Situation                               | Recommendation        |
|-----------------------------------------|-----------------------|
| Large studies, moderate event rates     | Two-stage is adequate |
| One or more studies have zero events    | One-stage preferred   |
| Small $`k`$ or sparse data              | One-stage preferred   |
| Baseline event rates vary substantially | One-stage preferred   |
| Computational speed matters             | Two-stage is faster   |

For most systematic reviews with reasonably sized studies and no
zero-event cells, the two approaches yield very similar estimates.

## Comparing estimates

``` r
#| eval: false
compare_models(fit_2s, fit_1s, labels = c("Two-stage", "One-stage"))
```

Non-trivial differences between the two approaches indicate that the
normal approximation is influencing the two-stage result. In that case,
the one-stage estimate should be reported as the primary analysis.

## References

Jackson D, Law M, Stijnen T, Viechtbauer W, White IR (2018). A
comparison of 7 random-effects models for meta-analyses that estimate
the summary odds ratio. *Statistics in Medicine*, 37(7), 1059–1085.
