# Model comparison (WAIC / LOO)

## Introduction

Formal model comparison in **bayesma** uses predictive accuracy
criteria. This vignette covers the widely applicable information
criterion (WAIC), LOO information criterion (LOO-IC), and
leave-one-study-out cross-validation (LOSO-CV), explaining when to use
each and how to interpret the output.

For the comparison workflow and function reference, see [Model
Comparison &
Diagnostics](https://blmoran.github.io/bayesma/articles/model-comparison-diagnostics.md).

## WAIC

WAIC (Watanabe, 2010) approximates the expected log predictive density
for new data:

``` math

\text{WAIC} = -2\left(\sum_{i=1}^{n} \log \overline{p}(y_i \mid \boldsymbol{\theta}) - \sum_{i=1}^{n} \text{Var}_{\text{post}}[\log p(y_i \mid \boldsymbol{\theta})]\right)
```

The first term is the in-sample log predictive density; the second term
penalises model complexity (effective number of parameters). Lower WAIC
is better.

WAIC is fully Bayesian, using the entire posterior distribution rather
than a point estimate. It can be computed directly from MCMC samples
without refitting.

``` r
fit$waic()
```

## LOO-IC

LOO-IC approximates leave-one-observation-out cross-validation using
Pareto-smoothed importance sampling (PSIS-LOO; Vehtari et al., 2017):

``` math

\text{LOO-IC} = -2 \sum_{i=1}^{n} \log p(y_i \mid \mathbf{y}_{-i})
```

PSIS-LOO is more reliable than WAIC for heavy-tailed posteriors and
provides per-observation diagnostics (Pareto $`\hat{k}`$).

``` r
fit$loo()
```

Pareto $`\hat{k}`$ values: - $`\hat{k} < 0.5`$: LOO reliable -
$`0.5 \leq \hat{k} < 0.7`$: LOO somewhat reliable -
$`\hat{k} \geq 0.7`$: LOO unreliable for this observation; use LOSO-CV

## LOSO-CV

LOSO-CV (leave-one-study-out) refits the model $`k`$ times, each time
withholding one study and predicting it from the remaining $`k-1`$
studies. It is the gold standard for meta-analytic model comparison
because:

- It is exact (no approximation).
- It is defined on the effect-size scale, enabling comparison across
  one-stage and two-stage models.
- It naturally handles influential studies (high-$`\hat{k}`$
  observations in LOO).

``` r
compare_models(model1, model2, criterion = "loso")
```

LOSO-CV is computationally expensive ($`k`$ additional fits per model).
Use it for the primary model comparison after candidates have been
screened with LOO-IC.

## Comparing models

``` r
comp <- compare_models(
  "Common effect" = fit_ce,
  "RE Gaussian"   = fit_re,
  "RE Student-t"  = fit_t,
  criterion       = "loso"
)

compare_table(comp)
compare_plot(comp)
```

### Interpreting compare_table()

| Column                 | Meaning                                     |
|------------------------|---------------------------------------------|
| `loso_crps`            | Mean LOSO-CRPS (lower = better)             |
| `delta_crps`           | Difference from best model                  |
| `se_delta`             | SE of the CRPS difference                   |
| `coverage_50/80/90/95` | Empirical PI coverage at each nominal level |

A model is preferred if its `loso_crps` is lowest AND its coverage is
well-calibrated (close to the nominal levels). A model with lower CRPS
but miscalibrated coverage should be investigated further.

### Comparing non-nested models

WAIC and LOO-IC differences between non-nested models (e.g.,
random-effects Gaussian vs random-effects Student-$`t`$) can be tested
using the standard error of the difference:

``` math

z = \frac{\Delta\text{LOO-IC}}{2 \cdot \hat{\sigma}_{\Delta}}
```

A large $`|z|`$ indicates the difference is unlikely to be sampling
noise. However, the null hypothesis (the two models are equally good)
should be interpreted cautiously: a small CRPS difference may be
practically negligible even if statistically significant.

## When model comparison is not appropriate

- **Nested models where the simpler model is the null.** Use Bayes
  factors (via
  [`robma()`](https://blmoran.github.io/bayesma/reference/robma.md))
  instead of LOO-IC for testing $`H_0: \tau = 0`$ vs $`H_1: \tau > 0`$.
- **Non-comparable likelihoods.** LOO-IC on different likelihoods (e.g.,
  log-OR vs log-RR) are not comparable; use LOSO-CV on a common effect
  scale.
- **Bias-corrected vs uncorrected models.** These models are answering
  different questions; model comparison is misleading. Use sensitivity
  analysis instead.
