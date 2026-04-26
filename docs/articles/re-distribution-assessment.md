# Assessment of random-effects distributions

## Introduction

Choosing between random-effects (RE) distributions — Gaussian,
Student-$`t`$, skew-normal, or mixture — should not be done by
inspection of the data alone. Visual assessment of a forest plot can
suggest departures from normality, but it conflates sampling variability
with true heterogeneity and is unreliable when $`k`$ is small.

**bayesma** provides two complementary strategies: formal predictive
model comparison via
[`compare_models()`](https://blmoran.github.io/bayesma/reference/compare_models.md),
and graphical posterior predictive assessment via the ECDF plot.

## Strategy 1: Leave-one-study-out model comparison

The preferred criterion for comparing RE distributions is
leave-one-study-out cross-validation (LOSO-CV). It asks: how well does
each model predict the withheld study’s effect, averaged over all
studies?

``` r
fit_normal  <- bayesma(data, re_dist = "normal")
fit_t       <- bayesma(data, re_dist = "t")
fit_skew    <- bayesma(data, re_dist = "skew_normal")

comparison <- compare_models(
  normal    = fit_normal,
  t         = fit_t,
  skew      = fit_skew
)

compare_table(comparison)
compare_plot(comparison)
```

[`compare_table()`](https://blmoran.github.io/bayesma/reference/compare_table.md)
reports the LOSO continuous ranked probability score (CRPS) for each
model, lower is better.
[`compare_plot()`](https://blmoran.github.io/bayesma/reference/compare_plot.md)
shows credible interval coverage at each nominal level.

LOSO-CV is preferable to WAIC or LOO-IC for RE distribution assessment
because:

- it is defined on the effect-size scale across both one-stage and
  two-stage models
- it directly measures predictive accuracy for new studies, which is
  what the RE distribution governs
- it is robust to the influential-observation problems that cause LOO-IC
  to fail in meta-analysis

## Strategy 2: ECDF plot

The empirical cumulative distribution function (ECDF) of the study-level
effects, overlaid with the posterior predictive CDF from each model,
shows whether the assumed RE distribution can generate data that look
like the observed effect distribution.

``` r
ecdf_model_plot(fit_normal, fit_t, fit_skew)
```

The ECDF plot does not replace predictive model comparison — it may look
reasonable for all models even when LOSO-CRPS clearly favours one. Use
it as a diagnostic to understand *why* models differ, not to make the
primary selection.

## Interpreting model comparison results

LOSO-CRPS differences of less than 0.1 on the log-OR scale (or
equivalent on other scales) are typically negligible. Prefer the simpler
Gaussian model unless there is a meaningful improvement in CRPS or a
compelling substantive reason to use a heavier-tailed distribution.

Coverage calibration is a useful secondary diagnostic. A well-calibrated
model should produce 50% prediction intervals that contain 50% of
withheld study effects, 80% intervals that contain 80%, and so on. A
model that systematically misses at all levels may be misspecified in
its distributional form.

## When $`k`$ is small

With $`k < 10`$, LOSO-CV has high variance and model comparison is
unreliable. In this situation:

- Default to Gaussian.
- Use a prior sensitivity analysis over $`\tau`$ rather than a
  distributional sensitivity analysis.
- Report that RE distribution uncertainty was not formally assessed.

## Prior sensitivity to the distribution choice

The RE distribution interacts with the prior on $`\tau`$. A
heavier-tailed distribution with the same $`\tau`$ prior will produce
wider prediction intervals. Before attributing a CRPS improvement to the
distributional form, verify that it persists across a range of $`\tau`$
priors.

## Recommendation

| $`k`$     | Strategy                                               |
|-----------|--------------------------------------------------------|
| $`< 10`$  | Gaussian; report as limitation                         |
| $`10–20`$ | Gaussian default; sensitivity check with Student-$`t`$ |
| $`> 20`$  | Formal LOSO-CV comparison; use ECDF as diagnostic      |
