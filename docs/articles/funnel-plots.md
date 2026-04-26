# Funnel plots

## Introduction

The funnel plot displays each study’s effect estimate against a measure
of its precision (typically the standard error). Under the common-effect
assumption with no systematic bias, study effects scatter symmetrically
around the pooled estimate, with larger scatter at lower precision. This
produces the characteristic inverted-funnel shape.

Asymmetry in the funnel — with studies clustered on one side, especially
at low precision — can indicate publication bias, small-study effects,
or genuine associations between study size and effect magnitude.

**bayesma** produces funnel plots via
[`funnel_plot()`](https://blmoran.github.io/bayesma/reference/funnel_plot.md),
with optional contour lines and Shi et al. (2020) pseudo-confidence
regions based on latent standard errors.

## Basic funnel plot

``` r
funnel_plot(fit)
```

By default this produces a scatter plot of $`y_i`$ (x-axis) vs $`s_i`$
(y-axis, inverted), with the y-axis running from 0 at the top to
$`\max(s_i)`$ at the bottom. The pooled effect line is overlaid.

## Contour-enhanced funnel plot

Contour lines at $`p = 0.10`$, $`0.05`$, and $`0.01`$ (two-sided) show
the statistical significance regions. Studies outside the $`p < 0.05`$
contour would be significant under the null of no effect.

``` r
funnel_plot(fit, contour = TRUE)
```

A cluster of studies just inside the $`p < 0.05`$ boundary — with few
studies outside it — is consistent with selective reporting of
significant results.

## Latent standard errors (Shi et al.)

Studies sometimes report effects on scales that mix signal and noise in
ways that make the standard error an unreliable precision proxy. Shi et
al. (2020) proposed plotting effects against the *latent* standard error
implied by the model, which separates publication-bias-related precision
from sampling variation.

``` r
funnel_plot(fit, se_type = "latent")
```

The latent SE funnel plot is particularly informative in the context of
the Bayesian Egger model.

## Subgroup colouring

When a grouping variable is available, studies can be coloured by
subgroup:

``` r
funnel_plot(fit, colour_by = "risk_of_bias_domain")
```

This is useful for investigating whether funnel asymmetry is driven by a
particular subgroup.

## Interpreting funnel asymmetry

Asymmetry has multiple possible causes:

| Cause | Pattern |
|----|----|
| Publication bias | Missing studies at low precision on one side |
| Small-study effects | Larger effects in smaller studies regardless of direction |
| Heterogeneity correlated with study size | Effects vary systematically with precision |
| Artefact (e.g., intervention fidelity) | Smaller, more controlled studies yield different effects |

A funnel plot is a starting point for investigation, not a conclusive
test. Formal asymmetry tests (Egger, [Bayesian Egger’s
test](https://blmoran.github.io/bayesma/articles/bayesian-egger-test.md))
and model-based adjustments (PET-PEESE, selection models) should follow.

## Supported effect measures

[`funnel_plot()`](https://blmoran.github.io/bayesma/reference/funnel_plot.md)
supports the following estimands, with automatic back-transformation for
display:

- Log-odds ratio (`OR`)
- Log-risk ratio (`RR`)
- Log-hazard ratio (`HR`)
- Log-incidence rate ratio (`IRR`)
- Mean difference (`MD`)
- Standardised mean difference (`SMD`)

## Worked example

``` r
fit_re <- bayesma(dat_binary, model_type = "random_effect")

funnel_plot(
  fit_re,
  contour   = TRUE,
  colour_by = "year_quartile",
  label     = FALSE
)
```
