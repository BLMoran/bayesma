# Prior Predictive Checks

## What is a prior predictive check?

A prior predictive check simulates data from the model using only the
prior distributions — before any observed data are incorporated. Because
the likelihood is bypassed entirely, the draws reflect what the model
considers plausible *a priori*, under nothing but the priors and the
assumed data-generating process.

The resulting distribution, $`p(y_\text{rep})`$, answers the question:
*if I knew nothing beyond my prior beliefs and this model structure,
what data would I expect to see?*

## What prior predictive checks are for

Prior predictive checks have two legitimate uses.

**Sanity-checking on the observable scale.** Priors are typically
specified on a transformed scale (log-OR, log-RR, Cohen’s $`d`$). It is
easy to choose a prior that looks reasonable in those units but implies
absurd data — e.g. a `normal(0, 10)` prior on log-OR implies that
studies routinely report ORs in the thousands. The prior predictive
translates the prior back into the scale of the actual data, making such
problems immediately visible.

**Encoding external knowledge.** When domain expertise or independent
reference data inform the prior — for example, empirical estimates of
between-study heterogeneity from a previous systematic review — the
prior predictive lets you verify the prior is consistent with that
knowledge before analysis begins.

## What prior predictive checks are not for

It is tempting to adjust priors iteratively until $`y_\text{rep}`$
resembles the observed data. **This is circular.** Doing so smuggles the
analysis dataset into the prior, so both the prior and the likelihood
carry information from the same data. The resulting posterior will be
overconfident: it appears to be updated by the data, but the data have
already shaped the prior.

Priors must be specified from sources that are genuinely independent of
the analysis dataset:

- Domain expertise and subject-matter knowledge
- Empirical estimates from previous meta-analyses or independent cohorts
- General plausibility constraints (effect sizes cannot be infinite;
  event rates must lie in $`[0, 1]`$)

The prior predictive check is then used to *verify* the chosen prior is
coherent — not to *select* it by matching the observed data.

## Running a prior predictive check

Pass `sample_prior = TRUE` to
[`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md).
The model is compiled and sampled in the usual way, but the likelihood
contribution is excluded from the model block so the MCMC explores the
prior distribution. The `data` argument is still required: arm sizes and
standard errors are used to keep the simulation on the same observable
scale as the real data.

``` r
prior_fit <- bayesma(
  data       = my_data,
  studyvar   = "study",
  event_ctrl = "events_c",
  event_int  = "events_i",
  n_ctrl     = "n_c",
  n_int      = "n_i",
  likelihood = "binomial",
  model_type = "random_effect",
  mu_prior   = normal(0, 1),
  tau_prior  = half_cauchy(0, 0.5),
  sample_prior = TRUE
)

pp_check(prior_fit)
```

The plot overlays the observed data $`y`$ (dark line) with draws from
$`p(y_\text{rep})`$ (light lines). The observed data are shown for
reference — to orient you on the scale of the data — not as a
calibration target.

## What to look for

Evaluate the prior predictive against *domain knowledge*, not against
the observed data.

**Effect scale.** On the log-OR scale, values of $`|\mu| > 2`$
correspond to ORs above ~7 or below ~0.14. If the prior predictive
assigns substantial mass to effects of this magnitude, ask whether that
is defensible for the clinical context. Most pharmacological
interventions fall in the range $`|\text{log-OR}|
< 1.5`$; a prior placing 30% of its mass outside this range is very
vague.

**Heterogeneity scale.** A $`\tau > 1`$ on the log-OR scale implies
studies differ by more than 2 log-ORs — unusual in most literatures.
Check whether the prior predictive for $`\tau`$ is consistent with the
degree of between-study variation you would expect from independent
knowledge of the intervention area.

**Study-level effects.** The prior predictive for individual study
effects $`\theta_i`$ should cover a plausible range. Very heavy tails
(studies implying ORs \> 100) suggest the prior is too diffuse.

## Specifying priors from external knowledge

Empirically-derived priors for $`\tau`$ are available from Turner et
al. (2015) for pharmacological and non-pharmacological interventions on
binary outcomes, and from Rhodes et al. (2015) for psychological
interventions on continuous outcomes. These are calibrated from large
collections of meta-analyses and provide a principled starting point
that is independent of any single analysis:

``` r
# Turner et al. (2015) median tau for pharmacological treatments (binary)
bayesma(
  data,
  model_type = "random_effect",
  mu_prior   = normal(0, 1),
  tau_prior  = half_normal(0, 0.32)
)
```

Using externally-calibrated priors removes the temptation to tune
against the analysis data and grounds the prior in the broader evidence
base for the intervention class.

## Comparing prior and posterior

Once the model has been fitted to the data,
[`ecdf_prior_plot()`](https://blmoran.github.io/bayesma/reference/ecdf_prior_plot.md)
overlays the prior predictive ECDF with the posterior ECDF. This is a
diagnostic for **prior-data conflict**: if the posterior sits far in the
tail of the prior predictive, the priors were inconsistent with what the
data showed. Moderate conflict is expected and healthy — it means the
data are updating the prior. Severe conflict may indicate a prior that
is too tight, a model misspecification, or anomalous data.

``` r
# Fit with the prior predictive
prior_fit    <- bayesma(data, ..., sample_prior = TRUE)

# Fit the full model
posterior_fit <- bayesma(data, ...)

# Compare
ecdf_prior_plot(posterior_fit)
```
