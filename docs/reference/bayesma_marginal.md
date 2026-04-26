# Compute a marginal estimand from a bayesma fit

Post-processes the posterior draws of a fitted
[`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md)
model to return a marginal estimand on the natural scale: risk
difference (RD/ARR), average treatment effect (ATE), average treatment
effect on the treated (ATT), or a conditional average treatment effect
(CATE).

## Usage

``` r
bayesma_marginal(fit, spec)
```

## Arguments

- fit:

  A `bayesma_fit` object (or the cmdstanr fit inside one).

- spec:

  The matching `bayesma_spec`.

## Value

A list with elements:

- `estimand`:

  The estimand label.

- `draws`:

  A numeric vector of posterior draws on the natural scale.

- `summary`:

  A tibble with median, 95\\ of being above zero.

## Details

For relative-effect estimands (`"OR"`, `"RR"`, `"HR"`, `"IRR"`, `"MD"`,
`"SMD"`) this function is unnecessary —
[`bayesma_output()`](https://blmoran.github.io/bayesma/reference/bayesma_output.md)
already returns the pooled effect on the appropriate scale.

## Methods by estimand and stage

- RD / ARR / ATE — binomial, one-stage:

  Computed via marginal standardisation (g-computation) over posterior
  draws of the per-study baseline logit (`gamma[s]`) and the pooled
  log-OR (`mu`), optionally shifted by study-level random effects
  (`epsilon[s]`). The per-study RD is
  `plogis(gamma[s] + mu + epsilon[s]) - plogis(gamma[s])`, averaged
  across studies weighted by the harmonic mean of arm sample sizes. This
  corresponds to a population-weighted ATE over the observed study mix
  and requires no external baseline assumption.

- RD / ARR / ATE — binomial, two-stage:

  Back-transforms posterior draws of the pooled log-OR using a baseline
  risk drawn per-iteration from a Beta distribution fitted
  (method-of-moments) to the observed control-arm event rates.
  Propagates baseline uncertainty into the posterior RD. A fixed scalar
  `baseline_risk` bypasses this and uses the supplied value directly
  (old behaviour).

- ATE — gaussian:

  Equivalent to MD. Returns the pooled posterior on the absolute scale.

- ATT:

  One-stage: weighted by intervention-arm sample size
  (`n_int / sum(n_int)`). Two-stage: same back-transform as ATE but with
  intervention-size-weighted baseline. Without IPD this is an
  arm-size-weighted ATE on the treated, not a true causal ATT —
  interpret with caution.

- CATE:

  Routes to
  [`meta_reg()`](https://blmoran.github.io/bayesma/reference/meta_reg.md)
  with `moderators = cate_covariate`. Reports the meta-regression
  effect; the user is expected to evaluate it at a specific covariate
  value downstream.
