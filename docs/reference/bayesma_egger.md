# Bayesian Egger's Test for Small-Study Effects

Implements a Bayesian approach to assessing small-study effects (funnel
plot asymmetry) in meta-analysis, based on Shi, Chu & Lin (2020). This
method controls false positive rates by using latent "true" standard
errors rather than sample standard errors in an Egger-type regression.

The standard Egger's test can have inflated false positive rates when
effect sizes (particularly odds ratios) are intrinsically correlated
with their standard errors. This Bayesian method addresses this
limitation by modeling the true within-study variances through a
hierarchical model.

## Usage

``` r
bayesma_egger(
  data,
  study,
  event_ctrl = NULL,
  event_int = NULL,
  n_ctrl,
  n_int,
  mean_ctrl = NULL,
  mean_int = NULL,
  sd_ctrl = NULL,
  sd_int = NULL,
  likelihood = c("binomial", "gaussian", "poisson"),
  heterogeneity = c("multiplicative", "additive"),
  alpha_prior = NULL,
  beta_prior = NULL,
  kappa_prior = NULL,
  gamma_prior = NULL,
  d_prior = NULL,
  tau_prior = NULL,
  credible_level = 0.9,
  chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  adapt_delta = 0.95,
  seed = 1234,
  ...
)
```

## Arguments

- data:

  A data frame containing the study data.

- study:

  Character. Column name for study identifiers.

- event_ctrl, event_int:

  Character. Column names for event counts in control and intervention
  groups. Required for binary outcomes.

- n_ctrl, n_int:

  Character. Column names for sample sizes in control and intervention
  groups.

- mean_ctrl, mean_int:

  Character. Column names for means in control and intervention groups.
  Required for continuous outcomes.

- sd_ctrl, sd_int:

  Character. Column names for standard deviations. Required for
  continuous outcomes.

- likelihood:

  Character. `"binomial"` (default), `"gaussian"`, or `"poisson"`.

- heterogeneity:

  Character. `"multiplicative"` (default) or `"additive"`. Specifies the
  heterogeneity structure in the Egger regression.

  - `"multiplicative"`: Assumes `Var(y_i) = kappa^2 * sigma_i^2`

  - `"additive"`: Assumes `Var(y_i) = sigma_i^2 + gamma^2`

- alpha_prior:

  Prior on the regression intercept. Default: `normal(0, 100)`.

- beta_prior:

  Prior on the regression slope (small-study effect parameter). Default:
  `normal(0, 100)`.

- kappa_prior:

  Prior on multiplicative heterogeneity SD (when
  `heterogeneity = "multiplicative"`). Default: `uniform(0, 2)`.

- gamma_prior:

  Prior on additive heterogeneity SD (when
  `heterogeneity = "additive"`). Default: `uniform(0, 2)`.

- d_prior:

  Prior on the overall log OR (mu). Default: `normal(0, 100)`.

- tau_prior:

  Prior on between-study heterogeneity SD. Default: `uniform(0, 2)`.

- credible_level:

  Numeric. Credible interval level for inference. Default: 0.90 (90%
  CrI, consistent with typical small-study effects testing).

- chains:

  Integer. Number of MCMC chains. Default: 4.

- iter_warmup:

  Integer. Number of warmup iterations per chain. Default: 2000.

- iter_sampling:

  Integer. Number of sampling iterations per chain. Default: 4000.

- adapt_delta:

  Numeric. Target acceptance rate. Default: 0.95.

- seed:

  Integer. Random seed for reproducibility. Default: 1234.

- ...:

  Additional arguments passed to
  [`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html).

## Value

A list of class `"bayesma_egger"` containing:

- fit:

  The cmdstanr fit object

- summary:

  Summary statistics for key parameters

- draws:

  Posterior draws as a data frame

- beta_summary:

  Summary of the slope parameter (small-study effects)

- conclusion:

  Character string interpreting the result

- meta:

  List of metadata including priors, settings, and effect sizes

- stan_code:

  The generated Stan code

- stan_data:

  The data passed to Stan

## Details

### Model Specification

The method uses a two-level hierarchical model:

**Level 1 (Meta-analysis model):** For binary outcomes with odds ratios:
\$\$n\_{ij1} \sim \text{Binomial}(n\_{ij\cdot}, p\_{ij})\$\$
\$\$\text{logit}(p\_{ij}) = \mu_i + \delta_i \cdot I(j=1)\$\$
\$\$\delta_i \sim N(d, \tau^2)\$\$

**Level 2 (Egger regression):** With multiplicative heterogeneity:
\$\$y_i \sim N(\alpha + \beta \sigma_i, \kappa^2 \sigma_i^2)\$\$

Or with additive heterogeneity: \$\$y_i \sim N(\alpha + \beta \sigma_i,
\sigma_i^2 + \gamma^2)\$\$

Where \\\sigma_i\\ is the latent "true" standard error derived from the
posterior distributions of \\p\_{i0}\\ and \\p\_{i1}\\.

### Interpretation

The slope parameter \\\beta\\ indicates the strength of small-study
effects:

- If the credible interval for \\\beta\\ excludes 0, there is evidence
  of small-study effects (funnel plot asymmetry).

- Positive \\\beta\\: Studies with smaller precision tend to show larger
  positive effects (potential missing negative studies).

- Negative \\\beta\\: Studies with smaller precision tend to show larger
  negative effects (potential missing positive studies).

## References

Shi L, Chu H, Lin L. (2020). A Bayesian approach to assessing
small-study effects in meta-analysis of a binary outcome with controlled
false positive rate. *Research Synthesis Methods*, 11(4):535-552.
[doi:10.1002/jrsm.1415](https://doi.org/10.1002/jrsm.1415)

## Examples

``` r
# Example with binary outcome data
result <- bayesma_egger(
  data = my_meta_data,
  study = "study_id",
  event_ctrl = "events_control",
  event_int = "events_treatment",
  n_ctrl = "n_control",
  n_int = "n_treatment",
  likelihood = "binomial"
)
#> Error: object 'my_meta_data' not found

# View results
print(result)
#> Error: object 'result' not found
summary(result)
#> Error: object 'result' not found

```
