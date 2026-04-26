# PRIMED: Preliminary Investigation of Meta-analytic Databases

Runs the four-step exploratory/preliminary data analysis workflow for
meta-analysis of dependent effect sizes as described in Pustejovsky,
Zhang, & Tipton (2026).

## Usage

``` r
primed(
  data,
  es_col,
  se_col,
  study_col,
  sample_col = NULL,
  n_col = NULL,
  moderators = NULL,
  es_type = c("SMD", "correlation"),
  df_col = NULL,
  rho_values = c(0.1, 0.3, 0.5, 0.7, 0.9),
  fence_multiplier = 3
)
```

## Arguments

- data:

  A data frame containing the meta-analytic database.

- es_col:

  Character. Name of the column containing effect size estimates.

- se_col:

  Character. Name of the column containing standard errors.

- study_col:

  Character. Name of the column identifying studies.

- sample_col:

  Character. Name of the column identifying samples within studies. If
  NULL (default), assumes one sample per study.

- n_col:

  Character. Name of the column containing total sample sizes. If NULL,
  sample size plots are skipped.

- moderators:

  Character vector of column names to examine as potential moderators.
  If NULL (default), no moderator analysis is performed.

- es_type:

  Character. Type of effect size: "SMD" (standardized mean difference)
  or "correlation". Affects how scaled SEs and weights are computed.
  Default is "SMD".

- df_col:

  Character. Name of the column containing degrees of freedom (used for
  scaled SE calculation when es_type = "SMD"). If NULL, scaled SEs are
  not computed.

- rho_values:

  Numeric vector. Assumed within-sample correlations for ISC weight
  calculations. Default is c(0.1, 0.3, 0.5, 0.7, 0.9).

- fence_multiplier:

  Numeric. Multiplier of the IQR for outlier fences in effect size
  density plots. Default is 3 (following Tukey's conventions as
  described in the paper).

## Value

A named list with elements:

- summary:

  A list of summary statistics about the database.

- plots:

  A named list of ggplot objects for each workflow step.

- tables:

  A named list of summary tables (tibbles).

## Details

The PRIMED workflow proceeds in four steps:

1.  **Data structure**: Counts observations at each level and describes
    the dependence structure (effects per sample, samples per study,
    sample size distributions).

2.  **Moderators**: Examines marginal distributions, missingness, and
    hierarchical (within- vs between-sample) structure of covariates.

3.  **Standard errors & weights**: Inspects SE distributions within
    samples, computes scaled SEs (for SMDs), and calculates ISC weights
    under varying assumed correlations.

4.  **Effect size distribution**: Visualises marginal and sample-level
    densities with outlier fences, and produces a hierarchical forest
    plot of dependent effect sizes.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- primed(
  data = my_meta_data,
  es_col = "g",
  se_col = "se",
  study_col = "study",
  sample_col = "sample_id",
  n_col = "n_total",
  moderators = c("intervention_type", "mean_age", "pct_female"),
  es_type = "SMD",
  df_col = "df"
)

# View all step-1 plots
results$plots$step1_es_per_sample
results$plots$step1_samples_per_study

# Access summary statistics
results$summary
} # }
```
