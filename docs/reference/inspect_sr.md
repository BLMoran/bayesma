# Run the INSPECT-SR Trustworthiness Assessment

Takes a data frame with one row per study and runs the automated Domain
4 checks (Carlisle's test, participant-number consistency, GRIM, p-value
verification). Manual items (D1, D2, D3, and the non-automated D4 items)
are read straight from the input. Domain-level and overall judgements
are derived per INSPECT-SR guidance (overall = most severe domain).

## Usage

``` r
inspect_sr(
  data,
  studyvar = study,
  bayes = FALSE,
  prior_prob_trustworthy = 0.9,
  pvalue_tolerance = 0.01,
  carlisle_method = "fisher",
  verbose = TRUE
)
```

## Arguments

- data:

  A data frame or tibble with one row per study. See **Expected
  columns** in the package vignette, or the bundled
  [inspect_sr_example](https://blmoran.github.io/bayesma/reference/inspect_sr_example.md)
  dataset for the exact layout.

- studyvar:

  Unquoted column name identifying the study (tidyeval). Defaults to
  `study`.

- bayes:

  Logical. If `FALSE` (default) produces frequentist pass/fail
  judgements. If `TRUE` produces Bayes factors and a posterior
  probability of trustworthiness.

- prior_prob_trustworthy:

  Numeric in (0, 1). Prior probability that each study is trustworthy,
  used only when `bayes = TRUE` (default 0.90).

- pvalue_tolerance:

  Numeric. Tolerance for the frequentist p-value check (default 0.01).

- carlisle_method:

  `"fisher"` (default) or `"ks"` — see
  [`carlisle_test()`](https://blmoran.github.io/bayesma/reference/carlisle_test.md).

- verbose:

  Logical. Print a summary to the console (default `TRUE`).

## Value

If `bayes = FALSE`: an object of class `inspect_sr` (a data frame with
columns `Study`, `D1`, `D2`, `D3`, `D4`, `Overall`), with per-study
details in `attr(x, "details")`.

If `bayes = TRUE`: an object of class `bayes_inspect_sr` (a data frame
with columns `Study`, `Prior`, `Posterior`, `Combined_BF`,
`Interpretation`), with individual Bayes factors in
`attr(x, "details")`.

## See also

[`inspect_sr_table()`](https://blmoran.github.io/bayesma/reference/inspect_sr_table.md)
for a per-check gt table;
[`inspect_plot()`](https://blmoran.github.io/bayesma/reference/inspect_plot.md)
for the traffic-light visualisation;
[`filter_trustworthy()`](https://blmoran.github.io/bayesma/reference/filter_trustworthy.md)
for filtering a meta-analysis dataset.

## Examples

``` r
if (FALSE) { # \dontrun{
data(inspect_sr_example)

# Frequentist
res <- inspect_sr(inspect_sr_example, studyvar = study)

# Bayesian
res_bayes <- inspect_sr(inspect_sr_example, studyvar = study, bayes = TRUE)
} # }
```
