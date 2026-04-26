# Per-Check INSPECT-SR Results Table

Unpacks the automated Domain 4 checks from an
[`inspect_sr()`](https://blmoran.github.io/bayesma/reference/inspect_sr.md)
result into a publication-ready
[`gt::gt()`](https://gt.rstudio.com/reference/gt.html) table, grouped by
study with visual separators.

## Usage

``` r
inspect_sr_table(
  x,
  study = NULL,
  check = c("grim", "pvalue", "carlisle", "n_consistency"),
  only_failed = FALSE
)
```

## Arguments

- x:

  An object returned by
  [`inspect_sr()`](https://blmoran.github.io/bayesma/reference/inspect_sr.md).

- study:

  Optional character vector of study names to restrict the table to
  (useful for per-study tabsets in Quarto). Default `NULL` = all
  studies.

- check:

  Character vector. Restrict to specific checks: `"grim"`, `"pvalue"`,
  `"carlisle"`, `"n_consistency"`. Default: all.

- only_failed:

  Logical. If `TRUE`, show only flagged rows. Default `FALSE`.

## Value

A `gt` table object.

## Examples

``` r
if (FALSE) { # \dontrun{
data(inspect_sr_example)
res <- inspect_sr(inspect_sr_example, verbose = FALSE)
inspect_sr_table(res)
inspect_sr_table(res, only_failed = TRUE)
inspect_sr_table(res, study = "Doe (1995)")
} # }
```
