# Create a Table for Model Comparison Results

Produces a publication-ready `gt` table with LOSO-CV metrics (when
available) and within-stage LOO diagnostics.

## Usage

``` r
compare_table(x, digits = 2, include_loo = TRUE, include_diagnostics = TRUE)
```

## Arguments

- x:

  A `bayesma_comparison` object.

- digits:

  Integer. Decimal places. Default 2.

- include_loo:

  Logical. Include within-stage LOO columns. Default `TRUE`.

- include_diagnostics:

  Logical. Include Pareto k column. Default `TRUE`.

## Value

A `gt` table object.

## Examples

``` r
if (FALSE) { # \dontrun{
compare_table(comparison)
compare_table(comparison) |> gt::gtsave("comparison.html")
} # }
```
