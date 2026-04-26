# Check Internal Consistency of Summary Statistics

Verifies relationships between reported summary statistics (CI symmetry,
SE vs CI width).

## Usage

``` r
check_statistics_consistency(
  estimate,
  ci_lower = NULL,
  ci_upper = NULL,
  se = NULL,
  log_scale = FALSE,
  ci_level = 0.95,
  tolerance = 0.1
)
```

## Arguments

- estimate:

  Numeric. Point estimate.

- ci_lower, ci_upper:

  Numeric or NULL. CI bounds.

- se:

  Numeric or NULL. Standard error.

- log_scale:

  Logical. Check on log scale (default FALSE).

- ci_level:

  Numeric. Confidence level (default 0.95).

- tolerance:

  Numeric. Tolerance (default 0.1).

## Value

A list with `consistent` and `checks`.
