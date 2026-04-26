# Verify a Reported P-Value

Recalculates a p-value from a reported test statistic and degrees of
freedom, then compares it to the reported p-value.

## Usage

``` r
verify_pvalue(
  test_type = c("t", "z", "chi_sq", "f"),
  statistic,
  df = NULL,
  reported_p,
  alternative = "two.sided",
  tolerance = 0.01
)
```

## Arguments

- test_type:

  Character. One of `"t"`, `"z"`, `"chi_sq"`, `"f"`.

- statistic:

  Numeric. The reported test statistic.

- df:

  Numeric. Degrees of freedom. For F-tests, a vector of length 2. Not
  required for z-tests.

- reported_p:

  Numeric. The reported p-value.

- alternative:

  Character. One of `"two.sided"` (default), `"less"`, `"greater"`.

- tolerance:

  Numeric. Acceptable absolute difference (default 0.01).

## Value

A list with components:

- consistent:

  Logical. TRUE if p-values match within tolerance.

- reported_p:

  The reported p-value.

- recalculated_p:

  The recalculated p-value.

- difference:

  Absolute difference.

- test_type:

  The test type used.

- statistic:

  The test statistic.

## Details

Implements INSPECT-SR check 4.9.

## Examples

``` r
# Consistent: chi-squared = 3.84, df = 1, p = 0.05
verify_pvalue("chi_sq", statistic = 3.84, df = 1, reported_p = 0.05)
#> $consistent
#> [1] TRUE
#> 
#> $reported_p
#> [1] 0.05
#> 
#> $recalculated_p
#> [1] 0.05004352
#> 
#> $difference
#> [1] 4.352125e-05
#> 
#> $test_type
#> [1] "chi_sq"
#> 
#> $statistic
#> [1] 3.84
#> 
```
