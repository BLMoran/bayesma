# GRIM Test (Granularity-Related Inconsistency of Means)

Tests whether a reported mean is mathematically possible given the
sample size and the number of decimal places reported. Applies to data
measured on an integer scale (e.g., Likert items, counts).

## Usage

``` r
grim_test(mean_value, n, decimals = NULL, tolerance = 1e-06)
```

## Arguments

- mean_value:

  Numeric. The reported mean.

- n:

  Integer. The sample size.

- decimals:

  Integer. Number of decimal places in the reported mean. If NULL
  (default), inferred from the reported value.

- tolerance:

  Numeric. Rounding tolerance for comparison (default 1e-6).

## Value

A list with components:

- consistent:

  Logical. TRUE if the mean is GRIM-consistent.

- mean_value:

  The tested mean.

- n:

  The sample size.

- decimals:

  Number of decimal places used.

## Details

The GRIM test (Brown & Heathers, 2017) checks whether a reported mean of
integer data is consistent with the reported sample size. For example,
with N = 20, a mean must be a multiple of 1/20 = 0.05. A reported mean
of 3.47 would be impossible.

This implements INSPECT-SR check 4.8.

## References

Brown NJL, Heathers JAJ (2017). The GRIM test: A simple technique
detects numerous anomalies in the reporting of results in psychology.
*Social Psychological and Personality Science*, 8(4), 363-369.

## Examples

``` r
# Possible mean: 52/20 = 2.60
grim_test(2.60, n = 20)
#> $consistent
#> [1] TRUE
#> 
#> $mean_value
#> [1] 2.6
#> 
#> $n
#> [1] 20
#> 
#> $decimals
#> [1] 1
#> 

# Impossible mean: no integer sum / 20 = 2.47
grim_test(2.47, n = 20)
#> $consistent
#> [1] FALSE
#> 
#> $mean_value
#> [1] 2.47
#> 
#> $n
#> [1] 20
#> 
#> $decimals
#> [1] 2
#> 
```
