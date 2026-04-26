# Bayesian GRIM Test

Computes a Bayes factor quantifying evidence that a reported mean of
integer-scale data is inconsistent with the reported sample size,
accounting for rounding uncertainty.

## Usage

``` r
bayes_grim_test(
  mean_value,
  n,
  decimals = NULL,
  max_items = 10,
  n_tolerance = 0
)
```

## Arguments

- mean_value:

  Numeric. The reported mean.

- n:

  Integer. The sample size.

- decimals:

  Integer. Number of decimal places in the reported mean. If NULL
  (default), inferred from the reported value.

- max_items:

  Integer. Maximum plausible value on the integer scale (default 10).
  Used to define the range of possible integer sums.

- n_tolerance:

  Integer. How many values around the reported N to consider as
  plausible (default 0, exact N only).

## Value

A list with components:

- bf_inconsistent:

  Numeric. Bayes factor in favour of inconsistency (fabrication).

- posterior_prob_inconsistent:

  Numeric. Posterior probability of inconsistency assuming equal prior
  odds.

- consistent_at_n:

  Logical. Classical GRIM result at exact N.

- consistent_nearby:

  Logical. GRIM-consistent at any nearby N.

- interpretation:

  Character. Evidence strength label.

## Details

Under H0 (genuine data), the mean must equal k/n for some integer k,
rounded to the reported decimal places. Under H1 (fabricated data), the
mean is drawn uniformly from the plausible range. BF_10 = P(data \| H1)
/ P(data \| H0).

## Examples

``` r
# Consistent mean: evidence for genuine data
bayes_grim_test(2.60, n = 20)
#> $bf_inconsistent
#> [1] 1.005
#> 
#> $posterior_prob_inconsistent
#> [1] 0.5012469
#> 
#> $consistent_at_n
#> [1] TRUE
#> 
#> $consistent_nearby
#> [1] TRUE
#> 
#> $interpretation
#> [1] "weak_evidence"
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

# Inconsistent mean: evidence for fabrication
bayes_grim_test(2.47, n = 20)
#> $bf_inconsistent
#> [1] Inf
#> 
#> $posterior_prob_inconsistent
#> [1] 1
#> 
#> $consistent_at_n
#> [1] FALSE
#> 
#> $consistent_nearby
#> [1] FALSE
#> 
#> $interpretation
#> [1] "strong_inconsistency"
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
