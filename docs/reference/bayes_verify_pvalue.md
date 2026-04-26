# Bayesian P-Value Verification

Computes a Bayes factor quantifying evidence that a reported p-value is
inconsistent with the reported test statistic.

## Usage

``` r
bayes_verify_pvalue(
  test_type = c("t", "z", "chi_sq", "f"),
  statistic,
  df = NULL,
  reported_p,
  alternative = "two.sided",
  rounding_sd = 0.005,
  fabrication_sd = 0.15
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

- rounding_sd:

  Numeric. SD of the rounding error model (default 0.005).

- fabrication_sd:

  Numeric. SD under the fabrication model (default 0.15).

## Value

A list with components:

- bf_inconsistent:

  Numeric. Bayes factor for inconsistency.

- posterior_prob_inconsistent:

  Numeric. Posterior probability.

- recalculated_p:

  Numeric. Recalculated p-value.

- discrepancy:

  Numeric. Absolute difference.

- interpretation:

  Character. Evidence strength label.

## Details

Discrepancy modelled as N(0, rounding_sd^2) under H0 (honest rounding)
and N(0, fabrication_sd^2) under H1 (fabrication/error).

## Examples

``` r
bayes_verify_pvalue("chi_sq", statistic = 3.84, df = 1, reported_p = 0.05)
#> $bf_inconsistent
#> [1] 0.03333459
#> 
#> $posterior_prob_inconsistent
#> [1] 0.03225925
#> 
#> $recalculated_p
#> [1] 0.05004352
#> 
#> $discrepancy
#> [1] 4.352125e-05
#> 
#> $interpretation
#> [1] "consistent"
#> 
#> $test_type
#> [1] "chi_sq"
#> 
#> $statistic
#> [1] 3.84
#> 
```
