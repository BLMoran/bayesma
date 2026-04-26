# Bayesian Carlisle's Test for Baseline Balance

Computes a Bayes factor comparing genuine randomisation (uniform
p-values) to fabrication (non-uniform).

## Usage

``` r
bayes_carlisle_test(p_values, prior_a = 1, prior_b = 1, n_grid = 200)
```

## Arguments

- p_values:

  Numeric vector. P-values from baseline comparisons.

- prior_a, prior_b:

  Numeric. Beta prior shape parameters (default 1, 1).

- n_grid:

  Integer. Grid size for numerical integration (default 200).

## Value

A list with components:

- bf_too_similar:

  BF for p-values biased towards 1.

- bf_too_different:

  BF for p-values biased towards 0.

- bf_nonuniform:

  Overall BF for non-uniformity.

- posterior_prob_fabrication:

  Posterior probability.

- posterior_mean_p:

  Mean of observed p-values.

- n_comparisons:

  Number of comparisons.

- interpretation:

  Evidence description.

## Examples

``` r
bayes_carlisle_test(c(0.45, 0.12, 0.78, 0.33, 0.91))
#> $bf_too_similar
#> [1] 0.07146978
#> 
#> $bf_too_different
#> [1] 0.05348527
#> 
#> $bf_nonuniform
#> [1] 0.06351762
#> 
#> $posterior_prob_fabrication
#> [1] 0.05972409
#> 
#> $posterior_mean_p
#> [1] 0.518
#> 
#> $n_comparisons
#> [1] 5
#> 
#> $interpretation
#> [1] "consistent_with_randomisation"
#> 
bayes_carlisle_test(c(0.93, 0.81, 0.95, 0.94, 0.85, 0.95))
#> $bf_too_similar
#> [1] 578.1149
#> 
#> $bf_too_different
#> [1] 0.04764676
#> 
#> $bf_nonuniform
#> [1] 287.6387
#> 
#> $posterior_prob_fabrication
#> [1] 0.9965355
#> 
#> $posterior_mean_p
#> [1] 0.905
#> 
#> $n_comparisons
#> [1] 6
#> 
#> $interpretation
#> [1] "strong_evidence_too_similar"
#> 
```
