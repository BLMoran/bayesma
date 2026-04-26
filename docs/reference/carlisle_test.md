# Carlisle's Test for Baseline Balance

Tests whether the distribution of p-values for baseline comparisons is
consistent with genuine randomisation.

## Usage

``` r
carlisle_test(p_values, method = c("fisher", "ks"))
```

## Arguments

- p_values:

  Numeric vector. P-values from baseline comparisons.

- method:

  Character. `"fisher"` (default) or `"ks"`.

## Value

A list with components:

- too_similar:

  Logical. Suspiciously well-balanced.

- too_different:

  Logical. Suspiciously imbalanced.

- combined_p:

  Combined p-value.

- n_comparisons:

  Number of comparisons.

- method:

  Method used.

- interpretation:

  `"plausible"`, `"too_similar"`, or `"too_different"`.

## Details

Implements INSPECT-SR check 4.3. Under genuine randomisation, baseline
p-values should be approximately uniform. Fabricated trials often show
implausibly well-matched groups (p-values near 1).

## References

Carlisle JB (2017). Data fabrication and other reasons for non-random
sampling in 5087 randomised, controlled trials in anaesthetic and
general medical journals. *Anaesthesia*, 72(8), 944-952.

## Examples

``` r
carlisle_test(c(0.45, 0.12, 0.78, 0.33, 0.91))
#> $too_similar
#> [1] FALSE
#> 
#> $too_different
#> [1] FALSE
#> 
#> $combined_p
#> [1] 0.443096
#> 
#> $n_comparisons
#> [1] 5
#> 
#> $method
#> [1] "fisher"
#> 
#> $interpretation
#> [1] "plausible"
#> 
carlisle_test(c(0.92, 0.88, 0.95, 0.91, 0.87))
#> $too_similar
#> [1] TRUE
#> 
#> $too_different
#> [1] FALSE
#> 
#> $combined_p
#> [1] 0.00016601
#> 
#> $n_comparisons
#> [1] 5
#> 
#> $method
#> [1] "fisher"
#> 
#> $interpretation
#> [1] "too_similar"
#> 
```
