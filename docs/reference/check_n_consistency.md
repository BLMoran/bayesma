# Check Participant Number Consistency

Verifies that reported participant numbers are internally consistent:
randomised = analysed + lost to follow-up.

## Usage

``` r
check_n_consistency(
  n_randomised_int,
  n_randomised_ctrl,
  n_analysed_int = NULL,
  n_analysed_ctrl = NULL,
  n_lost_int = NULL,
  n_lost_ctrl = NULL,
  n_randomised_total = NULL
)
```

## Arguments

- n_randomised_int, n_randomised_ctrl:

  Integer. Per-arm randomised counts.

- n_analysed_int, n_analysed_ctrl:

  Integer or NULL. Per-arm analysed.

- n_lost_int, n_lost_ctrl:

  Integer or NULL. Per-arm lost to follow-up.

- n_randomised_total:

  Integer or NULL. Total randomised.

## Value

A list with components:

- consistent:

  Logical. TRUE if all checks pass.

- checks:

  Data frame of individual checks.

- n_checks:

  Number of checks performed.

- n_failed:

  Number of failed checks.

## Details

Implements INSPECT-SR check 4.6. Missing values (NULL) are skipped.

## Examples

``` r
check_n_consistency(
  n_randomised_int = 100, n_randomised_ctrl = 100,
  n_analysed_int = 95, n_analysed_ctrl = 92,
  n_lost_int = 5, n_lost_ctrl = 8,
  n_randomised_total = 200
)
#> $consistent
#> [1] TRUE
#> 
#> $checks
#>                                        check expected observed pass
#> 1  Total randomised = Intervention + Control      200      200 TRUE
#> 2 Intervention: Randomised = Analysed + Lost      100      100 TRUE
#> 3      Control: Randomised = Analysed + Lost      100      100 TRUE
#> 4           Intervention: Lost <= Randomised      100        5 TRUE
#> 5                Control: Lost <= Randomised      100        8 TRUE
#> 
#> $n_checks
#> [1] 5
#> 
#> $n_failed
#> [1] 0
#> 
```
