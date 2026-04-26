# Bayesian Evidence for Meta-Regression Coefficients

Computes Bayesian measures of evidence for moderator coefficients,
including probability of direction, null range (Region of Practical
Equivalence), and credible intervals.

## Usage

``` r
coefficient_evidence(object, null_range = NULL, ci_level = 0.95)
```

## Arguments

- object:

  A `bayesma_reg` object.

- null_range:

  Numeric vector of length 2 defining the null range, or NULL (default)
  to skip null range calculation. For standardized coefficients,
  `c(-0.1, 0.1)` is a common choice.

- ci_level:

  Numeric. Credible interval width (default: 0.95).

## Value

A tibble with columns:

- `term`: Coefficient name

- `estimate`: Posterior median

- `std_error`: Posterior SD

- `ci_lower`, `ci_upper`: Credible interval bounds

- `pd`: Probability of direction (max of P(β\>0), P(β\<0))

- `direction`: Most probable direction ("positive" or "negative")

- `p_null`: Proportion of posterior inside null range (if specified

- `p_outside_null`: Proportion outside null range (if specified)

## Details

### Probability of Direction (pd)

The probability of direction is the proportion of the posterior on the
same side of zero as the median. It ranges from 0.5 (no evidence) to 1
(strong evidence for direction). A pd \> 0.95 is often considered
meaningful evidence, though this is not a formal threshold.

### null range (Region of Practical Equivalence)

The null range defines a range of values considered "practically
equivalent to zero". The proportion of the posterior inside the null
range indicates whether the effect is practically meaningful:

- `p_null` \> 0.95: Effect is practically zero (accept null)

- `p_null` \< 0.05: Effect is practically meaningful (reject null)

- Otherwise: Inconclusive

Common null range choices:

- `c(-0.1, 0.1)` for standardized coefficients (Cohen's d scale)

- `c(-0.05, 0.05)` for log-OR (corresponds to OR 0.95-1.05)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- meta_reg(data, studyvar = "author", yi = "yi", vi = "vi",
                mods = ~ year + quality)

# Basic evidence summary
coefficient_evidence(fit)

# With null range
coefficient_evidence(fit, null_range = c(-0.1, 0.1))
} # }
```
