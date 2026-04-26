# Compute null range probabilities

Compute null range probabilities

## Usage

``` r
compute_null_range_probs(mu_draws, null_range = NULL, effect_label = "log_or")
```

## Arguments

- mu_draws:

  Numeric vector of model-averaged posterior draws (on the log scale for
  ratio measures).

- null_range:

  Numeric vector of length 2. Can be specified on either the log scale
  or the natural scale for ratio measures. Auto-detected: if the range
  straddles 1 (e.g., c(0.9, 1.1)), it's assumed to be on the natural
  (OR/RR) scale and is log-transformed internally. If it straddles 0
  (e.g., c(-0.1, 0.1)), it's assumed to be on the log scale already.

- effect_label:

  Character. Effect scale identifier.

## Value

A list with p_negative, p_null, p_positive, null_range (on log scale),
null_range_natural (on natural scale).
