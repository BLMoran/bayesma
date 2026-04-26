# Extract draws from a bayesma object in bayesfoRest-compatible format

Reshapes the posterior draws from a bayesma fit into the long-format
tibble that bayesfoRest's forest.data_fn expects.

## Usage

``` r
extract_forest_draws(x)
```

## Arguments

- x:

  A bayesma object

## Value

A tibble of posterior draws in long format

## Details

Columns produced:

- Author (character): study label or "Pooled Effect" / "Prediction"

- b_Intercept (numeric): study-level effect draw (mu + random effect)

- r_Author (numeric): random effect deviation from pooled

- sd_Author\_\_Intercept (numeric): tau draw

- .chain, .iteration, .draw: MCMC identifiers
