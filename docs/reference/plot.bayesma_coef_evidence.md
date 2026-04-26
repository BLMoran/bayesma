# Plot Coefficient Evidence

Creates a forest-style plot of coefficient posteriors with credible
intervals and optional null range visualization.

## Usage

``` r
# S3 method for class 'bayesma_coef_evidence'
plot(x, include_intercept = FALSE, show_null_range = TRUE, ...)
```

## Arguments

- x:

  A `bayesma_coef_evidence` object.

- include_intercept:

  Logical. Include intercept in plot (default: FALSE).

- show_null_range:

  Logical. Show null range region if available (default: TRUE).

- ...:

  Additional arguments (unused).

## Value

A ggplot object.
