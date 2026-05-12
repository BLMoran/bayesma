# RoBMA bias-prior constructors

Constructors for the publication-bias priors used by
[`robma()`](https://blmoran.github.io/bayesma/reference/robma.md) and
[`robma_sensitivity()`](https://blmoran.github.io/bayesma/reference/robma_sensitivity.md).
Each returns a `robma_bias_prior` object.

## Usage

``` r
prior_bias(
  type = c("weight_function", "pet", "peese", "copas", "jung", "none"),
  parameters = list(),
  prior_weight = 1,
  ...
)

prior_weight_function(
  steps = c(0.025, 0.05),
  alpha = NULL,
  prior_weight = 1,
  sided = "one"
)

prior_pet(distribution = "normal", location = 0, scale = 1, prior_weight = 1)

prior_peese(distribution = "normal", location = 0, scale = 2, prior_weight = 1)

prior_copas(prior_weight = 1)

prior_jung(prior_weight = 1)

prior_no_bias(prior_weight = 1)
```

## Arguments

- type:

  Bias-prior family. One of `"weight_function"`, `"pet"`, `"peese"`,
  `"copas"`, `"jung"`, `"none"`.

- parameters:

  Named list of family-specific parameters.

- prior_weight:

  Numeric. Prior model weight in the bias-prior mixture.

- ...:

  Additional fields stored on the bias-prior object.

- steps:

  Numeric vector of p-value cutpoints for a step weight function.

- alpha:

  Numeric vector of Dirichlet concentration parameters; defaults to
  `rep(1, length(steps) + 1)`.

- sided:

  One of `"one"` or `"two"`. Determines whether the weight function is
  one- or two-sided.

- distribution:

  Distribution for the PET / PEESE slope prior. Currently `"cauchy"`.

- location, scale:

  Location and scale of the PET / PEESE slope prior.

## Value

A `robma_bias_prior` object.
