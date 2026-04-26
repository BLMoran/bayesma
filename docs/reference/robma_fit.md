# Compile and sample a RoBMA specification

Accepts the outputs of
[`robma_stan_code()`](https://blmoran.github.io/bayesma/reference/robma_stan_code.md)
and
[`robma_stan_data()`](https://blmoran.github.io/bayesma/reference/robma_stan_data.md).
Dispatches on `spec$method` to bridge sampling (with model averaging
across the component grid) or the spike-and-slab joint model.

## Usage

``` r
robma_fit(
  spec,
  code,
  stan_data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  seed = 1234,
  parallel = FALSE,
  quiet = FALSE,
  ...
)
```

## Arguments

- spec:

  A `bayesma_robma_spec`.

- code:

  A `bayesma_robma_stan_code`.

- stan_data:

  A `bayesma_robma_stan_data`.

- chains, iter_warmup, iter_sampling, adapt_delta, seed:

  MCMC settings.

- parallel:

  Logical. Parallelise sampling across component models (bridge only).

- quiet:

  Logical. Suppress progress messages.

- ...:

  Passed to `cmdstanr::CmdStanModel$sample()`.

## Value

An object of class `"bayesma_robma_fit"`.
