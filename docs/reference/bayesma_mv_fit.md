# Compile and sample a bivariate meta-analysis model

Compile and sample a bivariate meta-analysis model

## Usage

``` r
bayesma_mv_fit(
  spec,
  code = bayesma_mv_stan_code(spec),
  stan_data = bayesma_mv_stan_data(spec),
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  seed = 1234,
  ...
)
```

## Arguments

- spec:

  A `bayesma_mv_spec` object.

- code:

  A `bayesma_mv_stan_code` object. Defaults to
  `bayesma_mv_stan_code(spec)`.

- stan_data:

  A Stan data list. Defaults to `bayesma_mv_stan_data(spec)`.

- chains, iter_warmup, iter_sampling, adapt_delta, seed:

  MCMC settings.

- ...:

  Passed to `cmdstanr::CmdStanModel$sample()`.

## Value

An object of class `"bayesma_mv_fit"`.
