# Compile and sample a meta-regression model

Compile and sample a meta-regression model

## Usage

``` r
meta_reg_fit(
  spec,
  code = meta_reg_stan_code(spec),
  stan_data = meta_reg_stan_data(spec),
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

  A `meta_reg_spec` object.

- code:

  A `meta_reg_stan_code` object. Defaults to `meta_reg_stan_code(spec)`.

- stan_data:

  A Stan data list. Defaults to `meta_reg_stan_data(spec)`.

- chains, iter_warmup, iter_sampling, adapt_delta, seed:

  MCMC settings.

- ...:

  Passed to `cmdstanr::CmdStanModel$sample()`.

## Value

An object of class `"meta_reg_fit"`.
