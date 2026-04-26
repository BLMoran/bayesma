# Compile and sample an Egger's test model

Compile and sample an Egger's test model

## Usage

``` r
egger_fit(
  spec,
  code = egger_stan_code(spec),
  stan_data = egger_stan_data(spec),
  chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  adapt_delta = 0.95,
  seed = 1234,
  ...
)
```

## Arguments

- spec:

  An `egger_spec` object.

- code:

  An `egger_stan_code` object. Defaults to `egger_stan_code(spec)`.

- stan_data:

  A Stan data list. Defaults to `egger_stan_data(spec)`.

- chains, iter_warmup, iter_sampling, adapt_delta, seed:

  MCMC settings.

- ...:

  Passed to `cmdstanr::CmdStanModel$sample()`.

## Value

An object of class `"egger_fit"`.
