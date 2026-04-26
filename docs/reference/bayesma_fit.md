# Compile and sample a bayesma model

Takes the Stan code and data produced by earlier stages and runs
[`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html).
Handles the PET-PEESE decision rule (fits PET first, then switches to
PEESE when posterior evidence for an effect exceeds the threshold).

## Usage

``` r
bayesma_fit(
  spec,
  code = bayesma_stan_code(spec),
  stan_data = bayesma_stan_data(spec),
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

  A `bayesma_spec` object.

- code:

  A `bayesma_stan_code` object, or a character scalar Stan program.
  Defaults to `bayesma_stan_code(spec)`.

- stan_data:

  A Stan data list. Defaults to `bayesma_stan_data(spec)`.

- chains, iter_warmup, iter_sampling, adapt_delta, seed:

  MCMC settings.

- ...:

  Passed to `cmdstanr::CmdStanModel$sample()`.

## Value

A list with class `"bayesma_fit"` containing elements `fit` (the
cmdstanr `CmdStanMCMC`), `stan_code` (list), `stan_data` (list), and
(for PET-PEESE) `pet_peese` decision metadata.
