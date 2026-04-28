# Stan code — overview and conventions

## Introduction

**bayesma** fits all models in Stan via the **cmdstanr** interface. This
vignette explains how the interface works, documents the notation
conventions used across the Stan Code vignettes, and describes how to
extract and modify generated Stan code.

## How bayesma interfaces with Stan

Each call to
[`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md)
runs a six-stage internal pipeline:

    spec → code → data → fit → extract → output

1.  **spec** (`bayesma_spec()`): validates arguments and resolves model
    type, likelihood, stage, RE distribution, and prior specifications
    into a structured list.
2.  **code** (`bayesma_stan_code()`): generates the Stan program as a
    character string from the spec.
3.  **data** (`bayesma_stan_data()`): prepares the named list of data
    passed to `cmdstanr::cmdstan_model()$sample()`.
4.  **fit**: compiles the Stan model and runs MCMC. Compilation is
    cached.
5.  **extract** (`bayesma_extract()`): reshapes the posterior draws into
    tidy format.
6.  **output** (`bayesma_output()`): formats the posterior for display.

Any stage can be returned early for inspection:

``` r
spec  <- bayesma(data, model_type = "random_effect", return_stage = "spec")
code  <- bayesma(data, model_type = "random_effect", return_stage = "code")
sdata <- bayesma(data, model_type = "random_effect", return_stage = "data")
fit   <- bayesma(data, model_type = "random_effect", return_stage = "fit")
```

`return_stage = "code"` returns the Stan program as a character string;
print it with [`cat()`](https://rdrr.io/r/base/cat.html).

## Notation conventions

The Stan Code vignettes use the following notation consistently:

| Symbol       | Meaning                                   |
|--------------|-------------------------------------------|
| $`N`$        | Number of observations (rows in the data) |
| $`K`$        | Number of studies                         |
| $`y_i`$      | Observed effect in study $`i`$            |
| $`s_i`$      | Standard error of $`y_i`$                 |
| $`\theta_i`$ | Study-level true effect                   |
| $`\mu`$      | Population-level pooled effect            |
| $`\tau`$     | Between-study heterogeneity (SD)          |
| $`u_i`$      | Study-level random effect                 |
| $`z_i`$      | Non-centred auxiliary variable            |

Stan variable names follow this convention wherever possible.

## Non-centred parameterisation

All random-effects models in **bayesma** use the **non-centred
parameterisation** (NCP) for the study-level random effects:

``` stan
parameters {
  real mu;
  real<lower=0> tau;
  vector[K] z;          // standard normal auxiliary
}

transformed parameters {
  vector[K] u = tau * z;  // re-centred in transformed parameters
}
```

The NCP separates the scale ($`\tau`$) from the shape ($`z`$) of the
random effects. It is equivalent to the centred parameterisation

``` stan
parameters {
  real mu;
  real<lower=0> tau;
  vector[K] u;
}
model {
  u ~ normal(0, tau);
}
```

but samples more efficiently when $`\tau`$ is small relative to its
prior, because the geometry of the posterior near $`\tau = 0`$ is better
behaved in the NCP.

## Extracting and modifying Stan code

``` r
code <- bayesma(data, model_type = "random_effect", return_stage = "code")
cat(code)
```

The returned string is valid Stan. It can be edited and passed back via
`custom_model`:

``` r
modified <- stringr::str_replace(
  code,
  "target \\+= cauchy_lpdf\\(tau \\| 0, 0.5\\);",
  "target += half_normal_lpdf(tau | 0, 0.25);"
)

fit <- bayesma(data, custom_model = modified)
```

For larger modifications, write the custom code to a `.stan` file:

``` r
writeLines(code, "my_model.stan")
# ... edit my_model.stan ...
fit <- bayesma(data, custom_model = readLines("my_model.stan"))
```

## General Stan tips for meta-analysis

**Divergent transitions.** If `adapt_delta` warnings appear, increase
`adapt_delta`:

``` r
bayesma(data, adapt_delta = 0.99)
```

Persistent divergences after increasing `adapt_delta` indicate a
pathological posterior geometry, often caused by a funnel-shaped
posterior near $`\tau = 0`$. Switching to the NCP (already the default)
or placing a more informative prior on $`\tau`$ resolves most cases.

**Effective sample size.** Aim for bulk ESS \> 400 per parameter. Low
ESS for $`\tau`$ is common when $`k`$ is small and the posterior for
$`\tau`$ is near-zero; this does not indicate a bug but it does inflate
Monte Carlo error for tail quantities.

**Rhat.** All $`\hat{R} < 1.01`$ is the requirement for reliable
inference. Values above 1.05 indicate the chains have not mixed.

**Compilation caching.** Stan compiles each model to a C++ binary on
first use. Subsequent calls with the same model are fast. If the model
code changes (including whitespace in some cmdstanr versions),
recompilation is triggered.
