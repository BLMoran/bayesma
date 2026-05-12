# Prior distribution constructors

Constructors for prior distributions used by
[`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md),
[`robma()`](https://blmoran.github.io/bayesma/reference/robma.md),
[`meta_reg()`](https://blmoran.github.io/bayesma/reference/meta_reg.md),
and related fitting functions. Each returns an object of class
`bayesma_prior` that stores the family and its hyperparameters.

## Usage

``` r
normal(mean, sd)

half_normal(mean, sd)

cauchy(location, scale)

half_cauchy(location, scale)

half_student_t(df, location, scale)

exponential(rate)

uniform(lower, upper)

dirichlet(alpha)

beta(alpha, beta)

scaled_inv_chi_sq(df, scale)

lkj(eta)
```

## Arguments

- mean:

  Numeric. Mean of a normal or half-normal prior.

- sd:

  Numeric. Standard deviation of a normal or half-normal prior.

- location:

  Numeric. Location of a Cauchy or Student-t prior.

- scale:

  Numeric. Scale of a Cauchy, Student-t, or scaled inverse chi-squared
  prior.

- df:

  Numeric. Degrees of freedom for a Student-t or scaled inverse
  chi-squared prior.

- rate:

  Numeric. Rate parameter of an exponential prior.

- lower, upper:

  Numeric. Endpoints of a uniform prior.

- alpha:

  Numeric (vector for Dirichlet, scalar for Beta). Shape parameter(s).

- beta:

  Numeric. Second shape parameter of a Beta prior.

- eta:

  Numeric. Shape parameter of an LKJ correlation prior. `eta = 1` is
  uniform; `eta > 1` concentrates toward the identity; `eta < 1`
  concentrates toward perfect correlation.

## Value

A `bayesma_prior` object: a list with the family name and
hyperparameters.
