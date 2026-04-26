# Prior distribution constructors

#' Prior distribution constructors
#'
#' Constructors for prior distributions used by [bayesma()], [robma()],
#' [meta_reg()], and related fitting functions. Each returns an object of
#' class `bayesma_prior` that stores the family and its hyperparameters.
#'
#' @param mean Numeric. Mean of a normal or half-normal prior.
#' @param sd Numeric. Standard deviation of a normal or half-normal prior.
#' @param location Numeric. Location of a Cauchy or Student-t prior.
#' @param scale Numeric. Scale of a Cauchy, Student-t, or scaled inverse
#'   chi-squared prior.
#' @param df Numeric. Degrees of freedom for a Student-t or scaled inverse
#'   chi-squared prior.
#' @param rate Numeric. Rate parameter of an exponential prior.
#' @param lower,upper Numeric. Endpoints of a uniform prior.
#' @param alpha Numeric (vector for Dirichlet, scalar for Beta).
#'   Shape parameter(s).
#' @param beta Numeric. Second shape parameter of a Beta prior.
#' @param eta Numeric. Shape parameter of an LKJ correlation prior.
#'   `eta = 1` is uniform; `eta > 1` concentrates toward the identity;
#'   `eta < 1` concentrates toward perfect correlation.
#'
#' @return A `bayesma_prior` object: a list with the family name and
#'   hyperparameters.
#'
#' @name priors
NULL

#' @rdname priors
#' @export
normal <- function(mean, sd) {
  structure(list(family = "normal", mean = mean, sd = sd),
            class = "bayesma_prior")
}

#' @rdname priors
#' @export
half_normal <- function(mean, sd) {
  structure(list(family = "half_normal", mean = mean, sd = sd),
            class = "bayesma_prior")
}

#' @rdname priors
#' @export
half_cauchy <- function(location, scale) {
  structure(list(family = "half_cauchy", location = location, scale = scale),
            class = "bayesma_prior")
}

#' @rdname priors
#' @export
half_student_t <- function(df, location, scale) {
  structure(list(family = "half_student_t", df = df,
                 location = location, scale = scale),
            class = "bayesma_prior")
}

#' @rdname priors
#' @export
exponential <- function(rate) {
  structure(list(family = "exponential", rate = rate),
            class = "bayesma_prior")
}

#' @rdname priors
#' @export
uniform <- function(lower, upper) {
  structure(list(family = "uniform", lower = lower, upper = upper),
            class = "bayesma_prior")
}

#' @rdname priors
#' @export
dirichlet <- function(alpha) {
  structure(list(family = "dirichlet", alpha = alpha),
            class = "bayesma_prior")
}

#' @rdname priors
#' @export
beta <- function(alpha, beta) {
  structure(list(family = "beta", alpha = alpha, beta = beta),
            class = "bayesma_prior")
}

#' @rdname priors
#' @export
scaled_inv_chi_sq <- function(df, scale) {
  structure(list(family = "scaled_inv_chi_sq", df = df, scale = scale),
            class = "bayesma_prior")
}

#' @rdname priors
#' @export
lkj <- function(eta) {
  structure(list(family = "lkj", eta = eta),
            class = "bayesma_prior")
}


#' @noRd
format.bayesma_prior <- function(x, ...) {
  switch(x$family,
         normal             = glue::glue("N({x$mean}, {x$sd})"),
         half_normal        = glue::glue("HN({x$mean}, {x$sd})"),
         half_cauchy        = glue::glue("HC({x$location}, {x$scale})"),
         half_student_t     = glue::glue("t({x$df}, {x$location}, {x$scale})"),
         exponential        = glue::glue("Exp({x$rate})"),
         uniform            = glue::glue("U({x$lower}, {x$upper})"),
         dirichlet          = glue::glue("Dir({paste(x$alpha, collapse = ', ')})"),
         beta               = glue::glue("Beta({x$alpha}, {x$beta})"),
         scaled_inv_chi_sq  = glue::glue("Scaled-Inv-Chi2({x$df}, {x$scale})"),
         as.character(x$family)
  )
}

#' @noRd
print.bayesma_prior <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}
