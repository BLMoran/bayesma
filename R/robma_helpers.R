# RoBMA Helper Functions

#' RoBMA bias-prior constructors
#'
#' Constructors for the publication-bias priors used by [robma()] and
#' [robma_sensitivity()]. Each returns a `robma_bias_prior` object.
#'
#' @param type Bias-prior family. One of `"weight_function"`, `"pet"`,
#'   `"peese"`, `"copas"`, `"jung"`, `"none"`.
#' @param parameters Named list of family-specific parameters.
#' @param prior_weight Numeric. Prior model weight in the bias-prior mixture.
#' @param steps Numeric vector of p-value cutpoints for a step weight function.
#' @param alpha Numeric vector of Dirichlet concentration parameters; defaults
#'   to `rep(1, length(steps) + 1)`.
#' @param sided One of `"one"` or `"two"`. Determines whether the weight
#'   function is one- or two-sided.
#' @param distribution Distribution for the PET / PEESE slope prior. Currently
#'   `"cauchy"`.
#' @param location,scale Location and scale of the PET / PEESE slope prior.
#' @param ... Additional fields stored on the bias-prior object.
#'
#' @return A `robma_bias_prior` object.
#'
#' @name prior_bias
#' @export
prior_bias <- function(
    type = c("weight_function", "pet", "peese", "copas", "jung", "none"),
    parameters = list(),
    prior_weight = 1,
    ...
) {
  type <- rlang::arg_match(type)
  structure(
    list(type = type, parameters = parameters,
         prior_weight = prior_weight, ...),
    class = "robma_bias_prior"
  )
}

#' @rdname prior_bias
#' @export
prior_weight_function <- function(steps = c(0.025, 0.05),
                                  alpha = NULL, prior_weight = 1,
                                  sided = "one") {
  if (is.null(alpha)) alpha <- rep(1, length(steps) + 1)
  prior_bias("weight_function",
             parameters = list(steps = steps, alpha = alpha,
                               sided = sided),
             prior_weight = prior_weight)
}

#' @rdname prior_bias
#' @export
prior_pet <- function(distribution = "cauchy", location = 0, scale = 1,
                      prior_weight = 1) {
  prior_bias("pet",
             parameters = list(distribution = distribution,
                               location = location, scale = scale),
             prior_weight = prior_weight)
}

#' @rdname prior_bias
#' @export
prior_peese <- function(distribution = "cauchy", location = 0, scale = 5,
                        prior_weight = 1) {
  prior_bias("peese",
             parameters = list(distribution = distribution,
                               location = location, scale = scale),
             prior_weight = prior_weight)
}

#' @rdname prior_bias
#' @export
prior_copas <- function(prior_weight = 1) {
  prior_bias("copas", prior_weight = prior_weight)
}

#' @rdname prior_bias
#' @export
prior_jung <- function(prior_weight = 1) {
  prior_bias("jung", prior_weight = prior_weight)
}

#' @rdname prior_bias
#' @export
prior_no_bias <- function(prior_weight = 1) {
  prior_bias("none", prior_weight = prior_weight)
}


#' Default RoBMA prior set
#'
#' Returns the default list of priors for one of the three RoBMA dimensions
#' (effect, heterogeneity, bias). Used internally by [robma()] when the user
#' does not pass explicit priors.
#'
#' @param dimension One of `"effect"`, `"heterogeneity"`, `"bias"`.
#' @param null Logical. If `TRUE`, return the null-model priors (point at zero
#'   for effect / heterogeneity, no-bias for bias).
#' @param rescale Numeric. Multiplier on default scales.
#'
#' @return A list of `bayesma_prior` or `robma_bias_prior` objects.
#'
#' @export
robma_default_priors <- function(
    dimension = c("effect", "heterogeneity", "bias"),
    null = FALSE, rescale = 1
) {
  dimension <- rlang::arg_match(dimension)
  if (dimension == "effect") {
    if (null) list(list(distribution = "point", location = 0))
    else      list(normal(0, 1 * rescale))
  } else if (dimension == "heterogeneity") {
    if (null) list(list(distribution = "point", location = 0))
    else      list(half_student_t(3, 0, 2.5 * rescale))
  } else if (dimension == "bias") {
    if (null) {
      list(prior_no_bias())
    } else {
      list(
        prior_weight_function(steps = c(0.05),
                              prior_weight = 1/12, sided = "two"),
        prior_weight_function(steps = c(0.05, 0.10),
                              prior_weight = 1/12, sided = "two"),
        prior_weight_function(steps = c(0.05),
                              prior_weight = 1/12, sided = "one"),
        prior_weight_function(steps = c(0.025, 0.05),
                              prior_weight = 1/12, sided = "one"),
        prior_weight_function(steps = c(0.05, 0.5),
                              prior_weight = 1/12, sided = "one"),
        prior_weight_function(steps = c(0.025, 0.05, 0.5),
                              prior_weight = 1/12, sided = "one"),
        prior_pet(prior_weight = 1/4),
        prior_peese(prior_weight = 1/4)
      )
    }
  }
}


# Model grid construction

#' @noRd
build_model_grid <- function(
    priors_effect, priors_effect_null,
    priors_heterogeneity, priors_heterogeneity_null,
    priors_bias, priors_bias_null
) {
  models <- list()
  effect_entries <- c(
    lapply(priors_effect, function(p) list(prior = p, is_null = FALSE)),
    lapply(priors_effect_null, function(p) list(prior = p, is_null = TRUE))
  )
  hetero_entries <- c(
    lapply(priors_heterogeneity, function(p) list(prior = p, is_null = FALSE)),
    lapply(priors_heterogeneity_null, function(p) list(prior = p, is_null = TRUE))
  )
  bias_entries <- c(
    lapply(priors_bias, function(p) list(prior = p, is_null = FALSE)),
    lapply(priors_bias_null, function(p) list(prior = p, is_null = TRUE))
  )

  for (e in effect_entries) {
    for (h in hetero_entries) {
      for (b in bias_entries) {
        ew <- e$prior$prior_weight %||% 1
        hw <- h$prior$prior_weight %||% 1
        bw <- b$prior$prior_weight %||% 1
        e_lbl <- if (e$is_null) "H0" else "H1"
        h_lbl <- if (h$is_null) "FE" else "RE"
        if (b$is_null) {
          b_lbl <- "no bias"
        } else {
          b_type <- b$prior$type %||% "bias"
          if (b_type == "weight_function") {
            steps_str <- paste(b$prior$parameters$steps, collapse = ",")
            sided <- b$prior$parameters$sided %||% "one"
            side_tag <- if (sided == "two") "2s" else "1s"
            b_lbl <- paste0("wf-", side_tag, "(", steps_str, ")")
          } else {
            b_lbl <- toupper(b_type)
          }
        }
        label <- paste0(h_lbl, " / ", e_lbl, " / ", b_lbl)
        models <- c(models, list(list(
          effect_prior = e$prior, hetero_prior = h$prior,
          bias_prior = b$prior,
          is_effect_null = e$is_null, is_hetero_null = h$is_null,
          is_bias_null = b$is_null,
          prior_weight = ew * hw * bw, label = label
        )))
      }
    }
  }

  cli::cli_alert_info(paste0(
    "Model grid: ", length(models), " models from ",
    length(effect_entries), " effect x ",
    length(hetero_entries), " heterogeneity x ",
    length(bias_entries), " bias priors"))
  models
}


#' Compute null range probabilities
#'
#' @noRd
compute_null_range_probs <- function(mu_draws, null_range = NULL,
                                     effect_label = "log_or") {
  is_ratio <- effect_label %in% c("log_or", "log_rr")

  if (is.null(null_range)) {
    # Point null at exactly zero
    p_neg  <- mean(mu_draws < 0)
    p_null <- mean(mu_draws == 0)
    p_pos  <- mean(mu_draws > 0)
    nr_log <- c(0, 0)
    nr_nat <- if (is_ratio) c(1, 1) else c(0, 0)
  } else {
    stopifnot(length(null_range) == 2, null_range[1] <= null_range[2])

    if (is_ratio) {
      # Auto-detect scale:
      # If both values > 0 and range straddles 1 -> natural scale
      # If range straddles 0 -> log scale
      if (all(null_range > 0) && null_range[1] <= 1 && null_range[2] >= 1) {
        # Natural scale (e.g., OR: c(0.9, 1.1))
        nr_nat <- null_range
        nr_log <- log(null_range)
        cli::cli_alert_info(paste0(
          "Null range interpreted as natural scale: ",
          gsub("log_", "", effect_label), " in [",
          round(nr_nat[1], 3), ", ", round(nr_nat[2], 3), "] -> log scale [",
          round(nr_log[1], 4), ", ", round(nr_log[2], 4), "]"))
      } else {
        # Log scale (e.g., log OR: c(-0.1, 0.1))
        nr_log <- null_range
        nr_nat <- exp(null_range)
      }
    } else {
      # Non-ratio measure: null_range is on the effect scale directly
      nr_log <- null_range
      nr_nat <- null_range
    }

    p_neg  <- mean(mu_draws < nr_log[1])
    p_null <- mean(mu_draws >= nr_log[1] & mu_draws <= nr_log[2])
    p_pos  <- mean(mu_draws > nr_log[2])
  }

  list(
    p_negative = p_neg,
    p_null = p_null,
    p_positive = p_pos,
    null_range = nr_log,
    null_range_natural = nr_nat
  )
}

#' @noRd
compute_inclusion_bf <- function(included, post_probs, prior_weights,
                                 finite_mask = rep(TRUE, length(post_probs))) {
  incl_finite <- included & finite_mask
  excl_finite <- (!included) & finite_mask
  if (sum(incl_finite) == 0 || sum(excl_finite) == 0) return(NA_real_)
  prior_incl <- sum(prior_weights[included]) / sum(prior_weights)
  prior_excl <- 1 - prior_incl
  post_incl  <- sum(post_probs[included])
  post_excl  <- sum(post_probs[!included])
  if (prior_excl == 0 || prior_incl == 0) return(NA_real_)
  if (post_excl == 0) return(Inf)
  if (post_incl == 0) return(0)
  (post_incl / post_excl) / (prior_incl / prior_excl)
}

compute_posterior_probs <- function(log_mls, prior_weights) {
  log_prior   <- log(prior_weights / sum(prior_weights))
  log_unnorm  <- log_mls + log_prior
  finite_mask <- is.finite(log_unnorm)
  post_probs  <- rep(0, length(log_unnorm))
  if (any(finite_mask)) {
    log_max <- max(log_unnorm[finite_mask])
    post_probs[finite_mask] <- exp(log_unnorm[finite_mask] - log_max)
    post_probs <- post_probs / sum(post_probs)
  }
  post_probs
}


# Analytic / numerical log marginal likelihoods for null models

#' @noRd
compute_log_ml_ce_null <- function(y, se) {
  sum(stats::dnorm(y, mean = 0, sd = se, log = TRUE))
}

compute_log_ml_re_null <- function(y, se, tau_prior) {
  tau_log_prior_fn <- build_tau_log_prior(tau_prior)
  log_lik_given_tau <- function(tau_val) {
    sigma <- sqrt(tau_val^2 + se^2)
    sum(stats::dnorm(y, mean = 0, sd = sigma, log = TRUE))
  }
  n_grid <- 4000
  tau_grid <- seq(1e-6, 30, length.out = n_grid)
  dtau <- tau_grid[2] - tau_grid[1]
  log_vals <- vapply(tau_grid, function(tv) {
    log_lik_given_tau(tv) + tau_log_prior_fn(tv)
  }, numeric(1))
  finite_mask <- is.finite(log_vals)
  if (!any(finite_mask)) return(-Inf)
  log_max <- max(log_vals[finite_mask])
  integral <- sum(exp(log_vals[finite_mask] - log_max)) * dtau
  log_max + log(integral)
}

build_tau_log_prior <- function(prior) {
  family <- NULL; params <- NULL
  if (is.list(prior)) {
    family <- prior$family %||% prior$distribution %||%
      prior$dist %||% prior$type %||% NULL
    params <- prior$params %||% prior$parameters %||% NULL
    if (is.null(params)) {
      pn <- setdiff(names(prior),
                    c("family","distribution","dist","type","class","lower","upper","bounds"))
      if (length(pn) > 0) params <- prior[pn]
    }
  }
  if (!is.null(family)) family <- tolower(gsub("[^a-z0-9]", "_", family))

  if (!is.null(family) && grepl("half.*(student|t)", family)) {
    df <- params$df %||% params[[1]] %||% 3
    loc <- params$mu %||% params$loc %||% params[[2]] %||% 0
    sc <- params$sigma %||% params$scale %||% params[[3]] %||% 2.5
    return(function(x) {
      if (x < 0) return(-Inf)
      stats::dt((x - loc)/sc, df = df, log = TRUE) - log(sc) + log(2)
    })
  }
  if (!is.null(family) && grepl("half.*normal", family)) {
    loc <- params$mu %||% params$mean %||% params[[1]] %||% 0
    sc <- params$sigma %||% params$sd %||% params[[2]] %||% 1
    return(function(x) {
      if (x < 0) return(-Inf)
      stats::dnorm(x, mean = loc, sd = sc, log = TRUE) + log(2)
    })
  }
  if (!is.null(family) && grepl("half.*cauchy", family)) {
    loc <- params$mu %||% params$location %||% params[[1]] %||% 0
    sc <- params$sigma %||% params$scale %||% params[[2]] %||% 1
    return(function(x) {
      if (x < 0) return(-Inf)
      stats::dcauchy(x, location = loc, scale = sc, log = TRUE) + log(2)
    })
  }
  if (!is.null(family) && grepl("uniform", family)) {
    lo <- params$lower %||% params$min %||% params[[1]] %||% prior$lower %||% 0
    hi <- params$upper %||% params$max %||% params[[2]] %||% prior$upper %||% 10
    return(function(x) { if (x < lo || x > hi) -Inf else -log(hi - lo) })
  }
  if (!is.null(family) && grepl("exponential", family)) {
    rate <- params$rate %||% params$lambda %||% params[[1]] %||% 1
    return(function(x) { if (x < 0) -Inf else stats::dexp(x, rate, log = TRUE) })
  }
  cli::cli_alert_warning(paste0(
    "Could not identify tau prior family (", family %||% "NULL",
    "), using default half_student_t(3, 0, 2.5) for quadrature"))
  function(x) {
    if (x < 0) return(-Inf)
    stats::dt(x / 2.5, df = 3, log = TRUE) - log(2.5) + log(2)
  }
}


