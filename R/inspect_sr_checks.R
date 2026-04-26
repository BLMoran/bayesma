# ============================================================================
# INSPECT-SR: Core Statistical Check Functions (Frequentist + Bayesian)
# ============================================================================
# Automated checks for INSPECT-SR Domain 4 (Inspecting results in the study)
# Based on: Wilkinson et al. (2025) INSPECT-SR tool
#
# Each check has a frequentist version (binary pass/fail) and a Bayesian
# version (returns a Bayes factor quantifying evidence for/against
# trustworthiness).
# ============================================================================


# ============================================================================
# GRIM Test
# ============================================================================

#' GRIM Test (Granularity-Related Inconsistency of Means)
#'
#' Tests whether a reported mean is mathematically possible given the sample
#' size and the number of decimal places reported. Applies to data measured
#' on an integer scale (e.g., Likert items, counts).
#'
#' @param mean_value Numeric. The reported mean.
#' @param n Integer. The sample size.
#' @param decimals Integer. Number of decimal places in the reported mean.
#'   If NULL (default), inferred from the reported value.
#' @param tolerance Numeric. Rounding tolerance for comparison (default 1e-6).
#'
#' @return A list with components:
#'   \describe{
#'     \item{consistent}{Logical. TRUE if the mean is GRIM-consistent.}
#'     \item{mean_value}{The tested mean.}
#'     \item{n}{The sample size.}
#'     \item{decimals}{Number of decimal places used.}
#'   }
#'
#' @details
#' The GRIM test (Brown & Heathers, 2017) checks whether a reported mean of
#' integer data is consistent with the reported sample size. For example, with
#' N = 20, a mean must be a multiple of 1/20 = 0.05. A reported mean of 3.47
#' would be impossible.
#'
#' This implements INSPECT-SR check 4.8.
#'
#' @references
#' Brown NJL, Heathers JAJ (2017). The GRIM test: A simple technique detects
#' numerous anomalies in the reporting of results in psychology.
#' *Social Psychological and Personality Science*, 8(4), 363-369.
#'
#' @examples
#' # Possible mean: 52/20 = 2.60
#' grim_test(2.60, n = 20)
#'
#' # Impossible mean: no integer sum / 20 = 2.47
#' grim_test(2.47, n = 20)
#'
#' @export
grim_test <- function(mean_value, n, decimals = NULL, tolerance = 1e-6) {

  if (!is.numeric(mean_value) || !is.numeric(n)) {
    rlang::abort("Both `mean_value` and `n` must be numeric.")
  }
  if (n <= 0 || n != round(n)) {
    rlang::abort("`n` must be a positive integer.")
  }

  if (is.null(decimals)) {
    decimals <- count_decimals(mean_value)
  }

  half_bin <- 0.5 * 10^(-decimals)
  k_lower <- ceiling((mean_value - half_bin) * n)
  k_upper <- floor((mean_value + half_bin) * n)
  consistent <- (k_lower <= k_upper) && (k_lower >= 0)

  list(
    consistent = consistent,
    mean_value = mean_value,
    n = n,
    decimals = decimals
  )
}


#' Bayesian GRIM Test
#'
#' Computes a Bayes factor quantifying evidence that a reported mean of
#' integer-scale data is inconsistent with the reported sample size, accounting
#' for rounding uncertainty.
#'
#' @inheritParams grim_test
#' @param max_items Integer. Maximum plausible value on the integer scale
#'   (default 10). Used to define the range of possible integer sums.
#' @param n_tolerance Integer. How many values around the reported N to
#'   consider as plausible (default 0, exact N only).
#'
#' @return A list with components:
#'   \describe{
#'     \item{bf_inconsistent}{Numeric. Bayes factor in favour of
#'       inconsistency (fabrication).}
#'     \item{posterior_prob_inconsistent}{Numeric. Posterior probability of
#'       inconsistency assuming equal prior odds.}
#'     \item{consistent_at_n}{Logical. Classical GRIM result at exact N.}
#'     \item{consistent_nearby}{Logical. GRIM-consistent at any nearby N.}
#'     \item{interpretation}{Character. Evidence strength label.}
#'   }
#'
#' @details
#' Under H0 (genuine data), the mean must equal k/n for some integer k,
#' rounded to the reported decimal places. Under H1 (fabricated data), the
#' mean is drawn uniformly from the plausible range.
#' BF_10 = P(data | H1) / P(data | H0).
#'
#' @examples
#' # Consistent mean: evidence for genuine data
#' bayes_grim_test(2.60, n = 20)
#'
#' # Inconsistent mean: evidence for fabrication
#' bayes_grim_test(2.47, n = 20)
#'
#' @export
bayes_grim_test <- function(mean_value, n, decimals = NULL, max_items = 10,
                            n_tolerance = 0) {

  if (!is.numeric(mean_value) || !is.numeric(n)) {
    rlang::abort("Both `mean_value` and `n` must be numeric.")
  }
  if (n <= 0 || n != round(n)) {
    rlang::abort("`n` must be a positive integer.")
  }

  if (is.null(decimals)) {
    decimals <- count_decimals(mean_value)
  }

  rounding_unit <- 0.5 * 10^(-decimals)
  range_width <- max_items

  p_h1 <- (2 * rounding_unit) / range_width

  n_values <- seq(max(1, n - n_tolerance), n + n_tolerance)
  n_weights <- stats::dnorm(n_values, mean = n, sd = max(1, n_tolerance / 2))
  n_weights <- n_weights / sum(n_weights)

  p_h0 <- 0
  for (idx in seq_along(n_values)) {
    ni <- n_values[idx]
    max_sum <- ni * max_items
    k_lower <- ceiling((mean_value - rounding_unit) * ni)
    k_upper <- floor((mean_value + rounding_unit) * ni)
    k_lower <- max(0L, k_lower)
    k_upper <- min(max_sum, k_upper)
    matches <- max(0L, k_upper - k_lower + 1L)
    p_h0 <- p_h0 + n_weights[idx] * matches / (max_sum + 1)
  }

  if (p_h0 == 0) {
    bf_inconsistent <- Inf
    posterior_prob <- 1.0
  } else if (p_h1 == 0) {
    bf_inconsistent <- 0
    posterior_prob <- 0.0
  } else {
    bf_inconsistent <- p_h1 / p_h0
    posterior_prob <- bf_inconsistent / (1 + bf_inconsistent)
  }

  classical <- grim_test(mean_value, n, decimals = decimals)
  nearby_consistent <- any(vapply(n_values, function(ni) {
    grim_test(mean_value, ni, decimals = decimals)$consistent
  }, logical(1)))

  interpretation <- if (bf_inconsistent > 10) {
    "strong_inconsistency"
  } else if (bf_inconsistent > 3) {
    "moderate_inconsistency"
  } else if (bf_inconsistent > 1) {
    "weak_evidence"
  } else {
    "consistent"
  }

  list(
    bf_inconsistent = bf_inconsistent,
    posterior_prob_inconsistent = posterior_prob,
    consistent_at_n = classical$consistent,
    consistent_nearby = nearby_consistent,
    interpretation = interpretation,
    mean_value = mean_value,
    n = n,
    decimals = decimals
  )
}


# ============================================================================
# P-value verification
# ============================================================================

#' Verify a Reported P-Value
#'
#' Recalculates a p-value from a reported test statistic and degrees of freedom,
#' then compares it to the reported p-value.
#'
#' @param test_type Character. One of \code{"t"}, \code{"z"}, \code{"chi_sq"},
#'   \code{"f"}.
#' @param statistic Numeric. The reported test statistic.
#' @param df Numeric. Degrees of freedom. For F-tests, a vector of length 2.
#'   Not required for z-tests.
#' @param reported_p Numeric. The reported p-value.
#' @param alternative Character. One of \code{"two.sided"} (default),
#'   \code{"less"}, \code{"greater"}.
#' @param tolerance Numeric. Acceptable absolute difference (default 0.01).
#'
#' @return A list with components:
#'   \describe{
#'     \item{consistent}{Logical. TRUE if p-values match within tolerance.}
#'     \item{reported_p}{The reported p-value.}
#'     \item{recalculated_p}{The recalculated p-value.}
#'     \item{difference}{Absolute difference.}
#'     \item{test_type}{The test type used.}
#'     \item{statistic}{The test statistic.}
#'   }
#'
#' @details
#' Implements INSPECT-SR check 4.9.
#'
#' @examples
#' # Consistent: chi-squared = 3.84, df = 1, p = 0.05
#' verify_pvalue("chi_sq", statistic = 3.84, df = 1, reported_p = 0.05)
#'
#' @export
verify_pvalue <- function(test_type = c("t", "z", "chi_sq", "f"),
                          statistic,
                          df = NULL,
                          reported_p,
                          alternative = "two.sided",
                          tolerance = 0.01) {

  test_type <- rlang::arg_match(test_type)

  if (!is.numeric(statistic) || !is.numeric(reported_p)) {
    rlang::abort("`statistic` and `reported_p` must be numeric.")
  }

  recalculated_p <- switch(test_type,
                           "t" = {
                             if (is.null(df)) rlang::abort("`df` is required for t-tests.")
                             switch(alternative,
                                    "two.sided" = 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE),
                                    "less" = stats::pt(statistic, df = df, lower.tail = TRUE),
                                    "greater" = stats::pt(statistic, df = df, lower.tail = FALSE)
                             )
                           },
                           "z" = {
                             switch(alternative,
                                    "two.sided" = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
                                    "less" = stats::pnorm(statistic, lower.tail = TRUE),
                                    "greater" = stats::pnorm(statistic, lower.tail = FALSE)
                             )
                           },
                           "chi_sq" = {
                             if (is.null(df)) rlang::abort("`df` is required for chi-squared tests.")
                             stats::pchisq(statistic, df = df, lower.tail = FALSE)
                           },
                           "f" = {
                             if (is.null(df) || length(df) != 2) {
                               rlang::abort("`df` must be length 2 (df1, df2) for F-tests.")
                             }
                             stats::pf(statistic, df1 = df[1], df2 = df[2], lower.tail = FALSE)
                           }
  )

  difference <- abs(reported_p - recalculated_p)

  list(
    consistent = difference <= tolerance,
    reported_p = reported_p,
    recalculated_p = recalculated_p,
    difference = difference,
    test_type = test_type,
    statistic = statistic
  )
}


#' Bayesian P-Value Verification
#'
#' Computes a Bayes factor quantifying evidence that a reported p-value is
#' inconsistent with the reported test statistic.
#'
#' @inheritParams verify_pvalue
#' @param rounding_sd Numeric. SD of the rounding error model (default 0.005).
#' @param fabrication_sd Numeric. SD under the fabrication model (default 0.15).
#'
#' @return A list with components:
#'   \describe{
#'     \item{bf_inconsistent}{Numeric. Bayes factor for inconsistency.}
#'     \item{posterior_prob_inconsistent}{Numeric. Posterior probability.}
#'     \item{recalculated_p}{Numeric. Recalculated p-value.}
#'     \item{discrepancy}{Numeric. Absolute difference.}
#'     \item{interpretation}{Character. Evidence strength label.}
#'   }
#'
#' @details
#' Discrepancy modelled as N(0, rounding_sd^2) under H0 (honest rounding)
#' and N(0, fabrication_sd^2) under H1 (fabrication/error).
#'
#' @examples
#' bayes_verify_pvalue("chi_sq", statistic = 3.84, df = 1, reported_p = 0.05)
#'
#' @export
bayes_verify_pvalue <- function(test_type = c("t", "z", "chi_sq", "f"),
                                statistic,
                                df = NULL,
                                reported_p,
                                alternative = "two.sided",
                                rounding_sd = 0.005,
                                fabrication_sd = 0.15) {

  test_type <- rlang::arg_match(test_type)

  freq_result <- verify_pvalue(test_type, statistic, df, reported_p, alternative,
                               tolerance = Inf)
  recalculated_p <- freq_result$recalculated_p
  discrepancy <- abs(reported_p - recalculated_p)

  lik_h0 <- stats::dnorm(discrepancy, mean = 0, sd = rounding_sd)
  lik_h1 <- stats::dnorm(discrepancy, mean = 0, sd = fabrication_sd)
  bf_inconsistent <- lik_h1 / lik_h0
  posterior_prob <- bf_inconsistent / (1 + bf_inconsistent)

  interpretation <- if (bf_inconsistent > 10) {
    "strong_inconsistency"
  } else if (bf_inconsistent > 3) {
    "moderate_inconsistency"
  } else if (bf_inconsistent > 1) {
    "weak_evidence"
  } else {
    "consistent"
  }

  list(
    bf_inconsistent = bf_inconsistent,
    posterior_prob_inconsistent = posterior_prob,
    recalculated_p = recalculated_p,
    discrepancy = discrepancy,
    interpretation = interpretation,
    test_type = test_type,
    statistic = statistic
  )
}


# ============================================================================
# Carlisle's test
# ============================================================================

#' Carlisle's Test for Baseline Balance
#'
#' Tests whether the distribution of p-values for baseline comparisons is
#' consistent with genuine randomisation.
#'
#' @param p_values Numeric vector. P-values from baseline comparisons.
#' @param method Character. \code{"fisher"} (default) or \code{"ks"}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{too_similar}{Logical. Suspiciously well-balanced.}
#'     \item{too_different}{Logical. Suspiciously imbalanced.}
#'     \item{combined_p}{Combined p-value.}
#'     \item{n_comparisons}{Number of comparisons.}
#'     \item{method}{Method used.}
#'     \item{interpretation}{\code{"plausible"}, \code{"too_similar"}, or
#'       \code{"too_different"}.}
#'   }
#'
#' @details
#' Implements INSPECT-SR check 4.3. Under genuine randomisation, baseline
#' p-values should be approximately uniform. Fabricated trials often show
#' implausibly well-matched groups (p-values near 1).
#'
#' @references
#' Carlisle JB (2017). Data fabrication and other reasons for non-random
#' sampling in 5087 randomised, controlled trials in anaesthetic and general
#' medical journals. *Anaesthesia*, 72(8), 944-952.
#'
#' @examples
#' carlisle_test(c(0.45, 0.12, 0.78, 0.33, 0.91))
#' carlisle_test(c(0.92, 0.88, 0.95, 0.91, 0.87))
#'
#' @export
carlisle_test <- function(p_values, method = c("fisher", "ks")) {

  method <- rlang::arg_match(method)

  if (!is.numeric(p_values)) {
    rlang::abort("`p_values` must be a numeric vector.")
  }

  na_count <- sum(is.na(p_values))
  if (na_count > 0) {
    rlang::warn(paste0("Removed ", na_count, " NA p-value(s)."))
    p_values <- p_values[!is.na(p_values)]
  }

  if (length(p_values) < 2) {
    rlang::abort("At least 2 p-values are required for Carlisle's test.")
  }

  p_values <- pmax(pmin(p_values, 1 - 1e-10), 1e-10)
  k <- length(p_values)

  if (method == "fisher") {
    fisher_stat <- -2 * sum(log(p_values))
    p_too_different <- stats::pchisq(fisher_stat, df = 2 * k, lower.tail = FALSE)
    p_too_similar <- stats::pchisq(fisher_stat, df = 2 * k, lower.tail = TRUE)
    combined_p <- min(p_too_similar, p_too_different)
    too_similar <- p_too_similar < 0.05
    too_different <- p_too_different < 0.05
  } else {
    ks_result <- stats::ks.test(p_values, "punif")
    combined_p <- ks_result$p.value
    too_similar <- combined_p < 0.05 && stats::median(p_values) > 0.5
    too_different <- combined_p < 0.05 && stats::median(p_values) <= 0.5
  }

  interpretation <- if (too_similar) "too_similar"
  else if (too_different) "too_different"
  else "plausible"

  list(
    too_similar = too_similar,
    too_different = too_different,
    combined_p = combined_p,
    n_comparisons = k,
    method = method,
    interpretation = interpretation
  )
}


#' Bayesian Carlisle's Test for Baseline Balance
#'
#' Computes a Bayes factor comparing genuine randomisation (uniform p-values)
#' to fabrication (non-uniform).
#'
#' @inheritParams carlisle_test
#' @param prior_a,prior_b Numeric. Beta prior shape parameters (default 1, 1).
#' @param n_grid Integer. Grid size for numerical integration (default 200).
#'
#' @return A list with components:
#'   \describe{
#'     \item{bf_too_similar}{BF for p-values biased towards 1.}
#'     \item{bf_too_different}{BF for p-values biased towards 0.}
#'     \item{bf_nonuniform}{Overall BF for non-uniformity.}
#'     \item{posterior_prob_fabrication}{Posterior probability.}
#'     \item{posterior_mean_p}{Mean of observed p-values.}
#'     \item{n_comparisons}{Number of comparisons.}
#'     \item{interpretation}{Evidence description.}
#'   }
#'
#' @examples
#' bayes_carlisle_test(c(0.45, 0.12, 0.78, 0.33, 0.91))
#' bayes_carlisle_test(c(0.93, 0.81, 0.95, 0.94, 0.85, 0.95))
#'
#' @export
bayes_carlisle_test <- function(p_values,
                                prior_a = 1,
                                prior_b = 1,
                                n_grid = 200) {

  if (!is.numeric(p_values)) {
    rlang::abort("`p_values` must be a numeric vector.")
  }

  na_count <- sum(is.na(p_values))
  if (na_count > 0) {
    rlang::warn(paste0("Removed ", na_count, " NA p-value(s)."))
    p_values <- p_values[!is.na(p_values)]
  }

  if (length(p_values) < 2) {
    rlang::abort("At least 2 p-values are required.")
  }

  p_values <- pmax(pmin(p_values, 1 - 1e-10), 1e-10)
  k <- length(p_values)

  log_lik_h0 <- 0

  alpha_grid <- exp(seq(log(0.1), log(20), length.out = n_grid))
  beta_grid <- exp(seq(log(0.1), log(20), length.out = n_grid))

  log_marginals <- matrix(NA_real_, nrow = n_grid, ncol = n_grid)
  sum_log_p <- sum(log(p_values))
  sum_log_1mp <- sum(log(1 - p_values))

  for (i in seq_along(alpha_grid)) {
    for (j in seq_along(beta_grid)) {
      a <- alpha_grid[i]
      b <- beta_grid[j]
      log_marginals[i, j] <- (a - 1) * sum_log_p + (b - 1) * sum_log_1mp -
        k * lbeta(a, b)
    }
  }

  max_lm <- max(log_marginals)
  log_marginal_h1 <- max_lm + log(sum(exp(log_marginals - max_lm))) -
    log(n_grid^2)

  bf_nonuniform <- exp(log_marginal_h1 - log_lik_h0)

  # Directional BFs
  mask_similar <- outer(alpha_grid, beta_grid, function(a, b) a > b)
  bf_too_similar <- if (any(mask_similar)) {
    lm_sim <- log_marginals
    lm_sim[!mask_similar] <- -Inf
    max_ls <- max(lm_sim[mask_similar])
    exp(max_ls + log(sum(exp(lm_sim[mask_similar] - max_ls))) -
          log(sum(mask_similar)) - log_lik_h0)
  } else NA_real_

  mask_different <- outer(alpha_grid, beta_grid, function(a, b) b > a)
  bf_too_different <- if (any(mask_different)) {
    lm_diff <- log_marginals
    lm_diff[!mask_different] <- -Inf
    max_ld <- max(lm_diff[mask_different])
    exp(max_ld + log(sum(exp(lm_diff[mask_different] - max_ld))) -
          log(sum(mask_different)) - log_lik_h0)
  } else NA_real_

  posterior_prob <- bf_nonuniform / (1 + bf_nonuniform)

  interpretation <- if (bf_nonuniform > 10) {
    if (!is.na(bf_too_similar) && bf_too_similar > bf_too_different) {
      "strong_evidence_too_similar"
    } else {
      "strong_evidence_too_different"
    }
  } else if (bf_nonuniform > 3) {
    "moderate_evidence_nonuniform"
  } else if (bf_nonuniform > 1) {
    "weak_evidence_nonuniform"
  } else {
    "consistent_with_randomisation"
  }

  list(
    bf_too_similar = bf_too_similar,
    bf_too_different = bf_too_different,
    bf_nonuniform = bf_nonuniform,
    posterior_prob_fabrication = posterior_prob,
    posterior_mean_p = mean(p_values),
    n_comparisons = k,
    interpretation = interpretation
  )
}


# ============================================================================
# Participant number consistency
# ============================================================================

#' Check Participant Number Consistency
#'
#' Verifies that reported participant numbers are internally consistent:
#' randomised = analysed + lost to follow-up.
#'
#' @param n_randomised_int,n_randomised_ctrl Integer. Per-arm randomised counts.
#' @param n_analysed_int,n_analysed_ctrl Integer or NULL. Per-arm analysed.
#' @param n_lost_int,n_lost_ctrl Integer or NULL. Per-arm lost to follow-up.
#' @param n_randomised_total Integer or NULL. Total randomised.
#'
#' @return A list with components:
#'   \describe{
#'     \item{consistent}{Logical. TRUE if all checks pass.}
#'     \item{checks}{Data frame of individual checks.}
#'     \item{n_checks}{Number of checks performed.}
#'     \item{n_failed}{Number of failed checks.}
#'   }
#'
#' @details
#' Implements INSPECT-SR check 4.6. Missing values (NULL) are skipped.
#'
#' @examples
#' check_n_consistency(
#'   n_randomised_int = 100, n_randomised_ctrl = 100,
#'   n_analysed_int = 95, n_analysed_ctrl = 92,
#'   n_lost_int = 5, n_lost_ctrl = 8,
#'   n_randomised_total = 200
#' )
#'
#' @export
check_n_consistency <- function(n_randomised_int,
                                n_randomised_ctrl,
                                n_analysed_int = NULL,
                                n_analysed_ctrl = NULL,
                                n_lost_int = NULL,
                                n_lost_ctrl = NULL,
                                n_randomised_total = NULL) {

  checks <- data.frame(
    check = character(), expected = numeric(),
    observed = numeric(), pass = logical(),
    stringsAsFactors = FALSE
  )

  add_check <- function(check, expected, observed) {
    checks <<- rbind(checks, data.frame(
      check = check, expected = expected,
      observed = observed, pass = expected == observed,
      stringsAsFactors = FALSE
    ))
  }

  if (!is.null(n_randomised_total)) {
    add_check("Total randomised = Intervention + Control",
              n_randomised_int + n_randomised_ctrl, n_randomised_total)
  }
  if (!is.null(n_analysed_int) && !is.null(n_lost_int)) {
    add_check("Intervention: Randomised = Analysed + Lost",
              n_randomised_int, n_analysed_int + n_lost_int)
  }
  if (!is.null(n_analysed_ctrl) && !is.null(n_lost_ctrl)) {
    add_check("Control: Randomised = Analysed + Lost",
              n_randomised_ctrl, n_analysed_ctrl + n_lost_ctrl)
  }
  if (!is.null(n_lost_int)) {
    checks <- rbind(checks, data.frame(
      check = "Intervention: Lost <= Randomised",
      expected = n_randomised_int, observed = n_lost_int,
      pass = n_lost_int <= n_randomised_int, stringsAsFactors = FALSE
    ))
  }
  if (!is.null(n_lost_ctrl)) {
    checks <- rbind(checks, data.frame(
      check = "Control: Lost <= Randomised",
      expected = n_randomised_ctrl, observed = n_lost_ctrl,
      pass = n_lost_ctrl <= n_randomised_ctrl, stringsAsFactors = FALSE
    ))
  }

  list(
    consistent = nrow(checks) == 0 || all(checks$pass),
    checks = checks,
    n_checks = nrow(checks),
    n_failed = sum(!checks$pass)
  )
}


#' Check Internal Consistency of Summary Statistics
#'
#' Verifies relationships between reported summary statistics (CI symmetry,
#' SE vs CI width).
#'
#' @param estimate Numeric. Point estimate.
#' @param ci_lower,ci_upper Numeric or NULL. CI bounds.
#' @param se Numeric or NULL. Standard error.
#' @param log_scale Logical. Check on log scale (default FALSE).
#' @param ci_level Numeric. Confidence level (default 0.95).
#' @param tolerance Numeric. Tolerance (default 0.1).
#'
#' @return A list with `consistent` and `checks`.
#' @export
check_statistics_consistency <- function(estimate,
                                         ci_lower = NULL,
                                         ci_upper = NULL,
                                         se = NULL,
                                         log_scale = FALSE,
                                         ci_level = 0.95,
                                         tolerance = 0.1) {

  checks <- data.frame(
    check = character(), expected = numeric(),
    observed = numeric(), pass = logical(),
    stringsAsFactors = FALSE
  )

  z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)

  if (log_scale) {
    est <- log(estimate)
    lower <- if (!is.null(ci_lower)) log(ci_lower) else NULL
    upper <- if (!is.null(ci_upper)) log(ci_upper) else NULL
  } else {
    est <- estimate
    lower <- ci_lower
    upper <- ci_upper
  }

  if (!is.null(lower) && !is.null(upper)) {
    ci_midpoint <- (lower + upper) / 2
    checks <- rbind(checks, data.frame(
      check = "CI midpoint matches estimate",
      expected = est, observed = ci_midpoint,
      pass = abs(est - ci_midpoint) <= tolerance,
      stringsAsFactors = FALSE
    ))
  }

  if (!is.null(se) && !is.null(lower) && !is.null(upper)) {
    ci_width <- upper - lower
    expected_width <- 2 * z_crit * se
    checks <- rbind(checks, data.frame(
      check = "CI width consistent with SE",
      expected = expected_width, observed = ci_width,
      pass = abs(ci_width - expected_width) <= tolerance,
      stringsAsFactors = FALSE
    ))
  }

  list(
    consistent = nrow(checks) == 0 || all(checks$pass),
    checks = checks
  )
}


# ============================================================================
# Internal helpers
# ============================================================================

#' @noRd
count_decimals <- function(x) {
  x_str <- format(x, scientific = FALSE)
  if (!grepl("\\.", x_str)) return(0L)
  nchar(sub(".*\\.", "", x_str))
}
