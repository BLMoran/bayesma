# Internal helpers for bayesma_mv (Burke, Bujkiewicz & Riley 2018).
# Called by bayesma_mv_spec() and bayesma_mv_stan_code().

compute_effect_sizes_mv <- function(mean_ctrl, mean_int, sd_ctrl, sd_int,
                                    n_ctrl, n_int, likelihood) {
  switch(likelihood,
    gaussian = {
      yi  <- mean_int - mean_ctrl
      sp  <- sqrt(((n_ctrl - 1) * sd_ctrl^2 + (n_int - 1) * sd_int^2) /
                    (n_ctrl + n_int - 2))
      sei <- sp * sqrt(1 / n_ctrl + 1 / n_int)
      list(yi = yi, sei = sei, vi = sei^2, measure = "mean_diff")
    }
  )
}

resolve_priors_mv <- function(outcome_labels, mu_prior, tau_prior,
                              rho_between_prior) {
  mu_priors <- if (is.null(mu_prior)) {
    stats::setNames(purrr::map(outcome_labels, \(x) normal(0, 1)), outcome_labels)
  } else if (inherits(mu_prior, "bayesma_prior")) {
    stats::setNames(purrr::map(outcome_labels, \(x) mu_prior), outcome_labels)
  } else if (is.list(mu_prior)) {
    stats::setNames(
      purrr::map(outcome_labels, \(nm) mu_prior[[nm]] %||% normal(0, 1)),
      outcome_labels
    )
  } else {
    cli::cli_abort("{.arg mu_prior} must be a prior object or named list.",
                   call = rlang::caller_env())
  }

  tau_priors <- if (is.null(tau_prior)) {
    stats::setNames(purrr::map(outcome_labels, \(x) half_cauchy(0, 0.5)), outcome_labels)
  } else if (inherits(tau_prior, "bayesma_prior")) {
    stats::setNames(purrr::map(outcome_labels, \(x) tau_prior), outcome_labels)
  } else if (is.list(tau_prior)) {
    stats::setNames(
      purrr::map(outcome_labels, \(nm) tau_prior[[nm]] %||% half_cauchy(0, 0.5)),
      outcome_labels
    )
  } else {
    cli::cli_abort("{.arg tau_prior} must be a prior object or named list.",
                   call = rlang::caller_env())
  }

  if (is.null(rho_between_prior)) rho_between_prior <- uniform(-1, 1)
  if (!inherits(rho_between_prior, "bayesma_prior")) {
    cli::cli_abort("{.arg rho_between_prior} must be a prior object.",
                   call = rlang::caller_env())
  }

  list(mu = mu_priors, tau = tau_priors, rho_between = rho_between_prior)
}

generate_stan_code_one_stage_mv <- function(outcome_labels, priors) {
  data_block <- "data {
  int<lower=1> S;
  vector[S] y1;
  vector[S] y2;
  vector<lower=0>[S] se1;
  vector<lower=0>[S] se2;
  real<lower=-1, upper=1> rho_within;
}"

  rho_b_bnds <- emit_prior_bounds(priors$rho_between,
                                  default_lower = -1, default_upper = 1)

  par_block <- glue::glue("parameters {{
  real mu1;
  real mu2;
  real<lower=0> tau1;
  real<lower=0> tau2;
  real{rho_b_bnds} rho_between;
}}")

  tp_block <- "transformed parameters {
  cov_matrix[2] D;
  {
    matrix[2, 2] L;
    L[1, 1] = tau1;
    L[1, 2] = 0;
    L[2, 1] = rho_between * tau2;
    L[2, 2] = tau2 * sqrt(1 - square(rho_between));
    D = L * L';
  }
}"

  prior_lines <- c(
    emit_prior_target(priors$mu[[1]],        "mu1"),
    emit_prior_target(priors$mu[[2]],        "mu2"),
    emit_prior_target(priors$tau[[1]],       "tau1"),
    emit_prior_target(priors$tau[[2]],       "tau2"),
    emit_prior_target(priors$rho_between,    "rho_between")
  )
  prior_lines <- prior_lines[nzchar(prior_lines)]
  prior_str   <- paste0("  ", prior_lines, collapse = "\n")

  model_block <- glue::glue("model {{
{prior_str}

  for (s in 1:S) {{
    vector[2] y_s = [y1[s], y2[s]]';
    vector[2] mu_vec = [mu1, mu2]';
    matrix[2, 2] S_i;
    S_i[1, 1] = square(se1[s]);
    S_i[2, 2] = square(se2[s]);
    S_i[1, 2] = rho_within * se1[s] * se2[s];
    S_i[2, 1] = S_i[1, 2];
    matrix[2, 2] Sigma = S_i + D;
    target += multi_normal_lpdf(y_s | mu_vec, Sigma);
  }}
}}")

  gq_block <- "generated quantities {
  vector[2] mu_new;
  vector[S] log_lik;
  mu_new = multi_normal_rng([mu1, mu2]', D);
  for (s in 1:S) {
    vector[2] y_s = [y1[s], y2[s]]';
    vector[2] mu_vec = [mu1, mu2]';
    matrix[2, 2] S_i;
    S_i[1, 1] = square(se1[s]);
    S_i[2, 2] = square(se2[s]);
    S_i[1, 2] = rho_within * se1[s] * se2[s];
    S_i[2, 1] = S_i[1, 2];
    matrix[2, 2] Sigma = S_i + D;
    log_lik[s] = multi_normal_lpdf(y_s | mu_vec, Sigma);
  }
}"

  paste(data_block, par_block, tp_block, model_block, gq_block, sep = "\n\n")
}

generate_stan_code_two_stage_mv <- function(outcome_labels, priors) {
  data_block <- "data {
  int<lower=1> S;
  vector[S] y1;
  vector[S] y2;
  vector<lower=0>[S] se1;
  vector<lower=0>[S] se2;
  real<lower=-1, upper=1> rho_within;
}"

  rho_b_bnds <- emit_prior_bounds(priors$rho_between,
                                  default_lower = -1, default_upper = 1)

  par_block <- glue::glue("parameters {{
  real mu1;
  real mu2;
  real<lower=0> tau1;
  real<lower=0> tau2;
  real{rho_b_bnds} rho_between;
  matrix[2, S] z;
}}")

  tp_block <- "transformed parameters {
  matrix[2, S] theta;
  {
    matrix[2, 2] L;
    L[1, 1] = tau1;
    L[1, 2] = 0;
    L[2, 1] = rho_between * tau2;
    L[2, 2] = tau2 * sqrt(1 - square(rho_between));
    for (s in 1:S) {
      theta[1, s] = mu1 + L[1, 1] * z[1, s];
      theta[2, s] = mu2 + L[2, 1] * z[1, s] + L[2, 2] * z[2, s];
    }
  }
}"

  prior_lines <- c(
    emit_prior_target(priors$mu[[1]],        "mu1"),
    emit_prior_target(priors$mu[[2]],        "mu2"),
    emit_prior_target(priors$tau[[1]],       "tau1"),
    emit_prior_target(priors$tau[[2]],       "tau2"),
    emit_prior_target(priors$rho_between,    "rho_between")
  )
  prior_lines <- prior_lines[nzchar(prior_lines)]
  prior_str   <- paste0("  ", prior_lines, collapse = "\n")

  model_block <- glue::glue("model {{
{prior_str}
  to_vector(z) ~ std_normal();

  for (s in 1:S) {{
    vector[2] y_s = [y1[s], y2[s]]';
    vector[2] theta_s = [theta[1, s], theta[2, s]]';
    matrix[2, 2] S_i;
    S_i[1, 1] = square(se1[s]);
    S_i[2, 2] = square(se2[s]);
    S_i[1, 2] = rho_within * se1[s] * se2[s];
    S_i[2, 1] = S_i[1, 2];
    target += multi_normal_lpdf(y_s | theta_s, S_i);
  }}
}}")

  gq_block <- "generated quantities {
  cov_matrix[2] D;
  vector[2] mu_new;
  vector[S] log_lik;
  D[1, 1] = square(tau1);
  D[2, 2] = square(tau2);
  D[1, 2] = rho_between * tau1 * tau2;
  D[2, 1] = D[1, 2];
  mu_new = multi_normal_rng([mu1, mu2]', D);
  for (s in 1:S) {
    vector[2] y_s = [y1[s], y2[s]]';
    vector[2] theta_s = [theta[1, s], theta[2, s]]';
    matrix[2, 2] S_i;
    S_i[1, 1] = square(se1[s]);
    S_i[2, 2] = square(se2[s]);
    S_i[1, 2] = rho_within * se1[s] * se2[s];
    S_i[2, 1] = S_i[1, 2];
    log_lik[s] = multi_normal_lpdf(y_s | theta_s, S_i);
  }
}"

  paste(data_block, par_block, tp_block, model_block, gq_block, sep = "\n\n")
}
