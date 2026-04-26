

# Generate Stan code — ONE-STAGE

#' @noRd
generate_stan_code_one_stage <- function(likelihood, model_type, re_dist,
                                         priors,
                                         robust_config = list(enabled = FALSE),
                                         multi_arm_config = list(enabled = FALSE)) {

  is_re <- model_type == "random_effect"
  p <- priors
  use_robust <- isTRUE(robust_config$enabled)
  use_multi_arm <- isTRUE(multi_arm_config$enabled)

  # ==== FUNCTIONS ====
  fn_block <- emit_functions_block(re_dist, is_re)

  # ==== DATA ==== (one-stage specific: arm-level)
  data_block <- "data {\n  int<lower=1> N;\n  int<lower=1> S;"
  data_block <- switch(likelihood,
                       binomial = paste0(data_block,
                                         "\n  array[N] int<lower=0> events;\n  array[N] int<lower=1> n;"),
                       gaussian = paste0(data_block,
                                         "\n  vector[N] y;\n  vector<lower=0>[N] se;"),
                       poisson  = paste0(data_block,
                                         "\n  array[N] int<lower=0> events;\n  vector<lower=0>[N] exposure;")
  )
  data_block <- paste0(data_block, "
  array[N] int<lower=0, upper=1> treat;
  array[N] int<lower=1> study;")

  if (use_multi_arm) {
    data_block <- paste0(data_block, "
  // Multi-arm study structure
  int<lower=1> n_ma_studies;
  array[S] int<lower=1, upper=n_ma_studies> comp_to_ma;")
  }

  if (re_dist == "mixture" && is_re)
    data_block <- paste0(data_block, "
  int<lower=2> K;
  vector<lower=0>[K] prior_dirichlet_alpha;")
  data_block <- paste0(data_block, "\n}")

  # ==== PARAMETERS ====
  gamma_bnds <- emit_prior_bounds(p$gamma)
  mu_bnds    <- emit_prior_bounds(p$mu)
  par_block  <- glue::glue(
    "parameters {{\n  vector{gamma_bnds}[S] gamma;\n  real{mu_bnds} mu;")

  if (is_re) {
    par_block <- paste0(as.character(par_block),
                        emit_re_parameters(re_dist, p, include_rho = use_multi_arm))
  }

  if (use_robust) {
    rob_par <- emit_robust_parameters(robust_config)
    par_block <- paste0(as.character(par_block), "\n", rob_par$par)
  }
  par_block <- paste0(as.character(par_block), "\n}")

  # ==== TRANSFORMED PARAMETERS ==== (one-stage specific: epsilon)
  tp_block <- ""
  if (is_re && re_dist %in% c("normal", "t", "skew_normal")) {
    if (use_multi_arm) {
      # Multi-arm: apply within-study correlation using equicorrelation structure
      # epsilon_i = tau * (sqrt(|rho|) * sign(rho) * u_shared + sqrt(1-|rho|) * z_i)
      tp_block <- "transformed parameters {
  vector[S] epsilon;
  {
    // Apply within-study correlation for multi-arm studies
    for (m in 1:n_ma_studies) {
      // Find indices belonging to this MA study
      int n_arms_m = 0;
      for (s in 1:S) {
        if (comp_to_ma[s] == m) n_arms_m += 1;
      }

      if (n_arms_m == 1) {
        // Single comparison in this MA study: standard RE
        for (s in 1:S) {
          if (comp_to_ma[s] == m) {
            epsilon[s] = tau * z[s];
          }
        }
      } else {
        // Multiple comparisons: apply equicorrelation
        real sqrt_abs_rho = sqrt(fabs(rho));
        real sqrt_one_minus_abs_rho = sqrt(1.0 - fabs(rho));
        real sign_rho = rho >= 0 ? 1.0 : -1.0;

        // Use the first comparison's z as the shared component
        real u_shared = 0.0;
        int first_found = 0;
        for (s in 1:S) {
          if (comp_to_ma[s] == m && first_found == 0) {
            u_shared = z[s];
            first_found = 1;
          }
        }

        // Apply correlation structure
        for (s in 1:S) {
          if (comp_to_ma[s] == m) {
            epsilon[s] = tau * (sqrt_abs_rho * sign_rho * u_shared +
                                sqrt_one_minus_abs_rho * z[s]);
          }
        }
      }
    }
  }
}"
    } else {
      tp_block <- "transformed parameters {\n  vector[S] epsilon = tau * z;\n}"
    }
  }

  # ==== MODEL ====
  gamma_tgt <- emit_prior_target(p$gamma, "gamma")
  mu_tgt    <- emit_prior_target(p$mu, "mu")
  priors_code <- glue::glue(
    "  // Priors\n  {gamma_tgt}\n  {mu_tgt}")

  if (use_robust) {
    rob_prior_code <- emit_robust_priors(robust_config)
    priors_code <- paste0(as.character(priors_code), "\n  ", rob_prior_code)
  }

  if (is_re) {
    re_code <- emit_re_priors(re_dist, p, include_rho = use_multi_arm)
    priors_code <- paste0(as.character(priors_code), re_code)
  }

  # Likelihood (one-stage specific: arm-level)
  if (is_re && re_dist == "mixture") {
    # Mixture one-stage — robust wrapping not supported for mixture one-stage
    lik_code <- switch(likelihood,
                       binomial = "
  for (i in 1:N) {
    if (treat[i] == 0) {
      target += binomial_logit_lpmf(events[i:i] | n[i:i], rep_vector(gamma[study[i]], 1));
    } else {
      vector[K] lps;
      for (k in 1:K) {
        real eta_k = gamma[study[i]] + (mu + delta_k[k]) * treat[i];
        lps[k] = log(w[k]) + binomial_logit_lpmf(events[i:i] | n[i:i], rep_vector(eta_k, 1));
      }
      target += log_sum_exp(lps);
    }
  }",
                       gaussian = "
  for (i in 1:N) {
    if (treat[i] == 0) {
      target += normal_lpdf(y[i] | gamma[study[i]], se[i]);
    } else {
      vector[K] lps;
      for (k in 1:K) {
        real eta_k = gamma[study[i]] + (mu + delta_k[k]) * treat[i];
        lps[k] = log(w[k]) + normal_lpdf(y[i] | eta_k, se[i]);
      }
      target += log_sum_exp(lps);
    }
  }",
                       poisson = "
  for (i in 1:N) {
    if (treat[i] == 0) {
      target += poisson_log_lpmf(events[i:i] | rep_vector(log(exposure[i]) + gamma[study[i]], 1));
    } else {
      vector[K] lps;
      for (k in 1:K) {
        real eta_k = gamma[study[i]] + (mu + delta_k[k]) * treat[i];
        lps[k] = log(w[k]) + poisson_log_lpmf(events[i:i] | rep_vector(log(exposure[i]) + eta_k, 1));
      }
      target += log_sum_exp(lps);
    }
  }"
    )
  } else {
    eta_term <- if (is_re) {
      "gamma[study[i]] + mu * treat[i] + epsilon[study[i]] * treat[i]"
    } else {
      "gamma[study[i]] + mu * treat[i]"
    }

    if (use_robust && likelihood == "gaussian") {
      lik_code <- glue::glue("
  for (i in 1:N) {{
    real eta_i = {eta_term};
    real main_ll = normal_lpdf(y[i] | eta_i, se[i]);
    target += {emit_robust_likelihood('main_ll', 'y[i]', 'se[i]', robust_config)};
  }}")
    } else if (use_robust && likelihood == "binomial") {
      lik_code <- glue::glue("
  for (i in 1:N) {{
    real eta_i = {eta_term};
    real main_ll = binomial_logit_lpmf(events[i:i] | n[i:i], rep_vector(eta_i, 1));
    real p_hat = inv_logit(eta_i);
    real se_approx = sqrt(1.0 / (n[i] * p_hat * (1 - p_hat) + 0.001));
    real y_approx = eta_i;
    real outlier_ll = {if (is.finite(robust_config$df))
      glue::glue('student_t_lpdf(y_approx | {robust_config$df}, {robust_config$prior$mean}, se_approx * 3)')
    else
      glue::glue('normal_lpdf(y_approx | {robust_config$prior$mean}, se_approx * 3)')};
    target += log_sum_exp(log(pi_main) + main_ll, log1m(pi_main) + outlier_ll);
  }}")
    } else if (use_robust && likelihood == "poisson") {
      lik_code <- glue::glue("
  for (i in 1:N) {{
    real eta_i = {eta_term};
    real main_ll = poisson_log_lpmf(events[i:i] | rep_vector(log(exposure[i]) + eta_i, 1));
    real rate_hat = exp(eta_i) * exposure[i];
    real se_approx = sqrt(1.0 / (rate_hat + 0.001));
    real y_approx = eta_i;
    real outlier_ll = {if (is.finite(robust_config$df))
      glue::glue('student_t_lpdf(y_approx | {robust_config$df}, {robust_config$prior$mean}, se_approx * 3)')
    else
      glue::glue('normal_lpdf(y_approx | {robust_config$prior$mean}, se_approx * 3)')};
    target += log_sum_exp(log(pi_main) + main_ll, log1m(pi_main) + outlier_ll);
  }}")
    } else {
      lik_stmt <- switch(likelihood,
                         binomial = "target += binomial_logit_lpmf(events | n, eta);",
                         gaussian = "target += normal_lpdf(y | eta, se);",
                         poisson  = "target += poisson_log_lpmf(events | log(exposure) + eta);"
      )
      lik_code <- paste0(
        "  {\n",
        "    vector[N] eta;\n",
        "    for (i in 1:N) {\n",
        "      eta[i] = ", eta_term, ";\n",
        "    }\n",
        "    ", lik_stmt, "\n",
        "  }")
    }
  }

  model_block <- paste0("model {\n", as.character(priors_code),
                        "\n  // Likelihood\n",
                        as.character(lik_code), "\n}")

  # ==== GENERATED QUANTITIES ====
  gq <- switch(likelihood,
               binomial = "  real log_or = mu;\n  real or_pooled = exp(mu);",
               gaussian = "  real pooled_diff = mu;",
               poisson  = "  real log_rr = mu;\n  real rr_pooled = exp(mu);"
  )
  if (is_re) {
    gq <- paste0(gq, emit_re_extra_gq(re_dist, "one_stage"))

    if (use_multi_arm) {
      gq <- paste0(gq, "\n  real rho_est = rho;")
    }
  }

  if (use_robust) {
    gq <- paste0(gq, "\n  real pi_main_est = pi_main;")
  }

  # ---- log_lik for model comparison (LOO-CV, WAIC) ----
  if (is_re && re_dist %in% c("normal", "t", "skew_normal")) {
    eta_ll <- "gamma[study[i]] + mu * treat[i] + epsilon[study[i]] * treat[i]"
  } else {
    eta_ll <- "gamma[study[i]] + mu * treat[i]"
  }

  log_lik_gq <- switch(likelihood,
                       gaussian = glue::glue("
  vector[N] log_lik;
  for (i in 1:N) {{
    real eta_i = {eta_ll};
    log_lik[i] = normal_lpdf(y[i] | eta_i, se[i]);
  }}"),
                       binomial = glue::glue("
  vector[N] log_lik;
  for (i in 1:N) {{
    real eta_i = {eta_ll};
    log_lik[i] = binomial_lpmf(events[i] | n[i], inv_logit(eta_i));
  }}"),
                       poisson = glue::glue("
  vector[N] log_lik;
  for (i in 1:N) {{
    real eta_i = {eta_ll};
    log_lik[i] = poisson_log_lpmf(events[i] | log(exposure[i]) + eta_i);
  }}")
  )

  gq <- paste0(gq, as.character(log_lik_gq))

  # ---- y_rep for posterior predictive checks ----
  if (is_re && re_dist %in% c("normal", "t", "skew_normal")) {
    eta_re <- "gamma[study[i]] + mu * treat[i] + epsilon[study[i]] * treat[i]"
  } else {
    eta_re <- "gamma[study[i]] + mu * treat[i]"
  }

  y_rep_gq <- switch(likelihood,
                     gaussian = glue::glue("
  vector[N] y_rep;
  for (i in 1:N) {{
    real eta_i = {eta_re};
    y_rep[i] = normal_rng(eta_i, se[i]);
  }}"),
                     binomial = glue::glue("
  vector[N] y_rep;
  for (i in 1:N) {{
    real eta_i = {eta_re};
    y_rep[i] = binomial_rng(n[i], inv_logit(eta_i)) * 1.0;
  }}"),
                     poisson = glue::glue("
  vector[N] y_rep;
  for (i in 1:N) {{
    real eta_i = {eta_re};
    y_rep[i] = poisson_log_rng(log(exposure[i]) + eta_i) * 1.0;
  }}")
  )

  gq <- paste0(gq, as.character(y_rep_gq))

  # Add prediction interval for RE models
  if (is_re) {
    gq <- paste0(gq, as.character(emit_mu_new_gq(re_dist, p)))
  }

  gq_block <- glue::glue("generated quantities {{\n{gq}\n}}")

  # ==== ASSEMBLE ====
  assemble_stan_program(fn_block, data_block, par_block, tp_block,
                        model_block, gq_block)
}


# Generate Stan code — TWO-STAGE

#' @noRd
generate_stan_code_two_stage <- function(model_type, re_dist,
                                         use_t_likelihood, priors,
                                         robust_config = list(enabled = FALSE)) {

  is_re <- model_type == "random_effect"
  p <- priors
  use_robust <- isTRUE(robust_config$enabled)

  # ==== FUNCTIONS ====
  fn_block <- emit_functions_block(re_dist, is_re)

  # ==== DATA ==== (two-stage specific: study-level)
  data_block <- "data {
  int<lower=1> S;
  vector[S] y;
  vector<lower=0>[S] se;"
  if (use_t_likelihood)
    data_block <- paste0(data_block, "\n  vector<lower=1>[S] df;")
  if (re_dist == "mixture" && is_re)
    data_block <- paste0(data_block, "
  int<lower=2> K;
  vector<lower=0>[K] prior_dirichlet_alpha;")
  data_block <- paste0(data_block, "\n}")

  # ==== PARAMETERS ====
  mu_bnds <- emit_prior_bounds(p$mu)
  par_block <- glue::glue("parameters {{\n  real{mu_bnds} mu;")

  if (is_re) {
    par_block <- paste0(as.character(par_block),
                        emit_re_parameters(re_dist, p))
  }

  if (use_robust) {
    rob_par <- emit_robust_parameters(robust_config)
    par_block <- paste0(as.character(par_block), "\n", rob_par$par)
  }
  par_block <- paste0(as.character(par_block), "\n}")

  # ==== TRANSFORMED PARAMETERS ==== (two-stage specific: theta)
  tp_block <- ""
  if (is_re && re_dist %in% c("normal", "t", "skew_normal"))
    tp_block <- "transformed parameters {\n  vector[S] theta = mu + tau * z;\n}"

  # ==== MODEL ====
  mu_tgt <- emit_prior_target(p$mu, "mu")
  priors_code <- paste0("  // Priors\n  ", mu_tgt)

  if (use_robust) {
    rob_prior_code <- emit_robust_priors(robust_config)
    priors_code <- paste0(priors_code, "\n  ", rob_prior_code)
  }

  if (is_re) {
    re_code <- emit_re_priors(re_dist, p)
    priors_code <- paste0(priors_code, re_code)
  }

  # Likelihood (two-stage specific: study-level)
  if (is_re && re_dist == "mixture") {
    if (use_robust) {
      lik_code <- glue::glue("
  for (i in 1:S) {{
    vector[K] lps;
    for (k in 1:K)
      lps[k] = log(w[k]) + normal_lpdf(y[i] | mu_k[k], sqrt(square(tau_k[k]) + square(se[i])));
    real main_ll = log_sum_exp(lps);
    target += {emit_robust_likelihood('main_ll', 'y[i]', 'se[i]', robust_config)};
  }}")
    } else {
      lik_code <- "
  for (i in 1:S) {
    vector[K] lps;
    for (k in 1:K)
      lps[k] = log(w[k]) + normal_lpdf(y[i] | mu_k[k], sqrt(square(tau_k[k]) + square(se[i])));
    target += log_sum_exp(lps);
  }"
    }
  } else {
    loc <- if (is_re) "theta" else "rep_vector(mu, S)"
    if (use_robust) {
      if (use_t_likelihood) {
        lik_code <- glue::glue("
  for (i in 1:S) {{
    real main_ll = student_t_lpdf(y[i] | df[i], {loc}[i], se[i]);
    target += {emit_robust_likelihood('main_ll', 'y[i]', 'se[i]', robust_config)};
  }}")
      } else {
        lik_code <- glue::glue("
  for (i in 1:S) {{
    real main_ll = normal_lpdf(y[i] | {loc}[i], se[i]);
    target += {emit_robust_likelihood('main_ll', 'y[i]', 'se[i]', robust_config)};
  }}")
      }
    } else {
      if (use_t_likelihood) {
        lik_code <- glue::glue("
  for (i in 1:S) {{
    target += student_t_lpdf(y[i] | df[i], {loc}[i], se[i]);
  }}")
      } else {
        lik_code <- glue::glue(
          "\n  target += normal_lpdf(y | {loc}, se);")
      }
    }
  }

  model_block <- paste0("model {\n", priors_code,
                        "\n  // Likelihood\n",
                        as.character(lik_code), "\n}")

  # ==== GENERATED QUANTITIES ====
  gq <- "  real pooled = mu;"
  if (is_re) {
    gq <- paste0(gq, emit_re_extra_gq(re_dist, "two_stage"))
    gq <- paste0(gq, as.character(emit_mu_new_gq(re_dist, p)))
  }

  # Add robust GQ (two-stage specific: prob_outlier)
  if (use_robust) {
    gq <- paste0(gq, "\n  real pi_main_est = pi_main;")
    if (is_re && re_dist %in% c("normal", "t", "skew_normal")) {
      outlier_ll_gq <- if (is.finite(robust_config$df)) {
        glue::glue("student_t_lpdf(y[i] | {robust_config$df}, {robust_config$prior$mean}, se[i] * 3)")
      } else {
        glue::glue("normal_lpdf(y[i] | {robust_config$prior$mean}, se[i] * 3)")
      }
      gq <- paste0(gq, "
  vector[S] prob_outlier;
  for (i in 1:S) {
    real ll_main = normal_lpdf(y[i] | theta[i], se[i]);
    real ll_out = ", as.character(outlier_ll_gq), ";
    real lp_main = log(pi_main) + ll_main;
    real lp_out = log1m(pi_main) + ll_out;
    prob_outlier[i] = exp(lp_out - log_sum_exp(lp_main, lp_out));
  }")
    }
  }

  # ---- log_lik for model comparison (LOO-CV, WAIC) ----
  if (is_re && re_dist %in% c("normal", "t", "skew_normal")) {
    if (use_t_likelihood) {
      gq <- paste0(gq, "
  vector[S] log_lik;
  for (i in 1:S)
    log_lik[i] = student_t_lpdf(y[i] | df[i], theta[i], se[i]);")
    } else {
      gq <- paste0(gq, "
  vector[S] log_lik;
  for (i in 1:S)
    log_lik[i] = normal_lpdf(y[i] | theta[i], se[i]);")
    }
  } else if (is_re && re_dist == "mixture") {
    gq <- paste0(gq, "
  vector[S] log_lik;
  for (i in 1:S) {
    vector[K] lps;
    for (k in 1:K)
      lps[k] = log(w[k]) + normal_lpdf(y[i] | mu_k[k], sqrt(square(tau_k[k]) + square(se[i])));
    log_lik[i] = log_sum_exp(lps);
  }")
  } else {
    if (use_t_likelihood) {
      gq <- paste0(gq, "
  vector[S] log_lik;
  for (i in 1:S)
    log_lik[i] = student_t_lpdf(y[i] | df[i], mu, se[i]);")
    } else {
      gq <- paste0(gq, "
  vector[S] log_lik;
  for (i in 1:S)
    log_lik[i] = normal_lpdf(y[i] | mu, se[i]);")
    }
  }

  # ---- y_rep for posterior predictive checks ----
  if (is_re && re_dist %in% c("normal", "t", "skew_normal")) {
    gq <- paste0(gq, "
  vector[S] y_rep;
  for (i in 1:S)
    y_rep[i] = normal_rng(theta[i], se[i]);")
  } else if (is_re && re_dist == "mixture") {
    gq <- paste0(gq, "
  vector[S] y_rep;
  for (i in 1:S) {
    int comp = cluster[i];
    y_rep[i] = normal_rng(mu_k[comp], sqrt(square(tau_k[comp]) + square(se[i])));
  }")
  } else {
    gq <- paste0(gq, "
  vector[S] y_rep;
  for (i in 1:S)
    y_rep[i] = normal_rng(mu, se[i]);")
  }

  gq_block <- glue::glue("generated quantities {{\n{gq}\n}}")

  # ==== ASSEMBLE ====
  assemble_stan_program(fn_block, data_block, par_block, tp_block,
                        model_block, gq_block)
}


# Generate Stan code — TWO-STAGE NULL (mu fixed at 0)- for RoBMA

#' @noRd
generate_stan_code_two_stage_null <- function(model_type, re_dist, priors) {

  is_re <- model_type == "random_effect"
  p <- priors

  # ==== DATA ====
  data_block <- "data {
  int<lower=1> S;
  vector[S] y;
  vector<lower=0>[S] se;
}"

  # ==== PARAMETERS ====
  if (is_re) {
    re_params <- emit_re_parameters(re_dist, p)
    # Re-indent: emitter returns un-indented lines (glue trim strips leading spaces)
    re_params_indented <- gsub("(^|\n)", "\\1  ", re_params)
    par_block <- paste0("parameters {\n", re_params_indented, "\n}")
  } else {
    # Common effect with mu = 0: no free parameters
    # Use a dummy parameter so Stan can compile (Stan requires >= 1 param)
    par_block <- "parameters {
  real<lower=0, upper=0.001> dummy__;
}"
  }

  # ==== TRANSFORMED PARAMETERS ====
  if (is_re && re_dist %in% c("normal", "t")) {
    tp_block <- "transformed parameters {
  vector[S] theta = tau * z;  // mu = 0, so theta = 0 + tau * z
}"
  } else {
    tp_block <- ""
  }

  # ==== MODEL ====
  if (is_re) {
    re_prior_code <- emit_re_priors(re_dist, p)

    model_block <- paste0(
      "model {\n  // Priors\n",
      as.character(re_prior_code),
      "\n  // Likelihood (mu = 0)\n",
      "  target += normal_lpdf(y | theta, se);\n}"
    )
  } else {
    # Common effect, mu = 0: pure likelihood, no parameters to put priors on
    model_block <- "model {
  // Likelihood (mu = 0, common effect)
  target += normal_lpdf(y | rep_vector(0.0, S), se);
}"
  }

  # ==== GENERATED QUANTITIES ====
  # Null-specific: mu fixed at 0, mu_new uses 0 as centre
  gq_parts <- "  real mu = 0.0;"

  if (is_re) {
    gq_parts <- paste0(gq_parts, "
  vector[S] shrinkage = theta;
  real mu_new = normal_rng(0.0, tau);
  vector[S] y_rep;
  for (i in 1:S)
    y_rep[i] = normal_rng(theta[i], se[i]);")
  } else {
    gq_parts <- paste0(gq_parts, "
  vector[S] y_rep;
  for (i in 1:S)
    y_rep[i] = normal_rng(0.0, se[i]);")
  }

  gq_block <- paste0("generated quantities {\n", gq_parts, "\n}")

  # ==== ASSEMBLE ====
  assemble_stan_program("", data_block, par_block, tp_block, model_block, gq_block)
}


# Generate Stan code — BIAS-CORRECTED (Jung)

#' @noRd
generate_stan_code_bias_corrected <- function(priors) {

  p <- priors

  # B bounds from prior (uniform) or default
  b_bounds <- emit_prior_bounds(p$b, default_lower = 0)

  # tau: inverse variance parameterisation
  tau_tgt <- emit_prior_target(p$tau, "inv_var")
  mu_tgt  <- emit_prior_target(p$mu, "mu")
  b_tgt   <- emit_prior_target(p$b, "B")
  pb_tgt  <- emit_prior_target(p$p_bias, "p_bias")
  wb_tgt  <- emit_prior_target(p$w_bias, "w_bias_raw")

  data_block <- "data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se_y;
  int<lower=0, upper=1> use_known_bias;
  array[N] int<lower=0, upper=1> known_bias;
}"

  par_block <- paste0("parameters {
  real mu;
  real", b_bounds, " B;
  real<lower=0, upper=1> p_bias;
  real<lower=0> inv_var;
  vector[N] theta_raw_unbiased;
  vector[N] theta_raw_biased;
  vector<lower=0, upper=1>[N] w_bias_raw;
}")

  tp_block <- "transformed parameters {
  real tau = 1.0 / sqrt(inv_var);
  real mu_biased = mu + B;
  vector[N] theta_unbiased;
  vector[N] theta_biased;
  theta_unbiased = mu + tau * theta_raw_unbiased;
  for (n in 1:N) {
    real tau_biased_n = 1.0 / sqrt(inv_var * w_bias_raw[n]);
    theta_biased[n] = mu_biased + tau_biased_n * theta_raw_biased[n];
  }
}"

  model_block <- paste0("model {
  // Priors
  ", mu_tgt, "
  ", b_tgt, "
  ", pb_tgt, "
  ", tau_tgt, "
  ", wb_tgt, "
  target += std_normal_lpdf(theta_raw_unbiased);
  target += std_normal_lpdf(theta_raw_biased);
  // Likelihood
  for (n in 1:N) {
    if (use_known_bias) {
      if (known_bias[n] == 1) {
        target += normal_lpdf(y[n] | theta_biased[n], se_y[n]);
      } else {
        target += normal_lpdf(y[n] | theta_unbiased[n], se_y[n]);
      }
    } else {
      real lp_unbiased = log1m(p_bias)
        + normal_lpdf(y[n] | theta_unbiased[n], se_y[n]);
      real lp_biased = log(p_bias)
        + normal_lpdf(y[n] | theta_biased[n], se_y[n]);
      target += log_sum_exp(lp_unbiased, lp_biased);
    }
  }
}")

  gq_block <- "generated quantities {
  vector[N] y_rep;
  vector[N] prob_biased;
  real or_unbiased = exp(mu);
  real or_biased = exp(mu_biased);
  for (n in 1:N) {
    real lp_unbiased = log1m(p_bias)
      + normal_lpdf(y[n] | theta_unbiased[n], se_y[n]);
    real lp_biased = log(p_bias)
      + normal_lpdf(y[n] | theta_biased[n], se_y[n]);
    if (use_known_bias) {
      prob_biased[n] = known_bias[n];
    } else {
      prob_biased[n] = exp(lp_biased - log_sum_exp(lp_unbiased, lp_biased));
    }
    if (bernoulli_rng(prob_biased[n] > 0.5 ? prob_biased[n] : 0.0 + known_bias[n])) {
      y_rep[n] = normal_rng(theta_biased[n], se_y[n]);
    } else {
      y_rep[n] = normal_rng(theta_unbiased[n], se_y[n]);
    }
  }
}"

  paste(data_block, par_block, tp_block, model_block, gq_block,
        sep = "\n\n")
}

#' Generate Stan code for Bayesian Egger test with binomial likelihood
#' @noRd
generate_bayesian_egger_stan_binomial <- function(heterogeneity, priors) {

  # Prior specification helpers
  prior_to_stan <- function(p, param_name) {
    fam <- p$family
    switch(fam,
           normal = glue::glue("{param_name} ~ normal({p$mean}, {p$sd});"),
           half_normal = glue::glue("{param_name} ~ normal({p$mean}, {p$sd});"),
           uniform = glue::glue("{param_name} ~ uniform({p$lower}, {p$upper});"),
           half_cauchy = glue::glue("{param_name} ~ cauchy({p$location}, {p$scale});"),
           exponential = glue::glue("{param_name} ~ exponential({p$rate});"),
           glue::glue("{param_name} ~ normal(0, 100);")
    )
  }

  if (heterogeneity == "multiplicative") {
    stan_code <- glue::glue('
data {{
  int<lower=1> N;           // number of studies
  array[N] int<lower=0> n0; // sample size control
  array[N] int<lower=0> n1; // sample size treatment
  array[N] int<lower=0> r0; // events control
  array[N] int<lower=0> r1; // events treatment
  vector[N] y;              // observed log OR
  real<lower=0> upp_kappa;  // upper bound for kappa
  real<lower=0> upp_tau;    // upper bound for tau
}}

transformed data {{
  // Small constant to prevent division by zero
  real eps = 1e-6;
}}

parameters {{
  real alpha;                          // intercept
  real beta;                           // slope (small-study effect)
  real<lower=0, upper=upp_kappa> kappa; // multiplicative heterogeneity
  real d;                              // overall log OR
  real<lower=1e-4, upper=upp_tau> tau; // between-study SD (lower bounded for stability)
  vector[N] delta;                     // study-specific log ORs
  vector<lower=eps, upper=1-eps>[N] p0; // true event rate control (bounded away from 0/1)
}}

transformed parameters {{
  vector<lower=eps, upper=1-eps>[N] p1; // true event rate treatment
  vector<lower=0>[N] sigma;             // true within-study SE

  for (i in 1:N) {{
    // Compute p1 with bounds checking
    real logit_p1 = logit(p0[i]) + delta[i];
    // Constrain logit_p1 to avoid extreme values
    logit_p1 = fmax(fmin(logit_p1, 10.0), -10.0);
    p1[i] = inv_logit(logit_p1);

    // Compute sigma with numerical stability
    real var0 = 1.0 / (n0[i] * p0[i] * (1 - p0[i]));
    real var1 = 1.0 / (n1[i] * p1[i] * (1 - p1[i]));
    sigma[i] = sqrt(var0 + var1);
  }}
}}

model {{
  // Priors
  {prior_to_stan(priors$alpha, "alpha")}
  {prior_to_stan(priors$beta, "beta")}
  {prior_to_stan(priors$kappa, "kappa")}
  {prior_to_stan(priors$d, "d")}
  {prior_to_stan(priors$tau, "tau")}

  // Baseline rates - weakly informative, concentrated away from extremes
  p0 ~ beta(2, 2);

  // Hierarchical model for log ORs
  delta ~ normal(d, tau);

  // Data likelihood (binomial)
  r0 ~ binomial(n0, p0);
  r1 ~ binomial(n1, p1);

  // Egger regression with latent SEs (multiplicative heterogeneity)
  for (i in 1:N) {{
    y[i] ~ normal(alpha + beta * sigma[i], kappa * sigma[i]);
  }}
}}

generated quantities {{
  real d_est = d;
  real tau_est = tau;
  real beta_est = beta;
  vector[N] sigma_est = sigma;
}}
')
  } else {
    # Additive heterogeneity
    stan_code <- glue::glue('
data {{
  int<lower=1> N;           // number of studies
  array[N] int<lower=0> n0; // sample size control
  array[N] int<lower=0> n1; // sample size treatment
  array[N] int<lower=0> r0; // events control
  array[N] int<lower=0> r1; // events treatment
  vector[N] y;              // observed log OR
  real<lower=0> upp_gamma;  // upper bound for gamma
  real<lower=0> upp_tau;    // upper bound for tau
}}

transformed data {{
  // Small constant to prevent division by zero
  real eps = 1e-6;
}}

parameters {{
  real alpha;                          // intercept
  real beta;                           // slope (small-study effect)
  real<lower=0, upper=upp_gamma> gamma; // additive heterogeneity SD
  real d;                              // overall log OR
  real<lower=1e-4, upper=upp_tau> tau; // between-study SD (lower bounded for stability)
  vector[N] delta;                     // study-specific log ORs
  vector<lower=eps, upper=1-eps>[N] p0; // true event rate control (bounded away from 0/1)
}}

transformed parameters {{
  vector<lower=eps, upper=1-eps>[N] p1; // true event rate treatment
  vector<lower=0>[N] sigma;             // true within-study SE

  for (i in 1:N) {{
    // Compute p1 with bounds checking
    real logit_p1 = logit(p0[i]) + delta[i];
    // Constrain logit_p1 to avoid extreme values
    logit_p1 = fmax(fmin(logit_p1, 10.0), -10.0);
    p1[i] = inv_logit(logit_p1);

    // Compute sigma with numerical stability
    real var0 = 1.0 / (n0[i] * p0[i] * (1 - p0[i]));
    real var1 = 1.0 / (n1[i] * p1[i] * (1 - p1[i]));
    sigma[i] = sqrt(var0 + var1);
  }}
}}

model {{
  // Priors
  {prior_to_stan(priors$alpha, "alpha")}
  {prior_to_stan(priors$beta, "beta")}
  {prior_to_stan(priors$gamma, "gamma")}
  {prior_to_stan(priors$d, "d")}
  {prior_to_stan(priors$tau, "tau")}

  // Baseline rates - weakly informative, concentrated away from extremes
  p0 ~ beta(2, 2);

  // Hierarchical model for log ORs
  delta ~ normal(d, tau);

  // Data likelihood (binomial)
  r0 ~ binomial(n0, p0);
  r1 ~ binomial(n1, p1);

  // Egger regression with latent SEs (additive heterogeneity)
  for (i in 1:N) {{
    y[i] ~ normal(alpha + beta * sigma[i], sqrt(square(sigma[i]) + square(gamma)));
  }}
}}

generated quantities {{
  real d_est = d;
  real tau_est = tau;
  real beta_est = beta;
  real gamma_est = gamma;
  vector[N] sigma_est = sigma;
}}
')
  }

  stan_code
}


#' Generate Stan code for generic Bayesian Egger test (continuous/Poisson)
#' @noRd
generate_bayesian_egger_stan_generic <- function(heterogeneity, priors) {

  prior_to_stan <- function(p, param_name) {
    fam <- p$family
    switch(fam,
           normal = glue::glue("{param_name} ~ normal({p$mean}, {p$sd});"),
           half_normal = glue::glue("{param_name} ~ normal({p$mean}, {p$sd});"),
           uniform = glue::glue("{param_name} ~ uniform({p$lower}, {p$upper});"),
           half_cauchy = glue::glue("{param_name} ~ cauchy({p$location}, {p$scale});"),
           exponential = glue::glue("{param_name} ~ exponential({p$rate});"),
           glue::glue("{param_name} ~ normal(0, 100);")
    )
  }

  if (heterogeneity == "multiplicative") {
    stan_code <- glue::glue('
data {{
  int<lower=1> N;          // number of studies
  vector[N] y;             // observed effect sizes
  vector<lower=0>[N] se;   // observed standard errors
  real<lower=0> upp_kappa; // upper bound for kappa
}}

parameters {{
  real alpha;                           // intercept
  real beta;                            // slope (small-study effect)
  real<lower=0, upper=upp_kappa> kappa; // multiplicative heterogeneity
}}

model {{
  // Priors
  {prior_to_stan(priors$alpha, "alpha")}
  {prior_to_stan(priors$beta, "beta")}
  {prior_to_stan(priors$kappa, "kappa")}

  // Egger regression (multiplicative heterogeneity)
  for (i in 1:N) {{
    y[i] ~ normal(alpha + beta * se[i], kappa * se[i]);
  }}
}}

generated quantities {{
  real beta_est = beta;
}}
')
  } else {
    stan_code <- glue::glue('
data {{
  int<lower=1> N;          // number of studies
  vector[N] y;             // observed effect sizes
  vector<lower=0>[N] se;   // observed standard errors
  real<lower=0> upp_gamma; // upper bound for gamma
}}

parameters {{
  real alpha;                           // intercept
  real beta;                            // slope (small-study effect)
  real<lower=0, upper=upp_gamma> gamma; // additive heterogeneity SD
}}

model {{
  // Priors
  {prior_to_stan(priors$alpha, "alpha")}
  {prior_to_stan(priors$beta, "beta")}
  {prior_to_stan(priors$gamma, "gamma")}

  // Egger regression (additive heterogeneity)
  for (i in 1:N) {{
    y[i] ~ normal(alpha + beta * se[i], sqrt(square(se[i]) + square(gamma)));
  }}
}}

generated quantities {{
  real beta_est = beta;
  real gamma_est = gamma;
}}
')
  }

  stan_code
}

# Generate Stan code for one-stage meta-regression

#' @noRd
generate_stan_code_mreg_one_stage <- function(likelihood, model_type, priors, K) {

  is_re <- model_type == "random_effect"

  # ---- DATA BLOCK ----
  data_lines <- c(
    "  int<lower=1> N;           // number of observations (arms)",
    "  int<lower=1> S;           // number of studies",
    "  int<lower=1> K;           // number of moderators",
    "  array[N] int<lower=0, upper=1> treat;  // treatment indicator",
    "  array[N] int<lower=1, upper=S> study;  // study indicator",
    "  matrix[N, K] X;           // design matrix (arm-level)"
  )

  if (likelihood == "binomial") {
    data_lines <- c(data_lines,
                    "  array[N] int<lower=0> events;  // event counts",
                    "  array[N] int<lower=1> n;       // sample sizes")
  } else if (likelihood == "gaussian") {
    data_lines <- c(data_lines,
                    "  vector[N] y;                   // outcomes",
                    "  vector<lower=0>[N] se;         // standard errors")
  } else if (likelihood == "poisson") {
    data_lines <- c(data_lines,
                    "  array[N] int<lower=0> events;  // event counts",
                    "  vector<lower=0>[N] exposure;   // person-time or N")
  }

  data_block <- paste0("data {\n", paste(data_lines, collapse = "\n"), "\n}")

  # ---- PARAMETERS BLOCK ----
  tau_bounds <- if (is_re) emit_prior_bounds(priors$tau, default_lower = 0) else ""

  par_lines <- c(
    "  real mu;                  // intercept (pooled treatment effect)",
    "  vector[K] beta;           // regression coefficients",
    "  vector[S] gamma;          // study baselines"
  )

  if (is_re) {
    par_lines <- c(par_lines,
                   glue::glue("  real{tau_bounds} tau;     // residual heterogeneity"),
                   "  vector[S] z;              // standardized random effects")
  }

  par_block <- paste0("parameters {\n", paste(par_lines, collapse = "\n"), "\n}")

  # ---- TRANSFORMED PARAMETERS BLOCK ----
  if (is_re) {
    tp_block <- "transformed parameters {
  vector[S] epsilon;        // study-level random effects
  epsilon = tau * z;
}"
  } else {
    tp_block <- ""
  }

  # ---- MODEL BLOCK ----
  model_lines <- character()

  # Priors
  mu_tgt <- emit_prior_target(priors$mu, "mu")
  if (nzchar(mu_tgt)) model_lines <- c(model_lines, paste0("  ", mu_tgt))

  gamma_tgt <- emit_prior_target(priors$gamma, "gamma")
  if (nzchar(gamma_tgt)) model_lines <- c(model_lines, paste0("  ", gamma_tgt))

  # Beta priors
  unique_beta_priors <- unique(purrr::map_chr(priors$beta, format.bayesma_prior))
  if (length(unique_beta_priors) == 1) {
    first_prior <- priors$beta[[1]]
    beta_tgt <- emit_prior_target(first_prior, "beta")
    if (nzchar(beta_tgt)) model_lines <- c(model_lines, paste0("  ", beta_tgt))
  } else {
    for (k in seq_len(K)) {
      beta_tgt <- emit_prior_target(priors$beta[[k]], paste0("beta[", k, "]"))
      if (nzchar(beta_tgt)) model_lines <- c(model_lines, paste0("  ", beta_tgt))
    }
  }

  if (is_re) {
    tau_tgt <- emit_prior_target(priors$tau, "tau")
    if (nzchar(tau_tgt)) model_lines <- c(model_lines, paste0("  ", tau_tgt))
    model_lines <- c(model_lines, "  z ~ std_normal();")
  }

  # Likelihood
  if (is_re) {
    eta_expr <- "gamma[study[i]] + (mu + X[i] * beta + epsilon[study[i]]) * treat[i]"
  } else {
    eta_expr <- "gamma[study[i]] + (mu + X[i] * beta) * treat[i]"
  }

  model_lines <- c(model_lines, "  for (i in 1:N) {",
                   glue::glue("    real eta = {eta_expr};"))

  ll_stmt <- switch(likelihood,
                    binomial = "    events[i] ~ binomial_logit(n[i], eta);",
                    gaussian = "    y[i] ~ normal(eta, se[i]);",
                    poisson  = "    events[i] ~ poisson_log(log(exposure[i]) + eta);")

  model_lines <- c(model_lines, ll_stmt, "  }")

  model_block <- paste0("model {\n", paste(model_lines, collapse = "\n"), "\n}")

  # ---- GENERATED QUANTITIES BLOCK ----
  gq_lines <- character()

  if (is_re) {
    gq_lines <- c(gq_lines, "  real mu_new = normal_rng(mu, tau);")
  }

  # y_rep for posterior predictive checks
  gq_lines <- c(gq_lines, "  vector[N] y_rep;", "  for (i in 1:N) {")

  if (is_re) {
    gq_lines <- c(gq_lines,
                  glue::glue("    real eta_i = {eta_expr};"))
  } else {
    gq_lines <- c(gq_lines,
                  glue::glue("    real eta_i = {eta_expr};"))
  }

  y_rep_stmt <- switch(likelihood,
                       binomial = "    y_rep[i] = binomial_rng(n[i], inv_logit(eta_i)) * 1.0;",
                       gaussian = "    y_rep[i] = normal_rng(eta_i, se[i]);",
                       poisson  = "    y_rep[i] = poisson_log_rng(log(exposure[i]) + eta_i) * 1.0;")

  gq_lines <- c(gq_lines, y_rep_stmt, "  }")

  gq_block <- if (length(gq_lines) > 0) {
    paste0("generated quantities {\n", paste(gq_lines, collapse = "\n"), "\n}")
  } else {
    ""
  }

  # ---- ASSEMBLE ----
  parts <- c(data_block, par_block)
  if (nzchar(tp_block)) parts <- c(parts, tp_block)
  parts <- c(parts, model_block)
  if (nzchar(gq_block)) parts <- c(parts, gq_block)
  paste(parts, collapse = "\n\n")
}


# Generate Stan code for two-stage meta-regression

#' @noRd
generate_stan_code_mreg_two_stage <- function(model_type, use_t_likelihood,
                                              priors, K) {

  is_re <- model_type == "random_effect"

  # ---- DATA BLOCK ----
  data_lines <- c(
    "  int<lower=1> S;           // number of studies",
    "  int<lower=1> K;           // number of moderators",
    "  vector[S] y;              // observed effect sizes",
    "  vector<lower=0>[S] se;    // standard errors"
  )

  if (use_t_likelihood) {
    data_lines <- c(data_lines, "  vector<lower=0>[S] df;   // degrees of freedom")
  }

  data_lines <- c(data_lines, "  matrix[S, K] X;          // design matrix")

  data_block <- paste0("data {\n", paste(data_lines, collapse = "\n"), "\n}")

  # ---- PARAMETERS BLOCK ----
  tau_bounds <- if (is_re) emit_prior_bounds(priors$tau, default_lower = 0) else ""

  par_lines <- c(
    "  real mu;                  // intercept (pooled effect)",
    "  vector[K] beta;           // regression coefficients"
  )

  if (is_re) {
    par_lines <- c(par_lines,
                   glue::glue("  real{tau_bounds} tau;     // residual heterogeneity"),
                   "  vector[S] z;              // standardized random effects")
  }

  par_block <- paste0("parameters {\n", paste(par_lines, collapse = "\n"), "\n}")

  # ---- TRANSFORMED PARAMETERS BLOCK ----
  if (is_re) {
    tp_block <- "transformed parameters {
  vector[S] theta;          // true study effects
  theta = mu + X * beta + tau * z;
}"
  } else {
    tp_block <- "transformed parameters {
  vector[S] theta;          // true study effects
  theta = mu + X * beta;
}"
  }

  # ---- MODEL BLOCK ----
  model_lines <- character()

  # Priors
  mu_tgt <- emit_prior_target(priors$mu, "mu")
  if (nzchar(mu_tgt)) model_lines <- c(model_lines, paste0("  ", mu_tgt))

  # Beta priors
  unique_beta_priors <- unique(purrr::map_chr(priors$beta, format.bayesma_prior))
  if (length(unique_beta_priors) == 1) {
    first_prior <- priors$beta[[1]]
    beta_tgt <- emit_prior_target(first_prior, "beta")
    if (nzchar(beta_tgt)) model_lines <- c(model_lines, paste0("  ", beta_tgt))
  } else {
    for (k in seq_len(K)) {
      beta_tgt <- emit_prior_target(priors$beta[[k]], paste0("beta[", k, "]"))
      if (nzchar(beta_tgt)) model_lines <- c(model_lines, paste0("  ", beta_tgt))
    }
  }

  if (is_re) {
    tau_tgt <- emit_prior_target(priors$tau, "tau")
    if (nzchar(tau_tgt)) model_lines <- c(model_lines, paste0("  ", tau_tgt))
    model_lines <- c(model_lines, "  z ~ std_normal();")
  }

  # Likelihood
  if (use_t_likelihood) {
    model_lines <- c(model_lines,
                     "  for (i in 1:S) {",
                     "    target += student_t_lpdf(y[i] | df[i], theta[i], se[i]);",
                     "  }")
  } else {
    model_lines <- c(model_lines, "  y ~ normal(theta, se);")
  }

  model_block <- paste0("model {\n", paste(model_lines, collapse = "\n"), "\n}")

  # ---- GENERATED QUANTITIES BLOCK ----
  gq_lines <- character()

  if (is_re) {
    gq_lines <- c(gq_lines, "  real mu_new = normal_rng(mu, tau);")
  }

  if (use_t_likelihood) {
    gq_lines <- c(gq_lines,
                  "  vector[S] y_rep;",
                  "  for (i in 1:S) {",
                  "    y_rep[i] = student_t_rng(df[i], theta[i], se[i]);",
                  "  }")
  } else {
    gq_lines <- c(gq_lines,
                  "  vector[S] y_rep;",
                  "  for (i in 1:S) {",
                  "    y_rep[i] = normal_rng(theta[i], se[i]);",
                  "  }")
  }

  gq_block <- if (length(gq_lines) > 0) {
    paste0("generated quantities {\n", paste(gq_lines, collapse = "\n"), "\n}")
  } else {
    ""
  }

  # ---- ASSEMBLE ----
  parts <- c(data_block, par_block, tp_block, model_block)
  if (nzchar(gq_block)) parts <- c(parts, gq_block)
  paste(parts, collapse = "\n\n")
}
