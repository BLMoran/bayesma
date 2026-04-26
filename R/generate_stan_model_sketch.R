# ============================================================================
# Sketch: Composable Stan model generator for bayesma
#
# This file is a SKETCH / PROPOSAL — not ready to drop in as-is.
#
# PART 1: Block-builder pattern for the core generators (two-stage,
#          one-stage, meta-regression) — `generate_stan_model()`
#
# PART 2: Consolidated robma spec_component functions:
#          - `spec_component_pet_peese()`     (was 3 functions → 1)
#          - `spec_component_selection_weight()` (was 3 functions → 1)
#          - `spec_component_bias_corrected()`   (was 2 functions → 1)
# ============================================================================


# ---- Orchestrator ----

generate_stan_model <- function(stage = c("two_stage", "one_stage"),
                                model_type = c("random_effect", "common_effect"),
                                re_dist = c("normal", "t", "skew_normal", "mixture"),
                                priors,
                                likelihood = "gaussian",
                                use_t_likelihood = FALSE,
                                robust_config = list(enabled = FALSE),
                                multi_arm_config = list(enabled = FALSE),
                                null_model = FALSE,
                                mreg_config = list(enabled = FALSE)) {

  stage      <- rlang::arg_match(stage)
  model_type <- rlang::arg_match(model_type)
  re_dist    <- rlang::arg_match(re_dist)
  is_re      <- model_type == "random_effect"
  use_robust <- isTRUE(robust_config$enabled)
  use_multi_arm <- isTRUE(multi_arm_config$enabled)
  use_mreg   <- isTRUE(mreg_config$enabled)

  # Meta-regression currently only supports normal RE

  if (use_mreg && re_dist != "normal") {
    stop("Meta-regression currently only supports re_dist = 'normal'")
  }

  # Context object passed to all builders — avoids threading many args

  ctx <- list(
    stage           = stage,
    model_type      = model_type,
    re_dist         = re_dist,
    is_re           = is_re,
    priors          = priors,
    likelihood      = likelihood,
    use_t_likelihood = use_t_likelihood,
    use_robust      = use_robust,
    robust_config   = robust_config,
    use_multi_arm   = use_multi_arm,
    multi_arm_config = multi_arm_config,
    null_model      = null_model,
    use_mreg        = use_mreg,
    mreg_config     = mreg_config
  )

  # Build each block
  fn_block  <- build_functions_block(ctx)
  dat_block <- build_data_block(ctx)
  par_block <- build_parameters_block(ctx)
  tp_block  <- build_transformed_parameters_block(ctx)
  mod_block <- build_model_block(ctx)
  gq_block  <- build_generated_quantities_block(ctx)

  # Assemble — only include non-empty blocks
  parts <- c(fn_block, dat_block, par_block, tp_block, mod_block, gq_block)
  parts <- parts[nzchar(parts)]
  paste(parts, collapse = "\n\n")
}


# ---- Functions block ----

build_functions_block <- function(ctx) {
  if (ctx$re_dist == "skew_normal" && ctx$is_re) {
    "functions {
  real skew_normal_lpdf_custom(real x, real xi, real omega, real alpha) {
    real z = (x - xi) / omega;
    return normal_lpdf(z | 0, 1) - log(omega) + log(2) + normal_lcdf(alpha * z | 0, 1);
  }
}"
  } else {
    ""
  }
}


# ---- Data block ----

build_data_block <- function(ctx) {
  lines <- character()

  if (ctx$stage == "two_stage") {
    lines <- c(
      "  int<lower=1> S;",
      "  vector[S] y;",
      "  vector<lower=0>[S] se;"
    )
    if (ctx$use_t_likelihood) {
      lines <- c(lines, "  vector<lower=1>[S] df;")
    }
  } else {
    # one_stage
    lines <- c("  int<lower=1> N;", "  int<lower=1> S;")

    lines <- c(lines, switch(ctx$likelihood,
      binomial = c("  array[N] int<lower=0> events;",
                    "  array[N] int<lower=1> n;"),
      gaussian = c("  vector[N] y;",
                    "  vector<lower=0>[N] se;"),
      poisson  = c("  array[N] int<lower=0> events;",
                    "  vector<lower=0>[N] exposure;")
    ))

    lines <- c(lines,
      "  array[N] int<lower=0, upper=1> treat;",
      "  array[N] int<lower=1> study;"
    )

    if (ctx$use_multi_arm) {
      lines <- c(lines,
        "  // Multi-arm study structure",
        "  int<lower=1> n_ma_studies;",
        "  array[S] int<lower=1, upper=n_ma_studies> comp_to_ma;"
      )
    }
  }

  # Meta-regression: moderator count + design matrix
  if (ctx$use_mreg) {
    lines <- c(lines,
      "  int<lower=1> K;           // number of moderators",
      if (ctx$stage == "two_stage") {
        "  matrix[S, K] X;          // design matrix"
      } else {
        "  matrix[N, K] X;           // design matrix (arm-level)"
      }
    )
  }

  # Mixture data (both stages)
  if (ctx$re_dist == "mixture" && ctx$is_re) {
    lines <- c(lines,
      "  int<lower=2> K;",
      "  vector<lower=0>[K] prior_dirichlet_alpha;"
    )
  }

  paste0("data {\n", paste(lines, collapse = "\n"), "\n}")
}


# ---- Parameters block ----

build_parameters_block <- function(ctx) {
  p <- ctx$priors
  lines <- character()

  if (ctx$null_model && !ctx$is_re) {
    # Stan requires >= 1 parameter
    return("parameters {\n  real<lower=0, upper=0.001> dummy__;\n}")
  }

  # mu parameter (not present in null models)
  if (!ctx$null_model) {
    mu_bnds <- emit_prior_bounds(p$mu)
    if (ctx$stage == "one_stage") {
      gamma_bnds <- emit_prior_bounds(p$gamma)
      lines <- c(lines,
        glue::glue("  vector{gamma_bnds}[S] gamma;"),
        glue::glue("  real{mu_bnds} mu;")
      )
    } else {
      lines <- c(lines, glue::glue("  real{mu_bnds} mu;"))
    }
  }

  # Meta-regression coefficients
  if (ctx$use_mreg) {
    lines <- c(lines, "  vector[K] beta;           // regression coefficients")
  }

  # Random effects parameters
  if (ctx$is_re) {
    lines <- c(lines, build_re_parameters(ctx))
  }

  # Robust parameter
  if (ctx$use_robust) {
    rob_par <- emit_robust_parameters(ctx$robust_config)
    lines <- c(lines, rob_par$par)
  }

  paste0("parameters {\n", paste(lines, collapse = "\n"), "\n}")
}

# Sub-helper: random effects parameters
build_re_parameters <- function(ctx) {
  p <- ctx$priors
  tau_bnds <- emit_prior_bounds(p$tau, default_lower = 0)
  lines <- character()

  if (ctx$re_dist == "mixture") {
    # Mixture uses component-level parameters instead of tau + z
    if (ctx$stage == "one_stage") {
      lines <- c(lines, glue::glue("  real{tau_bnds} tau;"))
    }
    tk_bnds <- emit_prior_bounds(
      if (ctx$stage == "two_stage") p$mixture$tau_k else p$mixture$tau_k,
      default_lower = 0
    )
    # Two-stage uses mu_k (ordered); one-stage uses delta_k (ordered)
    comp_param <- if (ctx$stage == "two_stage") "mu_k" else "delta_k"
    lines <- c(lines,
      "  simplex[K] w;",
      glue::glue("  ordered[K] {comp_param};"),
      glue::glue("  vector{tk_bnds}[K] tau_k;")
    )
  } else {
    # normal / t / skew_normal all share tau + z
    lines <- c(lines, glue::glue("  real{tau_bnds} tau;"))

    # Multi-arm adds rho
    if (ctx$use_multi_arm) {
      rho_bnds <- emit_prior_bounds(p$rho, default_lower = -1, default_upper = 1)
      lines <- c(lines, glue::glue("  real{rho_bnds} rho;"))
    }

    # t adds nu
    if (ctx$re_dist == "t") {
      nu_bnds <- emit_prior_bounds(p$nu, default_lower = 2)
      lines <- c(lines, glue::glue("  real{nu_bnds} nu;"))
    }

    # skew_normal adds alpha_skew
    if (ctx$re_dist == "skew_normal") {
      lines <- c(lines, "  real alpha_skew;")
    }

    # All non-mixture RE models have z
    lines <- c(lines, "  vector[S] z;")
  }

  lines
}


# ---- Transformed parameters block ----

build_transformed_parameters_block <- function(ctx) {
  if (!ctx$is_re || ctx$re_dist == "mixture") return("")
  # Only normal/t/skew_normal have a transformed parameters block

  if (ctx$stage == "two_stage") {
    if (ctx$null_model) {
      "transformed parameters {\n  vector[S] theta = tau * z;\n}"
    } else if (ctx$use_mreg && ctx$is_re) {
      "transformed parameters {\n  vector[S] theta;\n  theta = mu + X * beta + tau * z;\n}"
    } else if (ctx$use_mreg) {
      "transformed parameters {\n  vector[S] theta;\n  theta = mu + X * beta;\n}"
    } else {
      "transformed parameters {\n  vector[S] theta = mu + tau * z;\n}"
    }
  } else {
    # one_stage
    if (ctx$use_multi_arm) {
      build_multi_arm_tp_block()
    } else {
      "transformed parameters {\n  vector[S] epsilon = tau * z;\n}"
    }
  }
}

build_multi_arm_tp_block <- function() {
  # This is long but structurally unique — keeping it as a standalone string
  # avoids fragile conditional assembly for complex Stan loop logic.
  "transformed parameters {
  vector[S] epsilon;
  {
    // Apply within-study correlation for multi-arm studies
    for (m in 1:n_ma_studies) {
      int n_arms_m = 0;
      for (s in 1:S) {
        if (comp_to_ma[s] == m) n_arms_m += 1;
      }

      if (n_arms_m == 1) {
        for (s in 1:S) {
          if (comp_to_ma[s] == m) {
            epsilon[s] = tau * z[s];
          }
        }
      } else {
        real sqrt_abs_rho = sqrt(fabs(rho));
        real sqrt_one_minus_abs_rho = sqrt(1.0 - fabs(rho));
        real sign_rho = rho >= 0 ? 1.0 : -1.0;

        real u_shared = 0.0;
        int first_found = 0;
        for (s in 1:S) {
          if (comp_to_ma[s] == m && first_found == 0) {
            u_shared = z[s];
            first_found = 1;
          }
        }

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
}


# ---- Model block ----

build_model_block <- function(ctx) {
  p <- ctx$priors

  # ---- Priors ----
  prior_lines <- character()

  if (ctx$null_model) {
    # Null model: no mu prior
  } else if (ctx$stage == "one_stage") {
    prior_lines <- c(
      emit_prior_target(p$gamma, "gamma"),
      emit_prior_target(p$mu, "mu")
    )
  } else {
    prior_lines <- emit_prior_target(p$mu, "mu")
  }

  # Meta-regression beta priors
  if (ctx$use_mreg) {
    prior_lines <- c(prior_lines, build_beta_priors(ctx))
  }

  if (ctx$use_robust) {
    prior_lines <- c(prior_lines, emit_robust_priors(ctx$robust_config))
  }

  if (ctx$is_re) {
    prior_lines <- c(prior_lines, build_re_priors(ctx))
  }

  # ---- Likelihood ----
  lik_code <- build_likelihood(ctx)

  paste0(
    "model {\n  // Priors\n  ",
    paste(prior_lines, collapse = "\n  "),
    "\n  // Likelihood\n",
    lik_code,
    "\n}"
  )
}

# Sub-helper: beta (meta-regression) priors
# Collapses to a single vectorised prior when all betas share the same prior
build_beta_priors <- function(ctx) {
  p <- ctx$priors
  unique_beta_priors <- unique(purrr::map_chr(p$beta, format.bayesma_prior))

  if (length(unique_beta_priors) == 1) {
    # All betas share the same prior — emit once as vector
    beta_tgt <- emit_prior_target(p$beta[[1]], "beta")
    if (nzchar(beta_tgt)) beta_tgt else character()
  } else {
    # Per-coefficient priors
    K <- length(p$beta)
    lines <- character()
    for (k in seq_len(K)) {
      beta_tgt <- emit_prior_target(p$beta[[k]], paste0("beta[", k, "]"))
      if (nzchar(beta_tgt)) lines <- c(lines, beta_tgt)
    }
    lines
  }
}

# Sub-helper: random effects priors
build_re_priors <- function(ctx) {
  p <- ctx$priors
  lines <- character()

  tau_tgt <- emit_prior_target(p$tau, "tau")
  lines <- c(lines, tau_tgt)

  if (ctx$use_multi_arm && !is.null(p$rho)) {
    lines <- c(lines, emit_prior_target(p$rho, "rho"))
  }

  switch(ctx$re_dist,
    normal = {
      lines <- c(lines, "target += std_normal_lpdf(z);")
    },
    t = {
      nu_tgt <- if (p$nu$family == "exponential") {
        glue::glue("target += exponential_lpdf(nu - 2 | {p$nu$rate});")
      } else {
        emit_prior_target(p$nu, "nu")
      }
      lines <- c(lines, nu_tgt,
        "for (i in 1:S)",
        "  target += student_t_lpdf(z[i] | nu, 0, 1);"
      )
    },
    skew_normal = {
      lines <- c(lines, emit_prior_target(p$alpha, "alpha_skew"),
        "for (i in 1:S)",
        "  target += skew_normal_lpdf_custom(z[i], 0, 1, alpha_skew);"
      )
    },
    mixture = {
      lines <- c(lines, "target += dirichlet_lpdf(w | prior_dirichlet_alpha);")
      comp_param <- if (ctx$stage == "two_stage") "mu_k" else "delta_k"
      comp_prior <- if (ctx$stage == "two_stage") p$mixture$mu_k else p$mixture$delta_k
      lines <- c(lines,
        "for (k in 1:K) {",
        paste0("  ", emit_prior_target(comp_prior, paste0(comp_param, "[k]"))),
        paste0("  ", emit_prior_target(p$mixture$tau_k, "tau_k[k]")),
        "}"
      )
    }
  )

  lines
}


# ---- Likelihood ----

build_likelihood <- function(ctx) {
  if (ctx$stage == "two_stage") {
    build_likelihood_two_stage(ctx)
  } else {
    build_likelihood_one_stage(ctx)
  }
}

build_likelihood_two_stage <- function(ctx) {
  # Two-stage is always Gaussian on summary statistics.
  # Variation comes from: t-likelihood, mixture RE, robust, null model.

  if (ctx$is_re && ctx$re_dist == "mixture") {
    return(build_likelihood_two_stage_mixture(ctx))
  }

  if (ctx$null_model) {
    loc <- if (ctx$is_re) "theta" else "rep_vector(0.0, S)"
  } else {
    loc <- if (ctx$is_re) "theta" else "rep_vector(mu, S)"
  }

  if (ctx$use_robust) {
    dist_fn <- if (ctx$use_t_likelihood) {
      "student_t_lpdf(y[i] | df[i], {loc}[i], se[i])"
    } else {
      "normal_lpdf(y[i] | {loc}[i], se[i])"
    }
    glue::glue("
  for (i in 1:S) {{
    real main_ll = {dist_fn};
    target += {emit_robust_likelihood('main_ll', 'y[i]', 'se[i]', ctx$robust_config)};
  }}")
  } else if (ctx$use_t_likelihood) {
    glue::glue("
  for (i in 1:S) {{
    target += student_t_lpdf(y[i] | df[i], {loc}[i], se[i]);
  }}")
  } else {
    glue::glue("\n  target += normal_lpdf(y | {loc}, se);")
  }
}

build_likelihood_two_stage_mixture <- function(ctx) {
  # Factored out because mixture likelihood is structurally different
  if (ctx$use_robust) {
    glue::glue("
  for (i in 1:S) {{
    vector[K] lps;
    for (k in 1:K)
      lps[k] = log(w[k]) + normal_lpdf(y[i] | mu_k[k], sqrt(square(tau_k[k]) + square(se[i])));
    real main_ll = log_sum_exp(lps);
    target += {emit_robust_likelihood('main_ll', 'y[i]', 'se[i]', ctx$robust_config)};
  }}")
  } else {
    "
  for (i in 1:S) {
    vector[K] lps;
    for (k in 1:K)
      lps[k] = log(w[k]) + normal_lpdf(y[i] | mu_k[k], sqrt(square(tau_k[k]) + square(se[i])));
    target += log_sum_exp(lps);
  }"
  }
}

build_likelihood_one_stage <- function(ctx) {
  # One-stage supports binomial/gaussian/poisson with optional robust wrapping.

  if (ctx$is_re && ctx$re_dist == "mixture") {
    return(build_likelihood_one_stage_mixture(ctx))
  }

  # Build eta expression — meta-regression adds X[i] * beta inside the
  # treatment-effect portion
  mu_part <- if (ctx$use_mreg) "mu + X[i] * beta" else "mu"
  re_part <- if (ctx$is_re) " + epsilon[study[i]]" else ""

  eta_term <- glue::glue(
    "gamma[study[i]] + ({mu_part}{re_part}) * treat[i]"
  )

  if (ctx$use_robust) {
    build_likelihood_one_stage_robust(ctx, eta_term)
  } else {
    # Vectorised standard likelihood
    lik_stmt <- switch(ctx$likelihood,
      binomial = "target += binomial_logit_lpmf(events | n, eta);",
      gaussian = "target += normal_lpdf(y | eta, se);",
      poisson  = "target += poisson_log_lpmf(events | log(exposure) + eta);"
    )
    paste0(
      "  {\n",
      "    vector[N] eta;\n",
      "    for (i in 1:N) {\n",
      "      eta[i] = ", eta_term, ";\n",
      "    }\n",
      "    ", lik_stmt, "\n",
      "  }"
    )
  }
}

build_likelihood_one_stage_robust <- function(ctx, eta_term) {
  # Each likelihood needs slightly different robust approximation.
  # Keeping these as explicit cases is clearer than trying to abstract further.
  rc <- ctx$robust_config

  outlier_dist <- if (is.finite(rc$df)) {
    glue::glue("student_t_lpdf(y_approx | {rc$df}, {rc$prior$mean}, se_approx * 3)")
  } else {
    glue::glue("normal_lpdf(y_approx | {rc$prior$mean}, se_approx * 3)")
  }

  switch(ctx$likelihood,
    gaussian = glue::glue("
  for (i in 1:N) {{
    real eta_i = {eta_term};
    real main_ll = normal_lpdf(y[i] | eta_i, se[i]);
    target += {emit_robust_likelihood('main_ll', 'y[i]', 'se[i]', rc)};
  }}"),
    binomial = glue::glue("
  for (i in 1:N) {{
    real eta_i = {eta_term};
    real main_ll = binomial_logit_lpmf(events[i:i] | n[i:i], rep_vector(eta_i, 1));
    real p_hat = inv_logit(eta_i);
    real se_approx = sqrt(1.0 / (n[i] * p_hat * (1 - p_hat) + 0.001));
    real y_approx = eta_i;
    real outlier_ll = {outlier_dist};
    target += log_sum_exp(log(pi_main) + main_ll, log1m(pi_main) + outlier_ll);
  }}"),
    poisson = glue::glue("
  for (i in 1:N) {{
    real eta_i = {eta_term};
    real main_ll = poisson_log_lpmf(events[i:i] | rep_vector(log(exposure[i]) + eta_i, 1));
    real rate_hat = exp(eta_i) * exposure[i];
    real se_approx = sqrt(1.0 / (rate_hat + 0.001));
    real y_approx = eta_i;
    real outlier_ll = {outlier_dist};
    target += log_sum_exp(log(pi_main) + main_ll, log1m(pi_main) + outlier_ll);
  }}")
  )
}

build_likelihood_one_stage_mixture <- function(ctx) {
  # Mixture one-stage — one function per likelihood (robust not supported here)
  switch(ctx$likelihood,
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
}


# ---- Generated quantities block ----

build_generated_quantities_block <- function(ctx) {
  lines <- character()

  # Pooled estimate
  if (ctx$null_model) {
    lines <- c(lines, "  real mu = 0.0;")
  } else if (ctx$stage == "one_stage") {
    lines <- c(lines, switch(ctx$likelihood,
      binomial = c("  real log_or = mu;", "  real or_pooled = exp(mu);"),
      gaussian = "  real pooled_diff = mu;",
      poisson  = c("  real log_rr = mu;", "  real rr_pooled = exp(mu);")
    ))
  } else {
    lines <- c(lines, "  real pooled = mu;")
  }

  # RE-specific quantities
  if (ctx$is_re) {
    lines <- c(lines, build_gq_re(ctx))
  }

  # Robust quantities
  if (ctx$use_robust) {
    lines <- c(lines, "  real pi_main_est = pi_main;")
    lines <- c(lines, build_gq_robust_outlier(ctx))
  }

  # Posterior predictive y_rep
  lines <- c(lines, build_gq_y_rep(ctx))

  # Prediction interval (mu_new) for RE models
  if (ctx$is_re) {
    lines <- c(lines, build_gq_mu_new(ctx))
  }

  paste0("generated quantities {\n", paste(lines, collapse = "\n"), "\n}")
}


# Sub-helpers for generated quantities

build_gq_re <- function(ctx) {
  lines <- character()

  if (ctx$stage == "two_stage") {
    # Two-stage: shrinkage estimates
    if (ctx$re_dist %in% c("normal", "t", "skew_normal")) {
      lines <- c(lines, "  vector[S] shrinkage = theta;")
    }
    if (ctx$re_dist == "t") {
      lines <- c(lines, "  real re_df = nu;")
    }
    if (ctx$re_dist == "mixture") {
      lines <- c(lines,
        "  vector[K] comp_means = mu_k;",
        "  vector[K] comp_taus = tau_k;",
        "  simplex[K] comp_weights = w;",
        build_gq_mixture_cluster("mu_k")
      )
    }
  } else {
    # One-stage: study-level effects
    if (ctx$re_dist %in% c("normal", "t", "skew_normal")) {
      lines <- c(lines, "  vector[S] study_effect = epsilon;")
    }
    if (ctx$re_dist == "t") {
      lines <- c(lines, "  real re_df = nu;")
    }
    if (ctx$re_dist == "mixture") {
      lines <- c(lines,
        "  vector[K] comp_deltas = delta_k;",
        "  vector[K] comp_taus = tau_k;",
        "  simplex[K] comp_weights = w;"
      )
    }
    if (ctx$use_multi_arm) {
      lines <- c(lines, "  real rho_est = rho;")
    }
  }

  if (ctx$re_dist == "skew_normal") {
    lines <- c(lines,
      "  real re_alpha = alpha_skew;",
      "  real delta_sk = alpha_skew / sqrt(1 + square(alpha_skew));",
      "  real b_sk = sqrt(2.0 / pi());",
      "  real re_skewness = (4.0 - pi()) / 2.0 * pow(b_sk * delta_sk, 3) /",
      "                     pow(1 - square(b_sk * delta_sk), 1.5);"
    )
  }

  lines
}

build_gq_mixture_cluster <- function(comp_mean_param) {
  # Cluster assignment for two-stage mixture
  c(
    "  array[S] int<lower=1, upper=K> cluster;",
    "  for (i in 1:S) {",
    "    vector[K] lps;",
    "    for (k in 1:K)",
    glue::glue("      lps[k] = log(w[k]) + normal_lpdf(y[i] | {comp_mean_param}[k], sqrt(square(tau_k[k]) + square(se[i])));"),
    "    cluster[i] = categorical_rng(softmax(lps));",
    "  }"
  )
}

build_gq_robust_outlier <- function(ctx) {
  # Only for two-stage non-mixture RE models
  if (ctx$stage != "two_stage") return(character())
  if (!ctx$is_re || !ctx$re_dist %in% c("normal", "t", "skew_normal")) return(character())

  rc <- ctx$robust_config
  outlier_ll <- if (is.finite(rc$df)) {
    glue::glue("student_t_lpdf(y[i] | {rc$df}, {rc$prior$mean}, se[i] * 3)")
  } else {
    glue::glue("normal_lpdf(y[i] | {rc$prior$mean}, se[i] * 3)")
  }

  c(
    "  vector[S] prob_outlier;",
    "  for (i in 1:S) {",
    "    real ll_main = normal_lpdf(y[i] | theta[i], se[i]);",
    paste0("    real ll_out = ", as.character(outlier_ll), ";"),
    "    real lp_main = log(pi_main) + ll_main;",
    "    real lp_out = log1m(pi_main) + ll_out;",
    "    prob_outlier[i] = exp(lp_out - log_sum_exp(lp_main, lp_out));",
    "  }"
  )
}

build_gq_y_rep <- function(ctx) {
  if (ctx$stage == "two_stage") {
    return(build_gq_y_rep_two_stage(ctx))
  }

  # One-stage: need eta computation (mirrors the model block logic)
  mu_part <- if (ctx$use_mreg) "mu + X[i] * beta" else "mu"
  re_part <- if (ctx$is_re && ctx$re_dist %in% c("normal", "t", "skew_normal")) {
    " + epsilon[study[i]]"
  } else {
    ""
  }
  eta_re <- glue::glue("gamma[study[i]] + ({mu_part}{re_part}) * treat[i]")

  switch(ctx$likelihood,
    gaussian = c(
      "  vector[N] y_rep;",
      "  for (i in 1:N) {",
      glue::glue("    real eta_i = {eta_re};"),
      "    y_rep[i] = normal_rng(eta_i, se[i]);",
      "  }"
    ),
    binomial = c(
      "  vector[N] y_rep;",
      "  for (i in 1:N) {",
      glue::glue("    real eta_i = {eta_re};"),
      "    y_rep[i] = binomial_rng(n[i], inv_logit(eta_i)) * 1.0;",
      "  }"
    ),
    poisson = c(
      "  vector[N] y_rep;",
      "  for (i in 1:N) {",
      glue::glue("    real eta_i = {eta_re};"),
      "    y_rep[i] = poisson_log_rng(log(exposure[i]) + eta_i) * 1.0;",
      "  }"
    )
  )
}

build_gq_y_rep_two_stage <- function(ctx) {
  if (ctx$use_mreg) {
    # Meta-regression: theta already includes mu + X*beta (+ tau*z if RE)
    loc <- "theta[i]"
  } else if (ctx$is_re && ctx$re_dist %in% c("normal", "t", "skew_normal")) {
    loc <- "theta[i]"
  } else if (ctx$is_re && ctx$re_dist == "mixture") {
    # Mixture uses cluster assignment
    return(c(
      "  vector[S] y_rep;",
      "  for (i in 1:S) {",
      "    int comp = cluster[i];",
      "    y_rep[i] = normal_rng(mu_k[comp], sqrt(square(tau_k[comp]) + square(se[i])));",
      "  }"
    ))
  } else if (ctx$null_model) {
    loc <- "0.0"
  } else {
    loc <- "mu"
  }

  c(
    "  vector[S] y_rep;",
    "  for (i in 1:S)",
    glue::glue("    y_rep[i] = normal_rng({loc}, se[i]);")
  )
}

build_gq_mu_new <- function(ctx) {
  switch(ctx$re_dist,
    normal = "  real mu_new = normal_rng(mu, tau);",
    t      = "  real mu_new = mu + tau * student_t_rng(nu, 0, 1);",
    skew_normal = c(
      "  real mu_new;",
      "  {",
      "    real u0 = normal_rng(0, 1);",
      "    real u1 = normal_rng(0, 1);",
      "    real dg = alpha_skew / sqrt(1 + square(alpha_skew));",
      "    real z_sn = dg * fabs(u0) + sqrt(1 - square(dg)) * u1;",
      "    mu_new = mu + tau * z_sn;",
      "  }"
    ),
    mixture = {
      comp_param <- if (ctx$stage == "two_stage") "mu_k" else "delta_k"
      c(
        "  real mu_new;",
        "  {",
        "    int comp = categorical_rng(w);",
        glue::glue("    mu_new = normal_rng({comp_param}[comp], tau_k[comp]);"),
        "  }"
      )
    }
  )
}


# ============================================================================
# PART 2: Consolidated robma spec_component functions
#
# These are standalone generators for the robma model-averaging framework.
# Each returns list(stan_code, stan_data). They share the two-axis pattern
# of has_re (RE vs FE) × null_model (mu estimated vs mu = 0).
# ============================================================================


# ============================================================================
# PET-PEESE component (consolidates 3 former functions)
#
# Replaces:
#   spec_component_pet_peese()     — H1, bias_transform = "inv_sqrt_n"
#   spec_component_pet_peese_h0()  — H0, bias_transform = "inv_sqrt_n"
#   spec_component_peese_h0()      — H0, bias_transform = "inv_n"
# ============================================================================

spec_component_pet_peese <- function(es, S, n_c, n_i, priors,
                                     has_re,
                                     null_model = FALSE,
                                     bias_transform = c("inv_sqrt_n", "inv_n")) {

  bias_transform <- rlang::arg_match(bias_transform)

  # ---- Priors ----
  # For H1 models, priors come from the priors list.
  # For H0 models, use sensible defaults (matching current hardcoded values).
  if (null_model) {
    mu_tgt   <- NULL
    tau_tgt  <- "target += student_t_lpdf(tau | 3, 0, 2.5)\n          - student_t_lccdf(0 | 3, 0, 2.5);"
    beta_tgt <- "target += normal_lpdf(beta_bias | 0, 5);"
    tau_bnds <- "<lower=0>"
  } else {
    mu_tgt   <- emit_prior_target(priors$mu, "mu")
    sp       <- priors$selection %||% list()
    if (is.null(sp$beta_bias)) sp$beta_bias <- normal(0, 5)
    beta_tgt <- emit_prior_target(sp$beta_bias, "beta_bias")
    if (has_re) {
      tau_tgt  <- emit_prior_target(priors$tau, "tau")
      tau_bnds <- emit_prior_bounds(priors$tau, default_lower = 0)
    }
  }

  # ---- Transformed data variable name and computation ----
  td_var  <- if (bias_transform == "inv_sqrt_n") "inv_sqrt_n" else "inv_n"
  td_expr <- if (bias_transform == "inv_sqrt_n") {
    "1.0 / sqrt(n_total[i])"
  } else {
    "1.0 / n_total[i]"
  }

  # ---- Build Stan code ----

  # Data block
  data_lines <- c(
    "  int<lower=1> N;",
    "  vector[N] y;",
    "  vector<lower=0>[N] se;",
    "  vector<lower=0>[N] n_total;"
  )
  data_block <- paste0("data {\n", paste(data_lines, collapse = "\n"), "\n}")

  # Transformed data block
  td_block <- paste0(
    "transformed data {\n",
    "  vector[N] ", td_var, ";\n",
    "  for (i in 1:N)\n",
    "    ", td_var, "[i] = ", td_expr, ";\n",
    "}"
  )

  # Parameters block
  par_lines <- character()
  if (!null_model) {
    par_lines <- c(par_lines, "  real mu;")
  }
  if (has_re) {
    par_lines <- c(par_lines, glue::glue("  real{tau_bnds} tau;"))
  }
  par_lines <- c(par_lines, "  real beta_bias;")
  par_block <- paste0("parameters {\n", paste(par_lines, collapse = "\n"), "\n}")

  # Model block
  model_lines <- character()
  if (!null_model) model_lines <- c(model_lines, paste0("  ", mu_tgt))
  if (has_re)      model_lines <- c(model_lines, paste0("  ", tau_tgt))
  model_lines <- c(model_lines, paste0("  ", beta_tgt))

  # Location term: mu (or 0) + beta_bias * bias_covariate
  mu_term <- if (null_model) "0.0" else "mu"
  loc_expr <- paste0(mu_term, " + beta_bias * ", td_var, "[i]")

  if (has_re) {
    model_lines <- c(model_lines,
      "  for (i in 1:N) {",
      "    real sigma_i = sqrt(square(tau) + square(se[i]));",
      glue::glue("    target += normal_lpdf(y[i] | {loc_expr}, sigma_i);"),
      "  }"
    )
  } else {
    model_lines <- c(model_lines,
      "  for (i in 1:N)",
      glue::glue("    target += normal_lpdf(y[i] | {loc_expr}, se[i]);")
    )
  }
  model_block <- paste0("model {\n", paste(model_lines, collapse = "\n"), "\n}")

  # Generated quantities block
  gq_lines <- character()
  if (null_model) {
    gq_lines <- c(gq_lines, "  real mu = 0.0;")
  } else {
    gq_lines <- c(gq_lines, "  real pooled = mu;")
  }
  gq_lines <- c(gq_lines, "  real bias_slope = beta_bias;")
  if (has_re && !null_model) {
    gq_lines <- c(gq_lines, "  real mu_new = normal_rng(mu, tau);")
  }
  gq_block <- paste0("generated quantities {\n",
                      paste(gq_lines, collapse = "\n"), "\n}")

  # Assemble
  stan_code <- paste(data_block, td_block, par_block, model_block, gq_block,
                     sep = "\n\n")

  stan_data <- list(N = S, y = es$yi, se = es$sei,
                    n_total = as.array(n_c + n_i))

  list(stan_code = stan_code, stan_data = stan_data)
}


# ============================================================================
# Selection weight component (consolidates 3 former functions)
#
# Replaces:
#   spec_component_selection_weight()     — H1 + RE
#   spec_component_selection_weight_fe()  — H1 + FE
#   spec_component_selection_weight_h0()  — H0 + RE or FE
#
# The normalisation constant in the weight-function likelihood has a general
# form across all 4 variants (has_re × null_model):
#
#   Phi((z_bounds[k] * se[i] - mu) / scale) - Phi((z_bounds[k+1] * se[i] - mu) / scale)
#
# where mu = 0 for null models, and scale = sigma[i] (RE) or se[i] (FE).
# The original code simplified these algebraically in some cases (e.g.,
# FE+H0 reduced to Phi(z_bounds[k]) - Phi(z_bounds[k+1])). The general
# form is mathematically equivalent and Stan computes the same result.
# ============================================================================

spec_component_selection_weight <- function(es, S, priors = NULL, p_cuts,
                                            has_re,
                                            null_model = FALSE) {

  K <- length(p_cuts) + 1L

  # ---- Priors (only used for H1 models) ----
  if (!null_model) {
    mu_tgt <- emit_prior_target(priors$mu, "mu")
    if (has_re) {
      tau_tgt  <- emit_prior_target(priors$tau, "tau")
      tau_bnds <- emit_prior_bounds(priors$tau, default_lower = 0)
    }
  } else if (has_re) {
    # H0 + RE: hardcoded half-t prior (matching current code)
    tau_tgt <- "target += student_t_lpdf(tau | 3, 0, 2.5)\n          - student_t_lccdf(0 | 3, 0, 2.5);"
    tau_bnds <- "<lower=0>"
  }

  # ---- Expressions that vary by has_re × null_model ----
  scale_expr <- if (has_re) "sigma[i]" else "se[i]"
  mu_expr    <- if (null_model) "0.0" else "mu"

  # Normalisation constant: Phi argument for z_bounds[k]
  # General form: (z_bounds[k] * se[i] - mu) / scale
  phi_arg_k   <- glue::glue("(z_bounds[k] * se[i] - {mu_expr}) / {scale_expr}")
  phi_arg_k1  <- glue::glue("(z_bounds[k+1] * se[i] - {mu_expr}) / {scale_expr}")

  # ---- Data block ----
  data_block <- "data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  int<lower=2> K;
  vector[K-1] p_cutoffs;
}"

  # ---- Transformed data block ----
  td_block <- "transformed data {
  vector[K+1] z_bounds;
  z_bounds[1] = 8.0;
  z_bounds[K+1] = -8.0;
  for (k in 1:(K-1))
    z_bounds[k+1] = inv_Phi(1.0 - p_cutoffs[k]);
}"

  # ---- Parameters block ----
  par_lines <- character()
  if (!null_model) {
    par_lines <- c(par_lines, "  real mu;")
  }
  if (has_re) {
    par_lines <- c(par_lines, glue::glue("  real{tau_bnds} tau;"))
  }
  par_lines <- c(par_lines, "  vector<lower=0.01, upper=0.99>[K-1] omega_raw;")
  par_block <- paste0("parameters {\n", paste(par_lines, collapse = "\n"), "\n}")

  # ---- Transformed parameters block ----
  tp_lines <- c(
    "  vector[K] omega;",
    if (has_re) "  vector[N] sigma;",
    "  omega[1] = 1.0;",
    "  for (k in 1:(K-1))",
    "    omega[k+1] = omega_raw[k];"
  )
  if (has_re) {
    tp_lines <- c(tp_lines,
      "  for (i in 1:N)",
      "    sigma[i] = sqrt(square(tau) + square(se[i]));"
    )
  }
  tp_block <- paste0("transformed parameters {\n",
                      paste(tp_lines, collapse = "\n"), "\n}")

  # ---- Model block ----
  model_lines <- character()

  if (!null_model) model_lines <- c(model_lines, paste0("  ", mu_tgt))
  if (has_re)      model_lines <- c(model_lines, paste0("  ", tau_tgt))

  model_lines <- c(model_lines,
    "  for (k in 1:(K-1))",
    "    target += beta_lpdf(omega_raw[k] | 1, 1);"
  )

  # Weight function likelihood — the core logic shared by all variants.
  # The weight assignment (log_w) is identical across all 4; only the
  # likelihood scale and normalisation constant Phi arguments differ.
  model_lines <- c(model_lines, as.character(glue::glue("  for (i in 1:N) {{
    real z_i = y[i] / se[i];
    real p_val = 1.0 - Phi(z_i);
    real log_w;
    if (p_val < p_cutoffs[1]) {{
      log_w = 0.0;
    }} else if (K == 2) {{
      log_w = log(omega[2]);
    }} else {{
      log_w = log(omega[K]);
      for (k in 1:(K-2)) {{
        if (p_val >= p_cutoffs[k] && p_val < p_cutoffs[k+1])
          log_w = log(omega[k+1]);
      }}
    }}
    target += normal_lpdf(y[i] | {mu_expr}, {scale_expr}) + log_w;
    {{
      real norm_c = 0;
      for (k in 1:K) {{
        real prob_k = Phi({phi_arg_k})
                    - Phi({phi_arg_k1});
        norm_c += omega[k] * fmax(prob_k, 1e-15);
      }}
      target += -log(fmax(norm_c, 1e-15));
    }}
  }}")))

  model_block <- paste0("model {\n", paste(model_lines, collapse = "\n"), "\n}")

  # ---- Generated quantities block ----
  gq_lines <- character()
  if (null_model) {
    gq_lines <- c(gq_lines, "  real mu = 0.0;")
  } else {
    gq_lines <- c(gq_lines, "  real pooled = mu;")
  }
  gq_lines <- c(gq_lines, "  vector[K] weights = omega;")
  if (has_re && !null_model) {
    gq_lines <- c(gq_lines, "  real mu_new = normal_rng(mu, tau);")
  }
  gq_block <- paste0("generated quantities {\n",
                      paste(gq_lines, collapse = "\n"), "\n}")

  # Assemble
  stan_code <- paste(data_block, td_block, par_block, tp_block,
                     model_block, gq_block, sep = "\n\n")

  stan_data <- list(N = S, y = es$yi, se = es$sei, K = K, p_cutoffs = p_cuts)

  list(stan_code = stan_code, stan_data = stan_data)
}


# ============================================================================
# Bias-corrected component (consolidates 2 former functions)
#
# Replaces:
#   spec_component_bias_corrected()  — H1, has_re delegates to
#                                      generate_stan_code_bias_corrected()
#                                      for the full RE model; FE is simpler.
#   spec_component_jung_h0()         — H0, mu = 0, bias shift B active
#
# The RE + H1 variant (generate_stan_code_bias_corrected) is substantially
# more complex than all other variants — it uses inverse-variance
# parameterisation, per-study bias weights (w_bias_raw), and separate
# biased/unbiased theta vectors. Rather than inlining that complexity here,
# the has_re + H1 case delegates to generate_stan_code_bias_corrected() as
# the original code does, keeping this function focused on the three simpler
# variants: FE+H1, RE+H0, FE+H0.
# ============================================================================

spec_component_bias_corrected <- function(es, S, priors,
                                          has_re,
                                          null_model = FALSE) {

  # ---- RE + H1: delegate to the full model (unchanged) ----
  if (has_re && !null_model) {
    stan_code <- generate_stan_code_bias_corrected(priors)
    stan_data <- list(N = S, y = es$yi, se_y = es$sei,
                      use_known_bias = 0L, known_bias = rep(0L, S))
    return(list(stan_code = stan_code, stan_data = stan_data))
  }

  # ---- Priors ----
  if (null_model) {
    # H0 models use hardcoded priors (matching current code)
    tau_tgt <- "target += student_t_lpdf(tau | 3, 0, 2.5)\n          - student_t_lccdf(0 | 3, 0, 2.5);"
    b_tgt   <- "target += uniform_lpdf(B | 0, 2);"
    pb_tgt  <- "target += beta_lpdf(p_bias | 1, 1);"
  } else {
    # FE + H1
    mu_tgt   <- emit_prior_target(priors$mu, "mu")
    b_prior  <- priors$bias %||% priors$b %||% uniform(0, 2)
    b_tgt    <- emit_prior_target(b_prior, "B")
    b_bnds   <- emit_prior_bounds(b_prior, default_lower = 0)
    pb_prior <- priors$p_bias %||% beta(1, 1)
    pb_tgt   <- emit_prior_target(pb_prior, "p_bias")
  }

  # ---- Location terms ----
  mu_term <- if (null_model) "0.0" else "mu"

  # ---- Data block ----
  # FE+H1 uses se_y and has known_bias; H0 variants use se and are simpler
  if (null_model) {
    data_block <- "data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
}"
  } else {
    data_block <- "data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se_y;
  int<lower=0, upper=1> use_known_bias;
  array[N] int<lower=0, upper=1> known_bias;
}"
  }

  # ---- Parameters block ----
  par_lines <- character()
  if (!null_model) {
    par_lines <- c(par_lines, "  real mu;",
                   glue::glue("  real{b_bnds} B;"))
  } else {
    par_lines <- c(par_lines, "  real<lower=0> B;")
  }
  par_lines <- c(par_lines, "  real<lower=0, upper=1> p_bias;")
  if (has_re) {
    if (null_model) {
      par_lines <- c(par_lines, "  real<lower=0> tau;")
    }
  }
  par_block <- paste0("parameters {\n", paste(par_lines, collapse = "\n"), "\n}")

  # ---- Model block ----
  model_lines <- character()
  if (!null_model) model_lines <- c(model_lines, paste0("  ", mu_tgt))
  model_lines <- c(model_lines, paste0("  ", b_tgt), paste0("  ", pb_tgt))
  if (has_re) model_lines <- c(model_lines, paste0("  ", tau_tgt))

  # SE variable name differs between H1 (se_y) and H0 (se)
  se_var <- if (null_model) "se" else "se_y"

  # Scale: sigma_i (with tau) or just se
  if (has_re) {
    scale_line <- glue::glue("    real sigma_i = sqrt(square(tau) + square({se_var}[i]));")
    scale_expr <- "sigma_i"
  } else {
    scale_line <- NULL
    scale_expr <- paste0(se_var, "[i]")
  }

  if (null_model) {
    # H0: mixture of unbiased (mu=0) and biased (mu=B) components
    model_lines <- c(model_lines, "  for (i in 1:N) {")
    if (!is.null(scale_line)) model_lines <- c(model_lines, paste0("  ", scale_line))
    model_lines <- c(model_lines,
      glue::glue("    real lp_u = log1m(p_bias) + normal_lpdf(y[i] | {mu_term}, {scale_expr});"),
      glue::glue("    real lp_b = log(p_bias + 1e-15) + normal_lpdf(y[i] | B, {scale_expr});"),
      "    target += log_sum_exp(lp_u, lp_b);",
      "  }"
    )
  } else {
    # FE + H1: supports known_bias indicator
    model_lines <- c(model_lines,
      "  for (i in 1:N) {",
      "    if (use_known_bias == 1) {",
      "      if (known_bias[i] == 1) {",
      glue::glue("        target += normal_lpdf(y[i] | {mu_term} + B, {scale_expr});"),
      "      } else {",
      glue::glue("        target += normal_lpdf(y[i] | {mu_term}, {scale_expr});"),
      "      }",
      "    } else {",
      glue::glue("      real lp_u = log1m(p_bias) + normal_lpdf(y[i] | {mu_term}, {scale_expr});"),
      glue::glue("      real lp_b = log(p_bias + 1e-15)"),
      glue::glue("                + normal_lpdf(y[i] | {mu_term} + B, {scale_expr});"),
      "      target += log_sum_exp(lp_u, lp_b);",
      "    }",
      "  }"
    )
  }
  model_block <- paste0("model {\n", paste(model_lines, collapse = "\n"), "\n}")

  # ---- Generated quantities block ----
  gq_lines <- character()
  if (null_model) {
    gq_lines <- c(gq_lines, "  real mu = 0.0;")
  } else {
    gq_lines <- c(gq_lines, "  real pooled = mu;")
  }
  gq_block <- paste0("generated quantities {\n",
                      paste(gq_lines, collapse = "\n"), "\n}")

  # Assemble
  stan_code <- paste(data_block, par_block, model_block, gq_block,
                     sep = "\n\n")

  # Stan data
  if (null_model) {
    stan_data <- list(N = S, y = es$yi, se = es$sei)
  } else {
    stan_data <- list(N = S, y = es$yi, se_y = es$sei,
                      use_known_bias = 0L, known_bias = rep(0L, S))
  }

  list(stan_code = stan_code, stan_data = stan_data)
}

