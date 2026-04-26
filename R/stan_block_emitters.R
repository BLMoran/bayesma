# Stan Block Emitter Functions
#
# Composable helper functions for generating Stan code blocks.
# Used by generate_stan_code_one_stage() and generate_stan_code_two_stage()
# to eliminate duplicated switch(re_dist) logic.


# ---- Functions block ----

#' Emit the Stan functions block (skew_normal custom lpdf if needed)
#' @param re_dist Character. Random effect distribution type.
#' @param is_re Logical. TRUE if random effects model.
#' @return Character. Stan functions block or "".
#' @noRd
emit_functions_block <- function(re_dist, is_re) {
  if (re_dist == "skew_normal" && is_re) {
    "functions {
  real skew_normal_lpdf_custom(real x, real xi, real omega, real alpha) {
    real z = (x - xi) / omega;
    return normal_lpdf(z | 0, 1) - log(omega) + log(2) + normal_lcdf(alpha * z | 0, 1);
  }
}"
  } else ""
}


# ---- RE parameter declarations ----

#' Emit random-effects parameter declarations for any stage
#'
#' Returns the Stan parameter lines to append inside the parameters block.
#' Detects one-stage (delta_k) vs two-stage (mu_k) mixture parameterisation
#' from the priors structure.
#'
#' @param re_dist Character. RE distribution: "normal", "t", "skew_normal", "mixture".
#' @param priors List. Prior objects (from resolve_priors).
#' @param include_rho Logical. Whether to include rho parameter (multi-arm).
#' @return Character. Stan parameter lines (no braces).
#' @noRd
emit_re_parameters <- function(re_dist, priors, include_rho = FALSE) {
  p <- priors
  tau_bnds <- emit_prior_bounds(p$tau, default_lower = 0)

  rho_line <- ""
  if (include_rho) {
    rho_bnds <- emit_prior_bounds(p$rho, default_lower = -1, default_upper = 1)
    rho_line <- as.character(glue::glue("\n  real{rho_bnds} rho;"))
  }

  as.character(switch(re_dist,
    normal = glue::glue(
      "\n  real{tau_bnds} tau;{rho_line}\n  vector[S] z;"),
    t = {
      nu_bnds <- emit_prior_bounds(p$nu, default_lower = 2)
      glue::glue(
        "\n  real{tau_bnds} tau;{rho_line}\n  real{nu_bnds} nu;\n  vector[S] z;")
    },
    skew_normal = glue::glue(
      "\n  real{tau_bnds} tau;{rho_line}\n  real alpha_skew;\n  vector[S] z;"),
    mixture = {
      tk_bnds <- emit_prior_bounds(p$mixture$tau_k, default_lower = 0)
      if (!is.null(p$mixture$delta_k)) {
        # One-stage mixture: has global tau + delta_k
        glue::glue(
          "\n  real{tau_bnds} tau;\n  simplex[K] w;\n  ordered[K] delta_k;\n  vector{tk_bnds}[K] tau_k;")
      } else {
        # Two-stage mixture: no global tau, uses mu_k
        glue::glue(
          "\n  simplex[K] w;\n  ordered[K] mu_k;\n  vector{tk_bnds}[K] tau_k;")
      }
    }
  ))
}


# ---- RE prior statements ----

#' Emit random-effects prior target statements for the model block
#'
#' Returns Stan model block code for RE distribution priors.
#' Detects one-stage vs two-stage mixture from priors structure.
#'
#' @param re_dist Character. RE distribution type.
#' @param priors List. Prior objects.
#' @param include_rho Logical. Include rho prior statement.
#' @return Character. Stan model block lines for RE priors.
#' @noRd
emit_re_priors <- function(re_dist, priors, include_rho = FALSE) {
  p <- priors
  tau_tgt <- emit_prior_target(p$tau, "tau")

  rho_prior_code <- ""
  if (include_rho && !is.null(p$rho)) {
    rho_prior_code <- paste0("\n  ", emit_prior_target(p$rho, "rho"))
  }

  as.character(switch(re_dist,
    normal = glue::glue(
      "\n  {tau_tgt}{rho_prior_code}\n  target += std_normal_lpdf(z);"),
    t = {
      nu_tgt <- if (p$nu$family == "exponential") {
        glue::glue(
          "target += exponential_lpdf(nu - 2 | {p$nu$rate});")
      } else emit_prior_target(p$nu, "nu")
      glue::glue(
        "\n  {tau_tgt}{rho_prior_code}\n  {nu_tgt}\n  for (i in 1:S)\n    target += student_t_lpdf(z[i] | nu, 0, 1);")
    },
    skew_normal = {
      a_tgt <- emit_prior_target(p$alpha, "alpha_skew")
      glue::glue(
        "\n  {tau_tgt}{rho_prior_code}\n  {a_tgt}\n  for (i in 1:S)\n    target += skew_normal_lpdf_custom(z[i], 0, 1, alpha_skew);")
    },
    mixture = {
      if (!is.null(p$mixture$delta_k)) {
        # One-stage mixture: has global tau
        dk_tgt <- emit_prior_target(p$mixture$delta_k, "delta_k[k]")
        tk_tgt <- emit_prior_target(p$mixture$tau_k, "tau_k[k]")
        glue::glue(
          "\n  {tau_tgt}\n  target += dirichlet_lpdf(w | prior_dirichlet_alpha);\n  for (k in 1:K) {{\n    {dk_tgt}\n    {tk_tgt}\n  }}")
      } else {
        # Two-stage mixture: no global tau
        mk_tgt <- emit_prior_target(p$mixture$mu_k, "mu_k[k]")
        tk_tgt <- emit_prior_target(p$mixture$tau_k, "tau_k[k]")
        glue::glue(
          "\n  target += dirichlet_lpdf(w | prior_dirichlet_alpha);\n  for (k in 1:K) {{\n    {mk_tgt}\n    {tk_tgt}\n  }}")
      }
    }
  ))
}


# ---- RE extra generated quantities ----

#' Emit RE-specific generated quantities (study effects, distribution extras)
#'
#' Returns Stan generated quantities code for study-level effects, degrees of
#' freedom, skewness, or mixture component summaries. Uses stage to select
#' appropriate variable names (epsilon/study_effect vs theta/shrinkage).
#'
#' For two-stage mixture, also includes the cluster assignment loop.
#'
#' @param re_dist Character. RE distribution type.
#' @param stage Character. "one_stage" or "two_stage".
#' @return Character. Stan generated quantities code.
#' @noRd
emit_re_extra_gq <- function(re_dist, stage) {
  if (stage == "two_stage") {
    switch(re_dist,
      normal = "\n  vector[S] shrinkage = theta;",
      t = paste0(
        "\n  vector[S] shrinkage = theta;",
        "\n  real re_df = nu;"),
      skew_normal = "
  vector[S] shrinkage = theta;
  real re_alpha = alpha_skew;
  real delta_sk = alpha_skew / sqrt(1 + square(alpha_skew));
  real b_sk = sqrt(2.0 / pi());
  real re_skewness = (4.0 - pi()) / 2.0 * pow(b_sk * delta_sk, 3) /
                     pow(1 - square(b_sk * delta_sk), 1.5);",
      mixture = "
  vector[K] comp_means = mu_k;
  vector[K] comp_taus = tau_k;
  simplex[K] comp_weights = w;
  array[S] int<lower=1, upper=K> cluster;
  for (i in 1:S) {
    vector[K] lps;
    for (k in 1:K)
      lps[k] = log(w[k]) + normal_lpdf(y[i] | mu_k[k], sqrt(square(tau_k[k]) + square(se[i])));
    cluster[i] = categorical_rng(softmax(lps));
  }"
    )
  } else {
    # one_stage
    switch(re_dist,
      normal = "\n  vector[S] study_effect = epsilon;",
      t = "\n  vector[S] study_effect = epsilon;\n  real re_df = nu;",
      skew_normal = "
  vector[S] study_effect = epsilon;
  real re_alpha = alpha_skew;
  real delta_sk = alpha_skew / sqrt(1 + square(alpha_skew));
  real b_sk = sqrt(2.0 / pi());
  real re_skewness = (4.0 - pi()) / 2.0 * pow(b_sk * delta_sk, 3) /
                     pow(1 - square(b_sk * delta_sk), 1.5);",
      mixture = "\n  vector[K] comp_deltas = delta_k;\n  vector[K] comp_taus = tau_k;\n  simplex[K] comp_weights = w;"
    )
  }
}


# ---- Prediction interval generated quantity ----

#' Emit mu_new prediction interval in generated quantities
#'
#' Generates Stan code for the predictive distribution of a new study effect.
#' Detects one-stage vs two-stage mixture from priors structure to use
#' delta_k (one-stage) vs mu_k (two-stage).
#'
#' @param re_dist Character. RE distribution type.
#' @param priors List. Prior objects.
#' @return Character. Stan generated quantities code for mu_new.
#' @noRd
emit_mu_new_gq <- function(re_dist, priors) {
  is_two_stage_mixture <- !is.null(priors$mixture$mu_k)

  switch(re_dist,
    normal = "\n  real mu_new = normal_rng(mu, tau);",
    t = "\n  real mu_new = mu + tau * student_t_rng(nu, 0, 1);",
    skew_normal = "
  real mu_new;
  {
    real u0 = normal_rng(0, 1);
    real u1 = normal_rng(0, 1);
    real dg = alpha_skew / sqrt(1 + square(alpha_skew));
    real z_sn = dg * fabs(u0) + sqrt(1 - square(dg)) * u1;
    mu_new = mu + tau * z_sn;
  }",
    mixture = {
      comp_loc <- if (is_two_stage_mixture) "mu_k" else "delta_k"
      paste0("
  real mu_new;
  {
    int comp = categorical_rng(w);
    mu_new = normal_rng(", comp_loc, "[comp], tau_k[comp]);
  }")
    }
  )
}


# ---- Stan program assembly ----

#' Assemble a complete Stan program from its blocks
#'
#' Concatenates non-empty Stan blocks with double-newline separators.
#'
#' @param fn_block Character. Functions block (may be "").
#' @param data_block Character. Data block.
#' @param par_block Character. Parameters block.
#' @param tp_block Character. Transformed parameters block (may be "").
#' @param model_block Character. Model block.
#' @param gq_block Character. Generated quantities block.
#' @return Character. Complete Stan program.
#' @noRd
assemble_stan_program <- function(fn_block, data_block, par_block,
                                  tp_block, model_block, gq_block) {
  parts <- character()
  if (nzchar(fn_block)) parts <- c(parts, fn_block)
  parts <- c(parts, data_block, as.character(par_block))
  if (nzchar(tp_block)) parts <- c(parts, tp_block)
  parts <- c(parts, as.character(model_block), as.character(gq_block))
  paste(parts, collapse = "\n\n")
}
