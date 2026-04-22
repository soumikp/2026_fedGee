###############################################################################
# Fed-GEE: Federated Generalized Estimating Equations
#
# Small-sample variance corrections (applied AFTER convergence):
#
#   correction = "none" — standard sandwich, no correction
#   correction = "KC"   — Kauermann-Carroll (2001): (I - H)^{-1/2} S
#                         Targets bias in the covariance matrix itself.
#                         Better calibrated coverage in moderate K.
#   correction = "MD"   — Mancl-DeRouen (2001): (I - H)^{-1} S
#                         More aggressive; better Type I error control
#                         when K is very small (K < 20).
#
# The correction is applied in score-space (p x p leverage) rather than
# observation-space (N_i x N_i hat matrix) to preserve privacy: sites
# transmit only p x p summary matrices, never adjusted residuals.
# See Supplement C, Remark 1 for theoretical justification.
#
# Sandwich clustering level:
#   sandwich_level = "site"    — one outer product per hospital (default)
#                                accounts for ALL within-hospital correlation
#   sandwich_level = "patient" — sum of per-patient outer products
#                                assumes patients independent within hospital
#
# Degrees of freedom for t-based inference:
#   site-level sandwich    -> df = K - p   (K = number of sites)
#   patient-level sandwich -> df = N - p   (N = total number of patients)
###############################################################################



###############################################################################
# 0. CORRELATION MATRIX HELPER
###############################################################################
get_Ri <- function(corstr, phi, n) {
  if (corstr == "independence" || n == 1)
    return(diag(1, n))
  
  if (corstr == "exchangeable") {
    Ri <- matrix(as.numeric(phi), n, n)
    diag(Ri) <- 1
    return(Ri)
  }
  
  if (corstr == "ar1") {
    exponent <- abs(outer(1:n, 1:n, "-"))
    return(as.numeric(phi)^exponent)
  }
  
  if (corstr == "unstructured") {
    Ri <- diag(1, n)
    if (length(phi) == (n * (n - 1) / 2)) {
      Ri[lower.tri(Ri)] <- phi
      Ri <- Ri + t(Ri) - diag(1, n)
    } else {
      return(diag(1, n))
    }
    return(Ri)
  }
  
  return(diag(1, n))
}

###############################################################################
# INTERNAL: symmetric matrix square-root inverse via eigen decomposition
#
# Computes (M)^{-power} for a symmetric positive-semi-definite matrix M.
#   power = 0.5  -> M^{-1/2}  (KC correction)
#   power = 1.0  -> M^{-1}    (MD correction)
#
# Eigenvalues are floored at `floor_val` before inversion to guard against
# numerical instability when (I - H) is near-singular (i.e. when a single
# site dominates the network). This is conservative: it shrinks the
# correction toward the identity rather than producing extreme inflations.
###############################################################################
.mat_pow_inv <- function(M, power = 0.5, floor_val = 1e-6) {
  eig <- tryCatch(
    eigen(M, symmetric = TRUE),
    error = function(e) NULL
  )
  if (is.null(eig)) return(diag(nrow(M)))
  vals_inv <- 1 / pmax(eig$values^power, floor_val)
  eig$vectors %*% diag(vals_inv, nrow = length(vals_inv)) %*% t(eig$vectors)
}

###############################################################################
# 1. PREPARATION: Local GLM fits for initial values + local GEE for alpha
#
# Site-level covariate handling:
#   Variables with zero within-site variance at EVERY site (e.g. teaching
#   hospital indicator) are "site-level constants". Local GLMs cannot
#   estimate these — they are dropped from the local formula and initialized
#   at 0. The iterative server update recovers them from between-site
#   variation when B_global = sum(B_i) is inverted. See Section 2.
#
# The (Intercept) is always constant within a site but is NOT site-constant
# in the relevant sense — it is excluded from the flagging logic.
###############################################################################
Prep_FedGEE <- function(data_list,
                        main_formula,
                        family_obj,
                        corstr,
                        id_col) {
  
  N_sites <- length(data_list)
  y_name  <- all.vars(main_formula)[1]
  
  # --- Identify site-level constant columns ---
  full_X_example <- model.matrix(main_formula, data = data_list[[1]])
  all_colnames   <- colnames(full_X_example)
  p_full         <- length(all_colnames)
  
  site_constant_flags <- matrix(FALSE, nrow = N_sites, ncol = p_full)
  colnames(site_constant_flags) <- all_colnames
  
  for (i in seq_len(N_sites)) {
    X_i    <- model.matrix(main_formula, data = data_list[[i]])
    col_var <- apply(X_i, 2, var)
    site_constant_flags[i, ] <- (col_var < 1e-15)
  }
  
  # A column is site-level constant only if it is constant at EVERY site
  site_level_constant <- apply(site_constant_flags, 2, all)
  site_level_constant["(Intercept)"] <- FALSE   # intercept is always "constant" but is not a site covariate
  
  constant_cols <- names(which(site_level_constant))
  varying_cols  <- setdiff(all_colnames, c("(Intercept)", constant_cols))
  
  if (length(constant_cols) > 0) {
    message("FedGEE Prep: Detected site-level covariates (zero within-site variance at all sites):\n  ",
            paste(constant_cols, collapse = ", "),
            "\n  Initialized at 0; estimated via between-site variation during server updates.")
  }
  
  # --- Reduced formula for local fits (excludes site-constant terms) ---
  if (length(varying_cols) > 0) {
    reduced_formula <- as.formula(paste(y_name, "~", paste(varying_cols, collapse = " + ")))
  } else {
    reduced_formula <- as.formula(paste(y_name, "~ 1"))
  }
  
  # --- Local GLM fits for initial beta values ---
  # GLM (not GEE) used here: only rough starting values needed; the
  # iterative protocol corrects them. GLM model-based SEs are always
  # invertible unlike sandwich SEs at small sites.
  glm_local <- map(data_list, function(df) {
    fit <- tryCatch(
      suppressWarnings(glm(reduced_formula, data = df, family = family_obj)),
      error = function(e) NULL
    )
    if (is.null(fit) || !fit$converged || any(is.na(coef(fit)))) return(NULL)
    fit
  })
  
  valid_idx <- !map_lgl(glm_local, is.null)
  glm_valid <- glm_local[valid_idx]
  n_failed  <- sum(!valid_idx)
  
  if (n_failed > 0)
    message(sprintf(
      "FedGEE Prep: %d / %d sites excluded from initialization (aliased coefs, non-convergence, or error). They will still participate in the iterative protocol.",
      n_failed, N_sites))
  
  # --- Meta-analyze local GLM estimates for initial values ---
  if (length(glm_valid) == 0) {
    message("FedGEE Prep: All local GLMs failed. Initializing all coefficients at 0.")
    has_intercept   <- "(Intercept)" %in% all_colnames
    expected_names  <- c(if(has_intercept) "(Intercept)" else NULL, varying_cols)
    beta_reduced    <- matrix(0, nrow = length(expected_names), ncol = 1)
    rownames(beta_reduced) <- expected_names
    se_meta_reduced <- rep(NA_real_, length(expected_names))
  } else {
    beta_hat_list <- map(glm_valid, ~ as.matrix(coef(.x)))
    vcov_list     <- map(glm_valid, ~ as.matrix(vcov(.x)))
    inv_vcov_list <- map(vcov_list, ~ tryCatch(solve(.x), error = function(e) NULL))
    
    ok <- !map_lgl(inv_vcov_list, is.null)
    if (sum(ok) == 0) {
      beta_reduced    <- as.matrix(coef(glm_valid[[1]]))
      se_meta_reduced <- sqrt(diag(vcov(glm_valid[[1]])))
    } else {
      beta_hat_list <- beta_hat_list[ok]
      inv_vcov_list <- inv_vcov_list[ok]
      den          <- Reduce(`+`, inv_vcov_list)
      num          <- Reduce(`+`, map2(inv_vcov_list, beta_hat_list, `%*%`))
      beta_reduced    <- solve(den, num)
      se_meta_reduced <- sqrt(diag(solve(den)))
    }
  }
  
  # --- Local GEE fits to extract working correlation alpha ---
  if (corstr != "independence") {
    alpha_list <- map(data_list, function(df) {
      df[[".id_var"]] <- df[[id_col]]
      fit <- tryCatch(
        geeglm(reduced_formula, data = df, family = family_obj,
               id = .id_var, corstr = corstr),
        error = function(e) NULL
      )
      if (is.null(fit)) return(0)
      a <- fit$geese$alpha
      if (length(a) == 0) 0 else a
    })
  } else {
    alpha_list <- map(data_list, ~ 0)
  }
  
  # --- Assemble full-dimension initial beta (site-constant cols = 0) ---
  reduced_names <- rownames(beta_reduced)
  if (is.null(reduced_names) && length(glm_valid) > 0) {
    reduced_names <- names(coef(glm_valid[[1]]))
  }
  
  initial_values <- matrix(0, nrow = p_full, ncol = 1)
  rownames(initial_values) <- all_colnames
  initial_se <- rep(NA_real_, p_full)
  names(initial_se) <- all_colnames
  
  for (nm in reduced_names) {
    if (nm %in% all_colnames) {
      i_full    <- which(all_colnames == nm)
      i_reduced <- which(reduced_names == nm)
      initial_values[i_full, 1] <- beta_reduced[i_reduced, 1]
      initial_se[i_full]        <- se_meta_reduced[i_reduced]
    }
  }
  
  # Guarantee every site has an alpha entry
  alpha_full <- vector("list", N_sites)
  for (i in seq_len(N_sites))
    alpha_full[[i]] <- if (i <= length(alpha_list)) alpha_list[[i]] else 0
  
  list(
    initial_values  = initial_values,
    initial_se      = initial_se,
    alpha_list      = alpha_full,
    valid_sites     = which(valid_idx),
    constant_cols   = constant_cols,
    reduced_formula = reduced_formula
  )
}

###############################################################################
# 2. SITE-LEVEL COMPUTATION
#
# Computes B_i (bread), S_i (score), and M_i (meat) for one site.
#
# correction argument controls the small-sample adjustment applied to the
# score vector before forming the meat:
#
#   "none"  M_i = S_i S_i'                           (standard)
#   "KC"    S_tilde = (I - H_ii)^{-1/2} S_i          (Kauermann-Carroll)
#   "MD"    S_tilde = (I - H_ii)^{-1}   S_i          (Mancl-DeRouen)
#
# where H_ii = B_i B_global^{-1} is the p x p score-space leverage for
# site i. B_global must be provided when correction != "none".
#
# At sandwich_level = "patient", the same logic applies per-patient using
# per-patient bread B_j and score S_j, with leverage H_jj = B_j B_global^{-1}.
###############################################################################
get_site_stats <- function(site_data,
                           beta_global,
                           alpha,
                           main_formula,
                           family_obj,
                           id_col,
                           corstr,
                           sandwich_level = "site",
                           correction     = "none",
                           B_global       = NULL) {
  
  stopifnot(correction %in% c("none", "KC", "MD"))
  stopifnot(sandwich_level %in% c("site", "patient"))
  
  beta_curr <- as.matrix(beta_global)
  y_name    <- all.vars(main_formula)[1]
  Full_X    <- model.matrix(main_formula, data = site_data)
  Full_y    <- site_data[[y_name]]
  p         <- ncol(Full_X)
  
  cluster_IDs    <- unique(site_data[[id_col]])
  n_clust        <- length(cluster_IDs)
  Bread_clusters <- vector("list", n_clust)
  Score_clusters <- vector("list", n_clust)
  
  idx <- 0L
  for (j in cluster_IDs) {
    idx     <- idx + 1L
    row_idx <- which(site_data[[id_col]] == j)
    X_j     <- Full_X[row_idx, , drop = FALSE]
    y_j     <- Full_y[row_idx]
    n_j     <- length(y_j)
    
    eta_j      <- as.vector(X_j %*% beta_curr)
    mu_j       <- family_obj$linkinv(eta_j)
    dmu_deta_j <- family_obj$mu.eta(eta_j)
    var_mu_j   <- pmax(family_obj$variance(mu_j), 1e-12)
    D_j        <- dmu_deta_j * X_j          # n_j x p
    A_half_j   <- sqrt(var_mu_j)
    
    if (corstr == "independence" || n_j == 1) {
      V_inv_diag <- 1 / var_mu_j
      V_inv_r    <- V_inv_diag * (y_j - mu_j)
      V_inv_D    <- V_inv_diag * D_j
    } else {
      R_j     <- get_Ri(corstr, alpha, n_j)
      V_j     <- (A_half_j %o% A_half_j) * R_j
      V_inv_j <- tryCatch(solve(V_j), error = function(e) NULL)
      if (is.null(V_inv_j)) next
      V_inv_r <- V_inv_j %*% (y_j - mu_j)
      V_inv_D <- V_inv_j %*% D_j
    }
    
    Bread_clusters[[idx]] <- crossprod(D_j, V_inv_D)   # B_j = D' V^{-1} D
    Score_clusters[[idx]] <- crossprod(D_j, V_inv_r)   # S_j = D' V^{-1} r
  }
  
  keep           <- !sapply(Bread_clusters, is.null)
  Bread_clusters <- Bread_clusters[keep]
  Score_clusters <- Score_clusters[keep]
  if (length(Bread_clusters) == 0) return(NULL)
  
  B_site <- Reduce(`+`, Bread_clusters)
  S_site <- Reduce(`+`, Score_clusters)
  
  # ---------------------------------------------------------------------------
  # Helper: apply score-space correction given a bread and score vector.
  # power = 0.5 for KC, power = 1.0 for MD.
  # ---------------------------------------------------------------------------
  .apply_correction <- function(B_unit, S_unit, B_global_inv, power) {
    H    <- B_unit %*% B_global_inv      # leverage: H_ii or H_jj
    I_H  <- diag(p) - H
    A    <- .mat_pow_inv(I_H, power = power)
    A %*% S_unit
  }
  
  # ---------------------------------------------------------------------------
  # Compute meat according to sandwich_level and correction
  # ---------------------------------------------------------------------------
  need_correction <- (correction != "none") && !is.null(B_global)
  power <- switch(correction, "KC" = 0.5, "MD" = 1.0, 0)   # 0 unused for "none"
  
  if (sandwich_level == "site") {
    
    if (need_correction) {
      B_global_inv <- tryCatch(solve(B_global), error = function(e) NULL)
      S_adj <- if (!is.null(B_global_inv)) {
        .apply_correction(B_site, S_site, B_global_inv, power)
      } else S_site
    } else {
      S_adj <- S_site
    }
    M_site <- tcrossprod(S_adj)
    
  } else {  # sandwich_level == "patient"
    
    if (need_correction) {
      B_global_inv <- tryCatch(solve(B_global), error = function(e) NULL)
      M_list <- lapply(seq_along(Score_clusters), function(i) {
        S_adj <- if (!is.null(B_global_inv)) {
          .apply_correction(Bread_clusters[[i]], Score_clusters[[i]], B_global_inv, power)
        } else Score_clusters[[i]]
        tcrossprod(S_adj)
      })
    } else {
      M_list <- lapply(Score_clusters, tcrossprod)
    }
    M_site <- Reduce(`+`, M_list)
  }
  
  list(
    Bread      = B_site,
    Score      = S_site,
    Meat       = M_site,
    n_clusters = sum(keep)
  )
}

###############################################################################
# 3. SERVER-LEVEL AGGREGATION
#
# Aggregates site outputs, performs Newton-Raphson update, computes
# sandwich variance, and derives degrees of freedom.
###############################################################################
aggregate_FedGEE <- function(site_results_list,
                             current_beta,
                             sandwich_level = "site") {
  
  param_names <- names(current_beta) %||% rownames(current_beta)
  valid       <- Filter(Negate(is.null), site_results_list)
  if (length(valid) == 0) return(NULL)
  
  B_total    <- Reduce(`+`, map(valid, "Bread"))
  S_total    <- Reduce(`+`, map(valid, "Score"))
  M_total    <- Reduce(`+`, map(valid, "Meat"))
  n_sites    <- length(valid)
  n_patients <- sum(map_dbl(valid, "n_clusters"))
  
  # Newton-Raphson update
  update_step <- tryCatch(solve(B_total, S_total), error = function(e) NULL)
  if (is.null(update_step)) return(NULL)
  new_beta <- current_beta + update_step
  
  p     <- length(new_beta)
  B_inv <- tryCatch(solve(B_total), error = function(e) NULL)
  if (is.null(B_inv)) return(NULL)
  
  vcov_mat <- B_inv %*% M_total %*% B_inv
  se_beta  <- setNames(sqrt(diag(vcov_mat)), param_names)
  
  # df = (independent clusters) - (parameters)
  if (sandwich_level == "site") {
    df_residual    <- n_sites - p
    total_clusters <- n_sites
  } else {
    df_residual    <- n_patients - p
    total_clusters <- n_patients
  }
  
  list(
    estimate       = new_beta,
    se             = se_beta,
    vcov           = vcov_mat,
    Bread          = B_total,
    Meat           = M_total,
    df_residual    = df_residual,
    total_clusters = total_clusters,
    n_sites        = n_sites,
    n_patients     = n_patients
  )
}

###############################################################################
# 4. TRAINING LOOP
#
# Two-pass design:
#   Pass 1 (iterative): standard sandwich (correction = "none") to find beta
#   Pass 2 (post-convergence): recompute meat with the chosen correction
#         using the final converged B_global
#
# Applying the correction during iteration would require B_global to be
# updated every step, creating a circular dependency. Since the correction
# affects only the variance estimate (not the point estimate), this
# separation is exact for the beta and conservative for the SE.
###############################################################################
train_FedGEE <- function(data_list,
                         initial_beta,
                         main_formula,
                         alpha,
                         family_obj,
                         corstr,
                         id_col,
                         n_iter         = 50,
                         tol            = 1e-8,
                         sandwich_level = "site",
                         correction     = "none",
                         verbose        = TRUE) {
  
  stopifnot(sandwich_level %in% c("site", "patient"))
  stopifnot(correction %in% c("none", "KC", "MD"))
  
  if (verbose) {
    cat("=== Fed-GEE ===\n")
    cat("  Family          :", family_obj$family, "|", family_obj$link, "\n")
    cat("  Working corr    :", corstr, "\n")
    cat("  Sandwich level  :", sandwich_level, "\n")
    cat("  SS correction   :", correction, "\n")
    cat("  Sites           :", length(data_list), "\n\n")
  }
  
  N_sites   <- length(data_list)
  beta_curr <- initial_beta
  converged <- FALSE
  history   <- list(beta_curr)
  
  for (k in seq_len(n_iter)) {
    
    # --- Pass 1 site computation: no correction during iteration ---
    site_outputs <- vector("list", N_sites)
    for (i in seq_len(N_sites)) {
      if (is.null(alpha[[i]])) next
      site_outputs[[i]] <- get_site_stats(
        site_data      = data_list[[i]],
        beta_global    = beta_curr,
        alpha          = alpha[[i]],
        main_formula   = main_formula,
        family_obj     = family_obj,
        id_col         = id_col,
        corstr         = corstr,
        sandwich_level = sandwich_level,
        correction     = "none",
        B_global       = NULL
      )
    }
    
    global_update <- aggregate_FedGEE(site_outputs, beta_curr, sandwich_level)
    if (is.null(global_update)) {
      if (verbose) cat("  Singularity in aggregation — stopping.\n")
      return(NULL)
    }
    
    beta_next <- global_update$estimate
    diff_norm <- sqrt(sum((beta_curr - beta_next)^2))
    if (verbose) cat(sprintf("  Iter %3d  | update norm = %.3e\n", k, diff_norm))
    
    beta_curr        <- beta_next
    history[[k + 1]] <- beta_curr
    
    if (diff_norm < tol) {
      converged <- TRUE
      if (verbose) cat("  Converged.\n\n")
      break
    }
    if (diff_norm > 1e4) {
      if (verbose) cat("  Divergence detected — stopping.\n")
      return(NULL)
    }
  }
  
  # --- Pass 2: recompute meat at converged beta with chosen correction ---
  B_global_final <- global_update$Bread
  
  site_outputs_final <- vector("list", N_sites)
  for (i in seq_len(N_sites)) {
    if (is.null(alpha[[i]])) next
    site_outputs_final[[i]] <- get_site_stats(
      site_data      = data_list[[i]],
      beta_global    = beta_curr,
      alpha          = alpha[[i]],
      main_formula   = main_formula,
      family_obj     = family_obj,
      id_col         = id_col,
      corstr         = corstr,
      sandwich_level = sandwich_level,
      correction     = correction,
      B_global       = B_global_final
    )
  }
  
  final_result <- aggregate_FedGEE(site_outputs_final, beta_curr, sandwich_level)
  if (is.null(final_result)) {
    warning("FedGEE: final variance pass failed; returning uncorrected sandwich.")
    final_result <- global_update
  }
  
  coef_names <- rownames(initial_beta) %||% names(initial_beta)
  
  structure(
    list(
      coefficients   = setNames(as.vector(beta_curr), coef_names),
      se             = final_result$se,
      vcov           = final_result$vcov,
      Bread          = final_result$Bread,
      Meat           = final_result$Meat,
      df_residual    = final_result$df_residual,
      total_clusters = final_result$total_clusters,
      n_sites        = final_result$n_sites,
      n_patients     = final_result$n_patients,
      iterations     = k,
      converged      = converged,
      history        = history,
      sandwich_level = sandwich_level,
      correction     = correction,
      corstr         = corstr
    ),
    class = "FedGEE"
  )
}

###############################################################################
# 5. PRINT METHOD
###############################################################################
#' Print Method for Federated GEE
#'
#' @param x An object of class \code{FedGEE}
#' @param ... Additional arguments
#' @export
print.FedGEE <- function(x, ...) {
  cat("Federated GEE\n")
  cat("  Working correlation :", x$corstr, "\n")
  cat("  Sandwich level      :", x$sandwich_level, "\n")
  cat("  SS correction       :", x$correction, "\n")
  cat("  Sites               :", x$n_sites,
      "| Patients:", x$n_patients, "\n")
  cat("  df residual         :", x$df_residual, "\n")
  cat("  Converged           :", x$converged,
      "in", x$iterations, "iterations\n")
  
  est <- x$coefficients
  se  <- x$se
  z   <- est / se
  
  if (!is.null(x$df_residual) && x$df_residual > 0) {
    pval    <- 2 * pt(-abs(z), df = x$df_residual)
    ci_mult <- qt(0.975, df = x$df_residual)
    cat(sprintf("  Inference           : t(%d)\n\n", x$df_residual))
  } else {
    pval    <- 2 * pnorm(-abs(z))
    ci_mult <- qnorm(0.975)
    cat("  Inference           : z\n\n")
  }
  
  tab <- data.frame(
    Estimate = round(est, 4),
    SE       = round(se, 4),
    CI.lo    = round(est - ci_mult * se, 4),
    CI.hi    = round(est + ci_mult * se, 4),
    `z/t`    = round(z, 3),
    p.value  = round(pval, 4),
    check.names = FALSE
  )
  print(tab)
  invisible(x)
}

###############################################################################
# 6. CONVENIENCE WRAPPER
#
# correction = "none" | "KC" | "MD"
###############################################################################
#' Federated Generalized Estimating Equations (FedGEE)
#'
#' Fits a federated generalized estimating equation model with small-sample variance corrections.
#'
#' @param data_list A list of data frames, one per site.
#' @param main_formula The model formula.
#' @param family_obj A family object (e.g., \code{binomial()}).
#' @param corstr Working correlation structure (e.g., "independence", "exchangeable").
#' @param id_col Character string naming the cluster/patient ID column.
#' @param sandwich_level Level to compute sandwich variance ("site" or "patient").
#' @param correction Small-sample variance correction ("none", "KC", or "MD").
#' @param n_iter Maximum number of iterations.
#' @param tol Convergence tolerance.
#' @param verbose Logical; if \code{TRUE}, prints fit progress.
#' @return An object of class \code{FedGEE}.
#' @import purrr
#' @import dplyr
#' @import geepack
#' @import Matrix
#' @import tibble
#' @importFrom stats as.formula coef gaussian glm pnorm pt qnorm qt var vcov
#' @importFrom methods as
#' @export
fedgee <- function(data_list,
                   main_formula,
                   family_obj     = binomial(link = "logit"),
                   corstr         = "independence",
                   id_col         = "pat_id",
                   sandwich_level = "site",
                   correction     = "KC",
                   n_iter         = 50,
                   tol            = 1e-8,
                   verbose        = TRUE) {
  
  prep <- Prep_FedGEE(data_list, main_formula, family_obj, corstr, id_col)
  
  train_FedGEE(
    data_list      = data_list,
    initial_beta   = prep$initial_values,
    main_formula   = main_formula,
    alpha          = prep$alpha_list,
    family_obj     = family_obj,
    corstr         = corstr,
    id_col         = id_col,
    n_iter         = n_iter,
    tol            = tol,
    sandwich_level = sandwich_level,
    correction     = correction,
    verbose        = verbose
  )
}