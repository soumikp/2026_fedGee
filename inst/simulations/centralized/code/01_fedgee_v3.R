###############################################################################
# Fed-GEE simulation code -- v3 (2026-09-27)
#
# Changes from v2 (sim_calibration):
#   1. "fg" is now TRUE Fay-Graubard at the SITE level (the federated cluster),
#      F_i = diag{(1 - min(0.75, diag(B_i B^{-1})))^{-1/2}} on the raw site score.
#      The v2 "fg" (observation-level HC2) is kept as "hc2obs".
#   2. Bell-McCaffrey df in score space (df_bm), design-only, from {B_i}.
#   3. FG's own df (dtildeH) via saws::saws(method = "d5") for the fg variant.
#   4. t(K - 1) reference added; the v2 Satterthwaite-labelled df is renamed
#      df_nu (realized-score effective-contributor count, kept for comparison).
#   5. site_z/site_t/site_satt collapse into one "site" variant (uncorrected
#      sandwich); the reference distribution is carried by `inference`.
###############################################################################
###############################################################################
# Federated GEE — v2 implementation
#
# Key change from v1: ONE call to fedgee() per replicate produces ALL variance
# corrections (uncorrected, KC, MD, FG, VD) and ALL inference variants (z, t-Kp,
# t-Satt). Point estimates are identical across these — only vcov / df differ.
#
# This restructure:
#   1. Cuts simulation runtime ~6x (was: re-fit per correction)
#   2. Centralizes correction logic in compute_all_variances()
#   3. Records per-site diagnostics needed for stress analysis
#   4. Adds Fay-Graubard (FG) and Variance Decomposition (VD) corrections
#
# Architecture:
#   Prep_FedGEE()        — local GLM init (UNCHANGED, returns site_betas now)
#   get_site_stats_v2()  — single cluster loop returns all primitives
#   train_FedGEE_v2()    — iterate to convergence; final pass collects all stats
#   compute_all_variances() — given primitives, return all variance variants
#   fedgee_v2()          — convenience wrapper
###############################################################################

library(geepack)
library(dplyr)
library(purrr)
library(tibble)
library(Matrix)

###############################################################################
# 0. CORRELATION MATRIX HELPER (unchanged)
###############################################################################
get_Ri <- function(corstr, phi, n) {
  if (corstr == "independence" || n == 1) return(diag(1, n))
  if (corstr == "exchangeable") {
    Ri <- matrix(as.numeric(phi), n, n); diag(Ri) <- 1; return(Ri)
  }
  if (corstr == "ar1") {
    exponent <- abs(matrix(1:n, nrow = n, ncol = n, byrow = TRUE) - (1:n))
    return(as.numeric(phi)^exponent)
  }
  if (corstr == "unstructured") {
    Ri <- diag(1, n)
    if (length(phi) == (n * (n - 1) / 2)) {
      Ri[lower.tri(Ri)] <- phi; Ri <- Ri + t(Ri) - diag(1, n)
    } else { return(diag(1, n)) }
    return(Ri)
  }
  return(diag(1, n))
}

###############################################################################
# Helper: matrix power via eigendecomposition. M MUST BE SYMMETRIC.
# Reports whether the eigenvalue floor was activated (near-singular I - G).
###############################################################################
.mat_pow_safe <- function(M, power = -0.5, tol = 1e-6) {
  M <- (M + t(M)) / 2                      # enforce exact symmetry
  eig <- tryCatch(eigen(M, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig)) return(list(mat = NULL, clipped = TRUE, min_eig = NA_real_))

  vals <- eig$values
  min_eig <- min(vals)
  clipped <- any(vals < tol)
  vals_safe <- pmax(vals, tol)

  mat <- eig$vectors %*% diag(vals_safe^power) %*% t(eig$vectors)
  list(mat = mat, clipped = clipped, min_eig = min_eig)
}

###############################################################################
# Helper: score-space leverage correction via the SIMILARITY TRANSFORM.
#
#   A_i = B^{1/2} (I - G_i)^{-power} B^{-1/2},   G_i = B^{-1/2} B_i B^{-1/2}
#
# The naive leverage H_i = B_i B^{-1} is a product of two symmetric matrices
# and is NOT symmetric, so eigen(I - H_i, symmetric = TRUE) silently discards
# the upper triangle and returns wrong eigenvalues. G_i is symmetric PSD with
# the SAME eigenvalues as H_i and satisfies sum_i G_i = I -- the identity that
# makes the correction exact under the working model M = phi * B.
#
# Returns the operator plus the clipping diagnostics the caller records.
###############################################################################
.leverage_op <- function(B_site, B_half, B_neghalf, power, tol = 0.05) {
  p   <- nrow(B_site)
  G_i <- B_neghalf %*% B_site %*% B_neghalf
  res <- .mat_pow_safe(diag(p) - G_i, power = -power, tol = tol)
  if (is.null(res$mat)) {
    return(list(mat = NULL, clipped = res$clipped, min_eig = res$min_eig))
  }
  list(mat     = B_half %*% res$mat %*% B_neghalf,
       clipped = res$clipped,
       min_eig = res$min_eig)
}

###############################################################################
# 1. PREPARATION (modified to also return per-site beta vectors for VD)
###############################################################################
Prep_FedGEE_v2 <- function(data_list, main_formula, family_obj, corstr, id_col) {

  N_sites <- length(data_list)
  y_name  <- all.vars(main_formula)[1]

  # 1. Identify site-constant columns
  full_X_example <- model.matrix(main_formula, data = data_list[[1]])
  all_colnames   <- colnames(full_X_example)
  p_full         <- length(all_colnames)

  site_constant_flags <- matrix(FALSE, nrow = N_sites, ncol = p_full)
  colnames(site_constant_flags) <- all_colnames

  for (i in seq_len(N_sites)) {
    X_i <- model.matrix(main_formula, data = data_list[[i]])
    site_constant_flags[i, ] <- (apply(X_i, 2, var) < 1e-15)
  }

  site_level_constant <- apply(site_constant_flags, 2, all)
  site_level_constant["(Intercept)"] <- FALSE
  constant_cols <- names(which(site_level_constant))
  varying_cols  <- setdiff(all_colnames, c("(Intercept)", constant_cols))

  # Map design-matrix columns back to source TERMS. model.matrix() expands a
  # factor `f` into columns f2, f3, ... which are NOT variables in the data, so
  # a formula built from column names fails with "object 'f2' not found" and
  # every local fit silently returns NULL. attr(,"assign") gives the generating
  # term for each column (0 = intercept); drop a term only if ALL its columns
  # are site-level constant.
  mm_assign   <- attr(full_X_example, "assign")
  term_labels <- attr(terms(main_formula, data = data_list[[1]]), "term.labels")

  varying_terms <- term_labels
  if (!is.null(mm_assign) && length(term_labels) > 0) {
    term_is_constant <- vapply(seq_along(term_labels), function(k) {
      cols_k <- which(mm_assign == k)
      length(cols_k) > 0 && all(site_level_constant[cols_k])
    }, logical(1))
    varying_terms <- term_labels[!term_is_constant]
  }

  if (length(varying_terms) > 0) {
    reduced_rhs     <- paste(varying_terms, collapse = " + ")
    reduced_formula <- as.formula(paste(y_name, "~", reduced_rhs))
  } else {
    reduced_formula <- as.formula(paste(y_name, "~ 1"))
  }

  # 2. Local GLM fits — capture per-site beta vectors AND meta-analytic init
  glm_local <- data_list %>%
    map(function(df) {
      fit <- tryCatch(
        suppressWarnings(glm(reduced_formula, data = df, family = family_obj)),
        error = function(e) NULL)
      if (is.null(fit) || !fit$converged || any(is.na(coef(fit)))) return(NULL)
      fit
    })

  valid_idx <- !map_lgl(glm_local, is.null)
  glm_valid <- glm_local[valid_idx]

  # Per-site beta on the REDUCED parameter space (for VD computation)
  site_betas_reduced <- map(glm_local, function(fit) {
    if (is.null(fit)) return(NULL)
    coef(fit)
  })

  if (length(glm_valid) == 0) {
    initial_values <- matrix(0, nrow = p_full, ncol = 1)
    rownames(initial_values) <- all_colnames
    return(list(
      initial_values = initial_values,
      alpha_list = as.list(rep(0, N_sites)),
      valid_sites = integer(0),
      constant_cols = constant_cols,
      reduced_formula = reduced_formula,
      site_betas_reduced = site_betas_reduced,
      reduced_colnames = NULL
    ))
  }

  # 3. Inverse-variance meta-analysis for initial values
  beta_hat_reduced <- map(glm_valid, ~ as.matrix(coef(.x)))
  var_beta_reduced <- map(glm_valid, ~ as.matrix(vcov(.x)))
  inv_var <- map(var_beta_reduced, ~ tryCatch(solve(.x), error = function(e) NULL))
  ok <- !map_lgl(inv_var, is.null)
  beta_hat_reduced <- beta_hat_reduced[ok]; inv_var <- inv_var[ok]

  if (length(inv_var) == 0) {
    beta_reduced <- as.matrix(coef(glm_valid[[1]]))
  } else {
    inv_var_beta <- map2(inv_var, beta_hat_reduced, ~ .x %*% .y)
    den_FE <- Reduce(`+`, inv_var); num_FE <- Reduce(`+`, inv_var_beta)
    beta_reduced <- solve(den_FE, num_FE)
  }

  # 4. Alpha extraction (only if non-independence)
  if (corstr != "independence") {
    gee_local <- data_list %>%
      map(function(df) {
        df[[".id_var"]] <- df[[id_col]]
        tryCatch(
          geeglm(formula = reduced_formula, data = df, family = family_obj,
                 id = .id_var, corstr = corstr),
          error = function(e) NULL)
      })
    alpha <- map(gee_local, function(m) {
      if (is.null(m)) return(0)
      a <- m$geese$alpha
      if (length(a) == 0) return(0); return(a)
    })
  } else {
    alpha <- map(data_list, ~ 0)
  }

  # 5. Assemble full initial beta — ZERO INIT
  # The Newton-Raphson loop is fast and converges in ~6 iterations from zero.
  # Meta-analytic init from local GLMs becomes unstable at high p (sparse local
  # coverage of binary covariates produces wild coefficient estimates that
  # blow up N-R on iteration 1). Zero init is robust across all regimes.
  reduced_names <- rownames(beta_reduced)
  if (is.null(reduced_names)) reduced_names <- names(coef(glm_valid[[1]]))
  
  initial_values <- matrix(0, nrow = p_full, ncol = 1)
  rownames(initial_values) <- all_colnames

  alpha_full <- vector("list", N_sites)
  for (i in seq_len(N_sites)) {
    alpha_full[[i]] <- if (i <= length(alpha)) alpha[[i]] else 0
  }

  list(
    initial_values     = initial_values,
    alpha_list         = alpha_full,
    valid_sites        = which(valid_idx),
    constant_cols      = constant_cols,
    reduced_formula    = reduced_formula,
    site_betas_reduced = site_betas_reduced,   # NEW: per-site betas on reduced space
    reduced_colnames   = reduced_names         # NEW: names for VD reconstruction
  )
}

###############################################################################
# 2. SITE-LEVEL COMPUTATION v2
#
# Returns ALL primitives needed for any correction. One pass through clusters.
# Caller decides which corrections to apply at the variance step.
#
# Returned object:
#   B_site         -- p x p site bread (sum over patients)
#   S_site         -- p x 1 site score  (sum over patients)
#   B_clusters     -- list of per-patient bread matrices
#   S_clusters     -- list of per-patient score vectors
#   H_diag_obs     -- observation-level FG leverage diagonals (computed only at
#                     final pass when B_global supplied)
#   n_clusters     -- number of patients with valid V_inv
#   n_obs          -- total observations at site
###############################################################################
get_site_stats_v2 <- function(site_data,
                              beta_global,
                              alpha,
                              main_formula,
                              family_obj,
                              id_col,
                              corstr,
                              B_global = NULL,        # only at final pass
                              compute_H_obs = FALSE) {  # only at final pass

  beta_curr <- as.matrix(beta_global)
  y_name    <- all.vars(main_formula)[1]
  Full_X    <- model.matrix(main_formula, data = site_data)
  Full_y    <- site_data[[y_name]]
  p         <- ncol(Full_X)
  cluster_IDs <- unique(site_data[[id_col]])

  Bread_clusters <- vector("list", length(cluster_IDs))
  Score_clusters <- vector("list", length(cluster_IDs))
  Score_clusters_fg <- vector("list", length(cluster_IDs))  # FG-adjusted
  H_diag_clusters   <- vector("list", length(cluster_IDs))  # for diagnostics

  # Pre-compute B_global^{-1} once (if FG/H_obs needed)
  B_global_inv <- NULL
  if (compute_H_obs && !is.null(B_global)) {
    B_global_inv <- tryCatch(solve(B_global), error = function(e) NULL)
  }

  total_obs <- 0L
  fg_any_clipped <- FALSE
  max_h_obs_local <- 0

  idx <- 0L
  for (j in cluster_IDs) {
    idx <- idx + 1L

    row_idx <- which(site_data[[id_col]] == j)
    X_j <- Full_X[row_idx, , drop = FALSE]
    y_j <- Full_y[row_idx]
    n_j <- length(y_j)

    eta_j      <- as.vector(X_j %*% beta_curr)
    mu_j       <- family_obj$linkinv(eta_j)
    dmu_deta_j <- family_obj$mu.eta(eta_j)
    var_mu_j   <- pmax(family_obj$variance(mu_j), 1e-12)

    D_j <- dmu_deta_j * X_j  # n_j x p
    A_half_j <- sqrt(var_mu_j)
    r_j <- y_j - mu_j

    if (corstr == "independence" || n_j == 1) {
      V_inv_diag <- 1 / var_mu_j
      V_inv_r <- V_inv_diag * r_j
      V_inv_D <- V_inv_diag * D_j
    } else {
      R_j <- get_Ri(corstr, alpha, n_j)
      V_j <- (A_half_j %o% A_half_j) * R_j
      V_inv_j <- tryCatch(solve(V_j), error = function(e) NULL)
      if (is.null(V_inv_j)) next
      V_inv_r <- V_inv_j %*% r_j
      V_inv_D <- V_inv_j %*% D_j
    }

    Bread_clusters[[idx]] <- crossprod(D_j, V_inv_D)
    Score_clusters[[idx]] <- crossprod(D_j, V_inv_r)
    total_obs <- total_obs + n_j

    # ---- FG: observation-level leverage and adjusted score ----
    # H_j = D_j B_global^{-1} D_j^T V_j^{-1}    (n_j x n_j)
    # diag(A B) = rowSums(A * t(B)) where A = D_j B_inv (n_j x p) and B = V_inv D_j (n_j x p)
    # Note: V_inv_D is ALREADY V_j^{-1} D_j (n_j x p), so we want
    #       diag(H_j) = diag( (D_j B_inv) (V_inv D_j)^T ... )
    # Actually H_j = D_j B_inv D_j^T V_j^{-1}, and
    # diag(H_j)_k = (D_j B_inv)_k . (V_inv D_j)_k? No — that's not right.
    # H_j has entries H_j[k,l] = (D_j B_inv)_k . D_j[l,] for V=I diagonal,
    # and more generally:
    #   H_j[k,l] = sum_a (D_j B_inv)[k,a] * (D_j^T V_j^{-1})[a,l]
    #            = sum_a DB_inv[k,a] * (V_inv D_j)[l,a]
    # Therefore diag(H_j)[k] = sum_a DB_inv[k,a] * (V_inv D_j)[k,a]
    #                       = rowSums(DB_inv * V_inv_D) <-- element-wise, NOT t()
    if (compute_H_obs && !is.null(B_global_inv)) {
      DB_inv   <- D_j %*% B_global_inv               # n_j x p
      h_diag_j <- rowSums(DB_inv * V_inv_D)          # length n_j  (note: V_inv_D, not t(V_inv_D))
      H_diag_clusters[[idx]] <- h_diag_j

      # FG-adjusted residual: r̃_kk = r_kk / sqrt(1 - h_kk)
      one_minus_h <- 1 - h_diag_j
      if (any(one_minus_h < 1e-6)) {
        fg_any_clipped <- TRUE
        one_minus_h <- pmax(one_minus_h, 1e-6)
      }
      max_h_obs_local <- max(max_h_obs_local, max(h_diag_j, na.rm = TRUE))

      # Recompute V_inv_r with FG-scaled residual to form FG cluster score
      r_j_fg <- r_j / sqrt(one_minus_h)
      if (corstr == "independence" || n_j == 1) {
        V_inv_r_fg <- (1 / var_mu_j) * r_j_fg
      } else {
        V_inv_r_fg <- V_inv_j %*% r_j_fg
      }
      Score_clusters_fg[[idx]] <- crossprod(D_j, V_inv_r_fg)
    }
  }

  keep <- !sapply(Bread_clusters, is.null)
  Bread_clusters <- Bread_clusters[keep]
  Score_clusters <- Score_clusters[keep]
  Score_clusters_fg <- Score_clusters_fg[keep]
  H_diag_clusters <- H_diag_clusters[keep]

  if (length(Bread_clusters) == 0) return(NULL)

  # Safe construction of S_site_fg: only if all clusters have FG scores
  S_site_fg <- NULL
  if (compute_H_obs && all(!sapply(Score_clusters_fg, is.null))) {
    S_site_fg <- Reduce(`+`, Score_clusters_fg)
  }

  list(
    B_site            = Reduce(`+`, Bread_clusters),
    S_site            = Reduce(`+`, Score_clusters),
    S_site_fg         = S_site_fg,
    B_clusters        = Bread_clusters,
    S_clusters        = Score_clusters,
    S_clusters_fg     = Score_clusters_fg,
    H_diag_clusters   = H_diag_clusters,
    n_clusters        = sum(keep),
    n_obs             = total_obs,
    fg_any_clipped    = fg_any_clipped,
    max_h_obs_site    = max_h_obs_local
  )
}

###############################################################################
# 3. SERVER UPDATE (used during iteration only)
#
# Light-weight aggregation for parameter updates. Variance/correction logic
# lives in compute_all_variances() instead.
###############################################################################
server_update <- function(site_results, current_beta) {
  valid <- Filter(Negate(is.null), site_results)
  if (length(valid) == 0) return(NULL)

  B_total <- Reduce(`+`, map(valid, "B_site"))
  S_total <- Reduce(`+`, map(valid, "S_site"))

  upd <- tryCatch(solve(B_total, S_total), error = function(e) NULL)
  if (is.null(upd)) return(NULL)

  list(
    new_beta = current_beta + upd,
    B_total  = B_total,
    S_total  = S_total
  )
}

###############################################################################
# 4. COMPUTE ALL VARIANCE/INFERENCE VARIANTS
#
# Given the final-pass site_results and final B_global, compute every variance
# correction we care about. Returns a list of variance variants, each with
# vcov, se, df_kp, df_satt, plus diagnostics.
#
# Variants returned:
#   patient_z       -- patient sandwich, z-based inference (df = Inf)
#   site_z          -- site sandwich, z (df = Inf)
#   site_t          -- site sandwich, t with df = K - p
#   site_satt       -- site sandwich, t with Bell-McCaffrey Satterthwaite df
#   kc              -- KC correction (I-H)^{-1/2} on site scores, t(satt)
#   md              -- MD correction (I-H)^{-1} on site scores, t(satt)
#   fg              -- FG correction with obs-level leverage, t(satt)
#   vd              -- variance decomposition: within(patient) + between(sites)
#
# Diagnostics:
#   min_eig_IH_kc   -- smallest eigenvalue of (I - H_ii) across sites (KC stress)
#   min_eig_IH_md   -- same for MD
#   max_obs_lev     -- max observation-level leverage diag (FG stress)
#   kc_clipped      -- # sites where KC eigenvalue clip activated
#   md_clipped      -- # sites where MD eigenvalue clip activated
#   n_sites_vd      -- # sites contributing to between-site variance in VD
###############################################################################
compute_all_variances <- function(site_results,
                                  beta_hat,
                                  B_global,
                                  site_betas_reduced = NULL,
                                  reduced_colnames = NULL,
                                  full_colnames = NULL) {

  beta_hat <- as.numeric(beta_hat)  # ensure plain vector for indexing
  valid <- Filter(Negate(is.null), site_results)
  K <- length(valid); p <- length(beta_hat)

  B_inv <- tryCatch(solve(B_global), error = function(e) NULL)
  if (is.null(B_inv)) return(NULL)

  # Symmetric square roots of B_global, formed once and reused by KC/MD.
  # Needed for the similarity transform G_i = B^{-1/2} B_i B^{-1/2}.
  B_half    <- .mat_pow_safe(B_global, power =  0.5, tol = 1e-12)$mat
  B_neghalf <- .mat_pow_safe(B_global, power = -0.5, tol = 1e-12)$mat
  if (is.null(B_half) || is.null(B_neghalf)) return(NULL)

  # ============================================================================
  # PATIENT-LEVEL MEAT (sum of per-patient outer products)
  # ============================================================================
  M_patient <- Reduce(`+`, lapply(valid, function(s) {
    Reduce(`+`, lapply(s$S_clusters, tcrossprod))
  }))
  vcov_patient <- B_inv %*% M_patient %*% B_inv
  n_patients <- sum(map_dbl(valid, "n_clusters"))

  # ============================================================================
  # SITE-LEVEL MEAT (uncorrected) — building block for KC/MD/FG
  # ============================================================================
  S_sites <- map(valid, "S_site")        # list of p x 1 score vectors
  M_site_uncorrected <- Reduce(`+`, lapply(S_sites, tcrossprod))
  vcov_site <- B_inv %*% M_site_uncorrected %*% B_inv

  # ============================================================================
  # KC CORRECTION: (I - H_ii)^{-1/2} on site scores
  # ============================================================================
  kc_min_eigs <- numeric(K); kc_clipped_count <- 0L
  S_kc_list <- vector("list", K)
  T_kc_list <- vector("list", K)   # v3: raw-score operator T_i, reused by the BM df
  for (i in seq_len(K)) {
    res <- .leverage_op(valid[[i]]$B_site, B_half, B_neghalf, power = 0.5, tol = 0.05)
    kc_min_eigs[i] <- res$min_eig
    if (res$clipped) kc_clipped_count <- kc_clipped_count + 1L
    T_kc_list[[i]] <- if (is.null(res$mat)) diag(p) else res$mat
    S_kc_list[[i]] <- T_kc_list[[i]] %*% S_sites[[i]]
  }
  M_kc <- Reduce(`+`, lapply(S_kc_list, tcrossprod))
  vcov_kc <- B_inv %*% M_kc %*% B_inv

  # ============================================================================
  # MD CORRECTION: (I - H_ii)^{-1} on site scores
  # ============================================================================
  md_min_eigs <- numeric(K); md_clipped_count <- 0L
  S_md_list <- vector("list", K)
  T_md_list <- vector("list", K)
  for (i in seq_len(K)) {
    res <- .leverage_op(valid[[i]]$B_site, B_half, B_neghalf, power = 1, tol = 0.05)
    md_min_eigs[i] <- res$min_eig
    if (res$clipped) md_clipped_count <- md_clipped_count + 1L
    T_md_list[[i]] <- if (is.null(res$mat)) diag(p) else res$mat
    S_md_list[[i]] <- T_md_list[[i]] %*% S_sites[[i]]
  }
  M_md <- Reduce(`+`, lapply(S_md_list, tcrossprod))
  vcov_md <- B_inv %*% M_md %*% B_inv

  # ============================================================================
  # FG CORRECTION (v3): Fay & Graubard (2001, Biometrics 57:1198), SITE level.
  #
  # The site is the independent cluster, so FG acts on the raw site score:
  #   S_fg_i = F_i S_i,   F_i = diag{ (1 - min(b, diag(Q_i)))^{-1/2} },
  #   Q_i = B_i B^{-1},   b = 0.75 (FG's default bound).
  # Needs only {S_i, B_i}, so it is computable under federation. Verified to
  # machine precision against saws::saws(method = "d5") and geesmv::GEE.var.fg.
  # ============================================================================
  FG_BOUND <- 0.75
  T_fg_list <- vector("list", K); fg_capped_count <- 0L
  for (i in seq_len(K)) {
    q_diag <- diag(valid[[i]]$B_site %*% B_inv)
    if (any(q_diag > FG_BOUND)) fg_capped_count <- fg_capped_count + 1L
    T_fg_list[[i]] <- diag((1 - pmin(FG_BOUND, q_diag))^(-0.5), nrow = p)
  }
  S_fg_sites <- lapply(seq_len(K), function(i) T_fg_list[[i]] %*% S_sites[[i]])
  M_fg <- Reduce(`+`, lapply(S_fg_sites, tcrossprod))
  vcov_fg <- B_inv %*% M_fg %*% B_inv

  # ============================================================================
  # HC2-OBS (v3 rename): the estimator the v2 code labelled "fg". It scales each
  # OBSERVATION residual by 1/sqrt(1 - h_kk). It is NOT Fay-Graubard; kept only
  # for continuity with the v2/calibration runs.
  # ============================================================================
  S_hc2obs_sites <- map(valid, "S_site_fg")
  if (any(map_lgl(S_hc2obs_sites, is.null))) {
    S_hc2obs_sites <- map(valid, "S_site")
    hc2obs_clipped_count <- NA_integer_
    max_h_obs <- NA_real_
  } else {
    hc2obs_clipped_count <- sum(map_lgl(valid, "fg_any_clipped"))
    max_h_obs <- max(map_dbl(valid, "max_h_obs_site"), na.rm = TRUE)
  }
  M_hc2obs <- Reduce(`+`, lapply(S_hc2obs_sites, tcrossprod))
  vcov_hc2obs <- B_inv %*% M_hc2obs %*% B_inv

  # ============================================================================
  # VARIANCE DECOMPOSITION: V_within (patient sandwich) + V_between (cross-site β)
  #
  # V_between = (1/(K_eff-1)) * sum over sites with valid local fit of
  #             (β̂_i - β̄)(β̂_i - β̄)^T,
  # then divided by K_eff to get variance of mean.
  # On the score scale: V_between contributes additively to vcov_patient.
  # ============================================================================
  vcov_vd <- vcov_patient  # default fallback
  n_sites_vd <- NA_integer_
  if (!is.null(site_betas_reduced) && !is.null(reduced_colnames) && !is.null(full_colnames)) {
    # Reconstruct full-dim beta_i from reduced beta_i (zero-pad site-constants)
    site_beta_mat <- matrix(NA_real_, nrow = length(site_betas_reduced), ncol = p)
    colnames(site_beta_mat) <- full_colnames
    for (i in seq_along(site_betas_reduced)) {
      bi <- site_betas_reduced[[i]]
      if (is.null(bi)) next
      for (nm in names(bi)) {
        if (nm %in% full_colnames) site_beta_mat[i, nm] <- bi[nm]
      }
    }
    # Use only sites with full reduced-space fits available
    has_full <- complete.cases(site_beta_mat[, intersect(reduced_colnames, full_colnames), drop = FALSE])
    valid_betas <- site_beta_mat[has_full, , drop = FALSE]
    n_sites_vd <- nrow(valid_betas)

    if (n_sites_vd >= 3) {
      # Trim site-level beta estimates that are wildly out of range (logistic
      # separation in tiny local fits produces |beta| > 100 sometimes). Use
      # winsorization at +/- 10 (logit scale -- |beta| > 10 means OR > 22000,
      # not a real effect, just a fitting artifact).
      valid_betas[is.na(valid_betas)] <- 0
      valid_betas[abs(valid_betas) > 10] <- NA
      # After trimming, replace remaining NAs with column means or beta_hat
      for (j in seq_len(p)) {
        col_j <- valid_betas[, j]
        if (all(is.na(col_j))) {
          valid_betas[, j] <- beta_hat[j]
        } else if (any(is.na(col_j))) {
          valid_betas[is.na(col_j), j] <- mean(col_j, na.rm = TRUE)
        }
      }
      beta_bar <- colMeans(valid_betas)
      dev <- sweep(valid_betas, 2, beta_bar, "-")
      V_between <- crossprod(dev) / (n_sites_vd - 1) / n_sites_vd
      vcov_vd <- vcov_patient + V_between
    }
  }

  # ============================================================================
  # DEGREES OF FREEDOM (v3). Five reference distributions per variant:
  #   z      : df = Inf
  #   t_k1   : df = K - 1   ("clusters minus cluster-level parameters")
  #   t_kp   : df = K - p   (Mancl-DeRouen / Li-Redden count)
  #   t_nu   : v2 "effective-contributor" count (sum q)^2 / sum q^2 with REALIZED
  #            q_ij. Kept only for comparison: it is a noisy statistic (~K/3
  #            under balance), not a Bell-McCaffrey df.
  #   t_bm   : Bell-McCaffrey df in score space (see bm_df below)
  #   t_fgdf : FG's own df (dtildeH, Fay & Graubard 2001) via saws, FG only
  # ============================================================================
  nu_df <- function(S_list) {
    sapply(seq_len(p), function(j) {
      contribs <- as.numeric(sapply(S_list, function(s) (B_inv[j, ] %*% s)^2))
      if (sum(contribs^2) < 1e-20) return(K - p)
      (sum(contribs))^2 / sum(contribs^2)
    })
  }

  # ----------------------------------------------------------------------------
  # Bell-McCaffrey df in score space.
  #
  # Working model = the Theorem-3 affine model: S_i(b) = U_i - B_i (b - b0),
  # U_i independent with Var(U_i) = tau B_i. Each variant's variance for
  # coefficient j is sum_i (a_i' S_i)^2 with a_i = T_i' B^{-1} e_j, where T_i is
  # the variant's raw-score operator (I, KC, MD, or FG). Under the working model
  # this is tau times a quadratic form in independent N(0, I) draws with K x K
  # kernel
  #     C_ik = 1{i = k} a_i' B_i a_i  -  (B_i a_i)' B^{-1} (B_k a_k),
  # so E = tau tr(C), Var = 2 tau^2 tr(C^2), and the Satterthwaite df is
  #     df_BM = tr(C)^2 / tr(C^2).
  # Design-only (no residuals) and computable at the center from {B_i}.
  # Checks (2026-09-27): tr(C_KC) = (B^{-1})_jj (Theorem-3 exactness);
  # Monte Carlo mean/var/df match; balanced sites give df = K - 1 exactly.
  # ----------------------------------------------------------------------------
  B_sites <- map(valid, "B_site")
  bm_df <- function(T_list) {
    sapply(seq_len(p), function(j) {
      a <- lapply(T_list, function(Ti) crossprod(Ti, B_inv[, j]))
      d <- vapply(seq_len(K), function(i) drop(crossprod(a[[i]], B_sites[[i]] %*% a[[i]])),
                  numeric(1))
      W <- vapply(seq_len(K), function(i) drop(B_sites[[i]] %*% a[[i]]), numeric(p))
      W <- matrix(W, nrow = p)
      C <- diag(d, nrow = K) - crossprod(W, B_inv %*% W)
      trC2 <- sum(C * C)
      if (!is.finite(trC2) || trC2 < 1e-300) return(NA_real_)
      sum(diag(C))^2 / trC2
    })
  }
  T_id_list <- rep(list(diag(p)), K)

  # FG's own df (dtildeH) from saws: needs only the site scores and breads.
  fg_df <- rep(NA_real_, p)
  if (requireNamespace("saws", quietly = TRUE)) {
    u_mat <- do.call(rbind, lapply(S_sites, function(s) as.numeric(s)))
    om    <- array(0, c(K, p, p))
    for (i in seq_len(K)) om[i, , ] <- B_sites[[i]]
    # saws' df calculation is costly (~1 s per coefficient at K = 50, p = 20),
    # so compute it only for the coefficients the caller reports. Set with
    # options(fedgee.df_coefs = <names>); NULL = all coefficients.
    want <- getOption("fedgee.df_coefs", NULL)
    j_set <- if (is.null(want) || is.null(full_colnames)) seq_len(p) else which(full_colnames %in% want)
    fg_df[j_set] <- vapply(j_set, function(j) {
      tm <- matrix(0, 1, p); tm[1, j] <- 1
      out <- tryCatch(suppressWarnings(
        saws::saws(list(coefficients = beta_hat, u = u_mat, omega = om),
                   test = tm, method = "d5", bound = FG_BOUND)),
        error = function(e) NULL)
      if (is.null(out)) NA_real_ else as.numeric(out$df)
    }, numeric(1))
  }

  df_kp_site    <- K - p
  df_k1_site    <- K - 1
  df_kp_patient <- n_patients - p

  # ============================================================================
  # ASSEMBLE OUTPUT
  # ============================================================================
  make_variant <- function(vcov_mat, df_nu = rep(NA_real_, p), df_bm = rep(NA_real_, p),
                           df_fg = rep(NA_real_, p), df_kp = df_kp_site, df_k1 = df_k1_site) {
    list(vcov = vcov_mat, se = sqrt(diag(vcov_mat)),
         df_k1 = df_k1, df_kp = df_kp, df_nu = df_nu, df_bm = df_bm, df_fg = df_fg)
  }

  list(
    variants = list(
      patient_z = make_variant(vcov_patient, df_kp = df_kp_patient, df_k1 = n_patients - 1),
      site      = make_variant(vcov_site,   nu_df(S_sites),        bm_df(T_id_list)),
      kc        = make_variant(vcov_kc,     nu_df(S_kc_list),      bm_df(T_kc_list)),
      md        = make_variant(vcov_md,     nu_df(S_md_list),      bm_df(T_md_list)),
      fg        = make_variant(vcov_fg,     nu_df(S_fg_sites),     bm_df(T_fg_list), df_fg = fg_df),
      hc2obs    = make_variant(vcov_hc2obs, nu_df(S_hc2obs_sites), bm_df(T_id_list)),
      vd        = make_variant(vcov_vd, df_kp = df_kp_patient)
    ),
    diagnostics = list(
      n_sites          = K,
      n_patients       = n_patients,
      site_size_cv     = { ns <- map_dbl(valid, "n_clusters"); sd(ns) / mean(ns) },
      max_leverage     = max(vapply(B_sites, function(b) max(diag(b %*% B_inv)), numeric(1))),
      kc_min_eig_min   = min(kc_min_eigs, na.rm = TRUE),
      kc_min_eig_med   = median(kc_min_eigs, na.rm = TRUE),
      md_min_eig_min   = min(md_min_eigs, na.rm = TRUE),
      md_min_eig_med   = median(md_min_eigs, na.rm = TRUE),
      kc_clipped       = kc_clipped_count,
      md_clipped       = md_clipped_count,
      fg_capped        = fg_capped_count,
      hc2obs_clipped   = hc2obs_clipped_count,
      max_obs_leverage = max_h_obs,
      n_sites_vd       = n_sites_vd
    )
  )
}

###############################################################################
# 5. TRAINING LOOP v2
###############################################################################
train_FedGEE_v2 <- function(data_list,
                            initial_beta,
                            main_formula,
                            alpha,
                            family_obj,
                            corstr,
                            id_col,
                            site_betas_reduced = NULL,
                            reduced_colnames   = NULL,
                            n_iter = 50,
                            tol    = 1e-8,
                            verbose = FALSE) {

  N_sites   <- length(data_list)
  beta_curr <- initial_beta
  p         <- length(beta_curr)
  full_colnames <- if (!is.null(rownames(initial_beta))) rownames(initial_beta) else names(beta_curr)
  converged <- FALSE
  k <- 0L

  # ---- Iterate (NO corrections — just score equation) ----
  for (k in 1:n_iter) {
    site_outputs <- vector("list", N_sites)
    for (i in 1:N_sites) {
      if (is.null(alpha[[i]])) next
      site_outputs[[i]] <- get_site_stats_v2(
        site_data = data_list[[i]], beta_global = beta_curr,
        alpha = alpha[[i]], main_formula = main_formula,
        family_obj = family_obj, id_col = id_col, corstr = corstr,
        B_global = NULL, compute_H_obs = FALSE)
    }

    upd <- server_update(site_outputs, beta_curr)
    if (is.null(upd)) {
      if (verbose) cat("Aggregation failure.\n")
      return(NULL)
    }

    diff_norm <- sqrt(sum((upd$new_beta - beta_curr)^2))
    if (verbose) cat(sprintf("  iter %d | norm %.2e\n", k, diff_norm))
    beta_curr <- upd$new_beta

    if (diff_norm < tol) { converged <- TRUE; break }
    if (diff_norm > 1e4) {
      if (verbose) cat("Divergence detected.\n"); return(NULL)
    }
  }

  # ---- Final pass: compute primitives + FG using converged B_global ----
  # Use upd$B_total (from last iteration's aggregation) as B_global throughout.
  # At convergence, this equals the Bread re-evaluated at beta_curr to within tol.
  B_global_final <- upd$B_total

  site_outputs_final <- vector("list", N_sites)
  for (i in 1:N_sites) {
    if (is.null(alpha[[i]])) next
    site_outputs_final[[i]] <- get_site_stats_v2(
      site_data = data_list[[i]], beta_global = beta_curr,
      alpha = alpha[[i]], main_formula = main_formula,
      family_obj = family_obj, id_col = id_col, corstr = corstr,
      B_global = B_global_final, compute_H_obs = TRUE)
  }

  # Compute all variance variants
  all_var <- compute_all_variances(
    site_results = site_outputs_final,
    beta_hat = beta_curr,
    B_global = B_global_final,
    site_betas_reduced = site_betas_reduced,
    reduced_colnames = reduced_colnames,
    full_colnames = full_colnames
  )
  if (is.null(all_var)) return(NULL)

  list(
    coefficients = setNames(as.vector(beta_curr), full_colnames),
    variants     = all_var$variants,
    diagnostics  = all_var$diagnostics,
    iterations   = k,
    converged    = converged,
    B_global     = B_global_final,
    corstr       = corstr
  )
}

###############################################################################
# 6. CONVENIENCE WRAPPER
###############################################################################
fedgee_v2 <- function(data_list,
                      main_formula,
                      family_obj = binomial(link = "logit"),
                      corstr     = "independence",
                      id_col     = "pat_id",
                      n_iter     = 50,
                      tol        = 1e-8,
                      verbose    = FALSE) {

  prep <- Prep_FedGEE_v2(data_list, main_formula, family_obj, corstr, id_col)

  fit <- train_FedGEE_v2(
    data_list = data_list,
    initial_beta = prep$initial_values,
    main_formula = main_formula,
    alpha = prep$alpha_list,
    family_obj = family_obj,
    corstr = corstr,
    id_col = id_col,
    site_betas_reduced = prep$site_betas_reduced,
    reduced_colnames = prep$reduced_colnames,
    n_iter = n_iter, tol = tol, verbose = verbose)

  fit
}

###############################################################################
# 7. HELPER: extract one variant as a flat data frame for simulation reporting
###############################################################################
extract_variant_row <- function(fit, variant_name, coef_name) {
  empty <- list(beta_hat = NA_real_, se = NA_real_, df_k1 = NA_real_, df_kp = NA_real_,
                df_nu = NA_real_, df_bm = NA_real_, df_fg = NA_real_)
  if (is.null(fit) || is.null(fit$variants[[variant_name]])) return(empty)
  v <- fit$variants[[variant_name]]
  idx <- which(names(fit$coefficients) == coef_name)
  if (length(idx) == 0) return(empty)
  list(
    beta_hat = unname(fit$coefficients[idx]),
    se       = unname(v$se[idx]),
    df_k1    = v$df_k1,
    df_kp    = v$df_kp,
    df_nu    = v$df_nu[idx],
    df_bm    = v$df_bm[idx],
    df_fg    = v$df_fg[idx]
  )
}
