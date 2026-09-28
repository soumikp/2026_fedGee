###############################################################################
# Shared internals for centralized and decentralized Fed-GEE.
###############################################################################

# Working correlation matrix R_i(alpha) for one cluster of size n.
get_Ri <- function(corstr, phi, n) {
  if (corstr == "independence" || n == 1L) {
    return(diag(1, n))
  }
  if (corstr == "exchangeable") {
    Ri <- matrix(as.numeric(phi)[1], n, n)
    diag(Ri) <- 1
    return(Ri)
  }
  if (corstr == "ar1") {
    return(as.numeric(phi)[1]^abs(outer(seq_len(n), seq_len(n), "-")))
  }
  if (corstr == "unstructured") {
    Ri <- diag(1, n)
    if (length(phi) == n * (n - 1) / 2) {
      Ri[lower.tri(Ri)] <- phi
      Ri <- Ri + t(Ri) - diag(1, n)
    }
    return(Ri)
  }
  diag(1, n)
}

###############################################################################
# Symmetric matrix power via eigendecomposition. M MUST BE SYMMETRIC (it is
# symmetrized here). Eigenvalues are floored at `tol` before exponentiation;
# the floor shrinks the correction toward the identity rather than producing
# extreme inflation when one site dominates. Reports whether the floor fired.
###############################################################################
.mat_pow_safe <- function(M, power = -0.5, tol = 1e-6) {
  M <- (M + t(M)) / 2
  if (!all(is.finite(M))) {
    return(list(mat = NULL, clipped = TRUE, min_eig = NA_real_))
  }
  eig <- tryCatch(eigen(M, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig)) {
    return(list(mat = NULL, clipped = TRUE, min_eig = NA_real_))
  }
  vals <- eig$values
  mat <- eig$vectors %*% diag(pmax(vals, tol)^power, nrow = length(vals)) %*%
    t(eig$vectors)
  list(mat = mat, clipped = any(vals < tol), min_eig = min(vals))
}

###############################################################################
# Score-space leverage correction via the SIMILARITY TRANSFORM.
#
#   A_i = B^{1/2} (I - G_i)^{-power} B^{-1/2},   G_i = B^{-1/2} B_i B^{-1/2}
#
# The naive leverage H_i = B_i B^{-1} is a product of two symmetric matrices
# and is NOT symmetric, so eigen(I - H_i, symmetric = TRUE) would silently
# discard the upper triangle. G_i is symmetric PSD with the SAME eigenvalues as
# H_i and satisfies sum_i G_i = I -- the identity that makes the KC correction
# exact under the working model M = phi * B.
#
#   power = 0.5 -> Kauermann-Carroll (unbiased under the working model)
#   power = 1.0 -> Mancl-DeRouen (over-corrects by 1 / (1 - g))
#
# B_half / B_neghalf are formed once by the caller and reused across sites.
# Returns the operator (identity on failure) plus clipping diagnostics.
###############################################################################
.leverage_op <- function(B_unit, B_half, B_neghalf, power, tol = 0.05) {
  p <- nrow(B_unit)
  G_i <- B_neghalf %*% B_unit %*% B_neghalf
  res <- .mat_pow_safe(diag(p) - G_i, power = -power, tol = tol)
  mat <- if (is.null(res$mat)) diag(p) else B_half %*% res$mat %*% B_neghalf
  list(mat = mat, clipped = res$clipped, min_eig = res$min_eig)
}

# Leverage power for a named correction.
.correction_power <- function(correction) {
  switch(correction,
    none = 0,
    KC = 0.5,
    MD = 1,
    stop("Unknown correction: ", correction)
  )
}

# Symmetrize, optionally ridge, and solve; NULL on failure.
.safe_solve <- function(A, b = NULL, ridge = 0) {
  A <- as.matrix(A)
  A <- (A + t(A)) / 2
  if (!all(is.finite(A))) {
    return(NULL)
  }
  if (!is.null(b) && !all(is.finite(b))) {
    return(NULL)
  }
  if (ridge > 0) A <- A + diag(ridge, nrow(A))
  tryCatch(if (is.null(b)) solve(A) else solve(A, b), error = function(e) NULL)
}

###############################################################################
# Per-cluster GEE bread and score at one site, at coefficient vector `beta`.
#
#   B_j = D_j' V_j^{-1} D_j,   S_j = D_j' V_j^{-1} (y_j - mu_j)
#
# Returns the site totals plus the per-cluster pieces (needed only for the
# patient-level sandwich), or NULL if no cluster could be evaluated.
###############################################################################
.cluster_stats <- function(site_data, beta, alpha, main_formula, family_obj,
                           id_col, corstr, keep_clusters = FALSE) {
  beta <- as.matrix(beta)
  y_name <- all.vars(main_formula)[1]
  X <- model.matrix(main_formula, data = site_data)
  y <- site_data[[y_name]]
  if (is.factor(y)) y <- as.numeric(y != levels(y)[1])
  ids <- site_data[[id_col]]
  rows_by_cluster <- split(seq_along(ids), factor(ids, levels = unique(ids)))

  n_c <- length(rows_by_cluster)
  B_list <- vector("list", n_c)
  S_list <- vector("list", n_c)

  for (idx in seq_len(n_c)) {
    rows <- rows_by_cluster[[idx]]
    X_j <- X[rows, , drop = FALSE]
    n_j <- length(rows)

    eta_j <- as.vector(X_j %*% beta)
    mu_j <- family_obj$linkinv(eta_j)
    var_mu_j <- pmax(family_obj$variance(mu_j), 1e-12)
    D_j <- family_obj$mu.eta(eta_j) * X_j
    r_j <- y[rows] - mu_j

    if (corstr == "independence" || n_j == 1L) {
      V_inv_r <- r_j / var_mu_j
      V_inv_D <- D_j / var_mu_j
    } else {
      A_half_j <- sqrt(var_mu_j)
      V_j <- (A_half_j %o% A_half_j) * get_Ri(corstr, alpha, n_j)
      V_inv_j <- tryCatch(solve(V_j), error = function(e) NULL)
      if (is.null(V_inv_j)) next
      V_inv_r <- V_inv_j %*% r_j
      V_inv_D <- V_inv_j %*% D_j
    }

    B_list[[idx]] <- crossprod(D_j, V_inv_D)
    S_list[[idx]] <- crossprod(D_j, V_inv_r)
  }

  keep <- !vapply(B_list, is.null, logical(1))
  if (!any(keep)) {
    return(NULL)
  }
  B_list <- B_list[keep]
  S_list <- S_list[keep]

  list(
    B_site = Reduce(`+`, B_list),
    S_site = Reduce(`+`, S_list),
    B_clusters = if (keep_clusters) B_list else NULL,
    S_clusters = if (keep_clusters) S_list else NULL,
    n_clusters = sum(keep),
    n_obs = sum(lengths(rows_by_cluster[keep]))
  )
}

###############################################################################
# Local working-correlation parameter alpha_i from a site-only geeglm fit on
# the locally estimable (reduced) formula. Falls back to 0 on failure.
###############################################################################
.local_alpha <- function(site_data, reduced_formula, family_obj, id_col, corstr) {
  if (corstr == "independence") {
    return(0)
  }
  # droplevels: geepack refuses (and prints) on factor levels absent locally.
  site_data <- droplevels(site_data[order(site_data[[id_col]]), , drop = FALSE])
  site_data[[".id_var"]] <- site_data[[id_col]]
  fit <- tryCatch(
    suppressWarnings(geepack::geeglm(
      formula = reduced_formula, data = site_data, family = family_obj,
      id = .id_var, corstr = corstr
    )),
    error = function(e) NULL
  )
  a <- if (is.null(fit)) numeric(0) else fit$geese$alpha
  a <- if (length(a) == 0 || !all(is.finite(a))) 0 else unname(a)
  .clamp_alpha(a, corstr, max(table(site_data[[id_col]])))
}

# Keep a local alpha inside the region where R_i(alpha) is positive definite.
# Small sites can return an exchangeable alpha at or below -1/(n - 1), which
# makes the working covariance singular or indefinite. Pull such values 5%
# back inside the boundary and flag it with attr(, "clamped").
.clamp_alpha <- function(a, corstr, n_max) {
  lo <- switch(corstr,
    exchangeable = if (n_max > 1) -0.95 / (n_max - 1) else -Inf,
    ar1 = -0.95,
    -Inf
  )
  hi <- if (corstr %in% c("exchangeable", "ar1")) 0.95 else Inf
  out <- pmin(pmax(a, lo), hi)
  if (any(out != a)) attr(out, "clamped") <- TRUE
  out
}

.warn_clamped <- function(alpha) {
  hit <- names(alpha)[vapply(alpha, function(a) isTRUE(attr(a, "clamped")), logical(1))]
  if (length(hit)) {
    warning(
      "Local working-correlation alpha was outside the positive-definite range ",
      "and was clamped at site(s): ", paste(hit, collapse = ", "), "."
    )
  }
}

# Build the reduced formula dropping every term whose columns are ALL flagged
# site-constant. Terms are mapped via attr(, "assign") so factors survive.
.reduced_formula <- function(main_formula, X_example, constant_cols, data) {
  tt <- terms(main_formula, data = data)
  term_labels <- attr(tt, "term.labels")
  assign <- attr(X_example, "assign")
  is_const <- colnames(X_example) %in% constant_cols
  drop_term <- vapply(seq_along(term_labels), function(k) {
    cols_k <- which(assign == k)
    length(cols_k) > 0 && all(is_const[cols_k])
  }, logical(1))
  keep <- term_labels[!drop_term]
  reformulate(
    termlabels = if (length(keep)) keep else "1",
    response = main_formula[[2L]],
    intercept = attr(tt, "intercept") == 1L,
    env = environment(main_formula)
  )
}
