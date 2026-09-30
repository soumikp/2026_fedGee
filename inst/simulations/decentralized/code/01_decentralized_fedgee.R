###############################################################################
# Decentralized Federated GEE (Dec-Fed-GEE)
#
# Estimation logic:
#   1. Prep_DecFedGEE() uses a strictly decentralized consensus meta-GLM warm
#      start. Each site fits a local GLM using the locally identifiable reduced
#      formula, constructs its model-based precision matrix P_i = V_i^{-1} and
#      precision-weighted coefficient vector q_i = P_i beta_i, and exchanges
#      only these summaries through L_init rounds of neighbor-only consensus.
#      Site i initializes beta from its consensus summaries as P_tilde_i^{-1}
#      q_tilde_i. Site-level coefficients are initialized at zero and estimated
#      later through the full decentralized GEE iterations. Local working-
#      correlation parameters alpha_i are estimated using only site i's data.
#      Prep_FedGEE() is retained only for backward compatibility with
#      centralized Fed-GEE code.
#   2. At each outer iteration:
#        a. Run L_beta rounds of neighbor-only consensus on beta.
#        b. At each site's consensus beta, compute its local GEE Score and Bread
#           using get_site_stats(), reused without modification.
#        c. Run L_S rounds of Score consensus and L_B rounds of Bread consensus.
#        d. Update beta once by the consensus Fisher-scoring direction.
#        e. Monitor update error, neighbor-consensus error, and the norm of the
#           consensus Fisher-scoring direction.
#   3. After stopping, run a final L_beta-round beta consensus.
#
# Inference logic:
#   1. Recompute each site's local Bread and Meat at the fixed final beta.
#   2. Run L_B rounds of consensus on Bread and L_S rounds on Meat.
#   3. Construct the decentralized sandwich covariance from the consensus
#      averages. The default sandwich level is patient, so the Meat is the sum
#      of patient-level score outer products within each site.
#
# Main controls in decentralized_fedgee():
#   L_init  number of consensus rounds for meta-GLM precision summaries
#   L_beta  number of consensus rounds for beta
#   L_S     number of consensus rounds for Score and final Meat
#   L_B     number of consensus rounds for Bread
#
# WEIGHT MATRIX W:
#   1. Pass W directly; or
#   2. build hub/ring/VISN weights through the structure argument; or
#   3. use an object named W in .GlobalEnv for backward compatibility.
#   Rows and columns of W must follow the site order in data_list.
#
# SITE-LEVEL COVARIATES:
#   Pass site_constant_cols = c(...) to decentralized_fedgee() to name any
#   covariate that is constant within every site by study design (e.g. a
#   hospital-region indicator). This is supplied as metadata rather than
#   auto-detected across data_list, because auto-detecting "constant at every
#   site" requires the same kind of pre-loop, all-sites-at-once visibility
#   that the decentralized protocol is designed to avoid. Site-level columns
#   are initialized at 0 in Prep_DecFedGEE() and estimated purely through the
#   consensus + Fisher-scoring iterations (between-site variation).
###############################################################################
library(geepack)
library(dplyr)
library(purrr)
library(tibble)
library(Matrix)

###############################################################################
# 0. WEIGHT MATRIX CONSTRUCTION (hub-and-spoke / ring / VISN)
#
# Metropolis-Hastings rule using SELF-INCLUSIVE degree: d_i' = (# neighbors) + 1
# (the self-loop counts as one connection), matching Fig. 1 of Gu & Chen (2024):
#   w_ij = 1 / max(d_i', d_j')   for edges (i != j)
#   w_ii = 1 - sum_{j != i} w_ij
#
# Usage:
#   build_weight_matrix(K = 5, structure = "ring", K_neigh = 2)
#   build_weight_matrix(K = 5, structure = "hub", hub = 1)
#   build_weight_matrix(K = 8, structure = "visn",
#                        region_id = rep(1:4, each = 2),
#                        hub_sites = c(1, 3, 5, 7))
#   build_weight_matrix(K = 8, structure = "complete")
###############################################################################

# A must have 1 on the diagonal (self-loop) and 1 for each edge.
# Degree = rowSums(A), self-inclusive by construction.
mh_weights <- function(A) {
  K <- nrow(A)
  stopifnot(nrow(A) == ncol(A))
  stopifnot(all(diag(A) == 1))       # self-loops required
  stopifnot(isSymmetric(A))
  
  deg <- rowSums(A)                  # self-inclusive degree
  
  W <- matrix(0, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      if (i != j && A[i, j] == 1) {
        W[i, j] <- 1 / max(deg[i], deg[j])
      }
    }
  }
  diag(W) <- 1 - rowSums(W)
  
  stopifnot(isSymmetric(W))
  stopifnot(all(abs(rowSums(W) - 1) < 1e-10))
  W
}

# Hub-and-spoke: `hub` connects to every other site; spokes connect only to hub.
build_adj_hub <- function(K, hub = 1) {
  A <- diag(1, K)
  others <- setdiff(seq_len(K), hub)
  A[hub, others] <- 1
  A[others, hub] <- 1
  A
}

# Ring: each site connects to K_neigh nearest neighbors (circularly, K_neigh/2
# on each side). K_neigh must be even and < K.
build_adj_ring <- function(K, K_neigh = 2) {
  stopifnot(K_neigh %% 2 == 0, K_neigh < K)
  A <- diag(1, K)
  half <- K_neigh / 2
  for (i in seq_len(K)) {
    for (d in seq_len(half)) {
      j1 <- ((i - 1 + d) %% K) + 1
      j2 <- ((i - 1 - d) %% K) + 1
      A[i, j1] <- 1
      A[i, j2] <- 1
    }
  }
  A
}

# VISN: sites within the same region are fully connected to each other
# (including their region's hub); regional hubs are additionally fully
# connected to each other. Non-hub sites in different regions have no edge.
#   region_id: length-K vector, which region each site belongs to
#   hub_sites: indices (into 1:K) of each region's hub
build_adj_visn <- function(region_id, hub_sites) {
  K <- length(region_id)
  A <- diag(1, K)
  for (r in unique(region_id)) {
    idx <- which(region_id == r)
    A[idx, idx] <- 1                 # fully connected within region
  }
  A[hub_sites, hub_sites] <- 1       # hubs fully connected to each other
  A
}

# Complete graph: every site connects to every other site.
# Resulting W has every entry equal to 1/K (diagonal included).
build_adj_complete <- function(K) {
  matrix(1, K, K)
}

# Main entry point: build W for a given topology.
build_weight_matrix <- function(K,
                                structure = c("hub", "ring", "visn", "complete"),
                                hub = 1,
                                K_neigh = 2,
                                region_id = NULL,
                                hub_sites = NULL,
                                verbose = TRUE) {
  structure <- match.arg(structure)
  
  A <- switch(structure,
              hub      = build_adj_hub(K, hub = hub),
              ring     = build_adj_ring(K, K_neigh = K_neigh),
              complete = build_adj_complete(K),
              visn = {
                if (is.null(region_id) || is.null(hub_sites)) {
                  stop("VISN structure requires `region_id` (length K) and `hub_sites`.")
                }
                stopifnot(length(region_id) == K)
                build_adj_visn(region_id, hub_sites)
              }
  )
  
  W <- mh_weights(A)
  
  if (verbose) message(sprintf("[%s] K = %d", structure, K))
  
  list(W = W, A = A, structure = structure)
}

###############################################################################
# 0b. BASIC HELPERS
###############################################################################

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
    exponent <- abs(
      matrix(seq_len(n), nrow = n, ncol = n, byrow = TRUE) - seq_len(n)
    )
    return(as.numeric(phi)[1]^exponent)
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
# Symmetric matrix power, with an eigenvalue floor for numerical safety.
# M MUST be symmetric; the caller enforces this.
###############################################################################
.sym_pow_dec <- function(M, power, floor_val = 1e-10) {
  M <- (M + t(M)) / 2
  if (!all(is.finite(M))) return(NULL)
  eig <- tryCatch(eigen(M, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig)) return(NULL)
  vals <- pmax(eig$values, floor_val)^power
  eig$vectors %*% diag(vals, nrow = length(vals)) %*% t(eig$vectors)
}

###############################################################################
# Score-space leverage operator via the similarity transform.
#
#   G_i = B^{-1/2} B_i B^{-1/2}   symmetric, same eigenvalues as B_i B^{-1},
#                                 and sum_i G_i = I
#   A_i = B^{1/2} (I - G_i)^{-power} B^{-1/2}
#
# power = 0.5 -> Kauermann-Carroll;  power = 1.0 -> Mancl-DeRouen.
# Returns NULL if B is not invertible, signalling the caller to fall back to
# the uncorrected meat.
#
# NOTE: in the decentralized setting each site supplies its OWN consensus
# estimate of B (K * Bbar_i^(L)), so sum_i G_i = I holds only to O(rho^L).
###############################################################################
.dec_leverage_op <- function(B_unit, B_global, power = 0.5, floor_val = 1e-6) {
  p <- nrow(B_unit)
  B_half    <- .sym_pow_dec(B_global,  0.5)
  B_neghalf <- .sym_pow_dec(B_global, -0.5)
  if (is.null(B_half) || is.null(B_neghalf)) return(NULL)
  G_i  <- B_neghalf %*% B_unit %*% B_neghalf
  core <- .sym_pow_dec(diag(p) - G_i, -power, floor_val = floor_val)
  if (is.null(core)) return(NULL)
  B_half %*% core %*% B_neghalf
}

safe_solve <- function(A, b = NULL, ridge = 0) {
  A <- as.matrix(A)
  A <- (A + t(A)) / 2

  if (!all(is.finite(A))) return(NULL)
  if (!is.null(b) && !all(is.finite(b))) return(NULL)
  
  if (ridge > 0) {
    A <- A + diag(ridge, nrow(A))
  }
  
  tryCatch(
    {
      if (is.null(b)) solve(A) else solve(A, b)
    },
    error = function(e) NULL
  )
}

resolve_weight_matrix <- function(W, data_list, tolerance = 1e-8) {
  if (is.null(W)) {
    if (exists("W", envir = .GlobalEnv, inherits = FALSE)) {
      W <- get("W", envir = .GlobalEnv, inherits = FALSE)
    } else {
      stop(
        "Weight matrix W has not been defined. Define W in the simulation ",
        "script or pass W = W to fedgee()."
      )
    }
  }
  
  W <- as.matrix(W)
  K <- length(data_list)
  site_names <- names(data_list)
  if (is.null(site_names)) site_names <- as.character(seq_len(K))
  
  if (!all(dim(W) == c(K, K))) {
    stop(sprintf("W must be a %d x %d matrix because data_list has %d sites.",
                 K, K, K))
  }
  
  # If W is named, reorder it to match data_list.
  if (!is.null(rownames(W)) && !is.null(colnames(W)) &&
      setequal(rownames(W), site_names) && setequal(colnames(W), site_names)) {
    W <- W[site_names, site_names, drop = FALSE]
  }
  
  if (any(!is.finite(W))) stop("W contains non-finite entries.")
  if (any(W < -tolerance)) stop("W contains negative weights.")
  
  row_error <- max(abs(rowSums(W) - 1))
  if (row_error > tolerance) {
    stop("Rows of W must sum to 1. Maximum row-sum error: ", signif(row_error, 4))
  }
  
  col_error <- max(abs(colSums(W) - 1))
  if (col_error > 1e-6) {
    warning(
      "W is row-stochastic but not doubly stochastic. Ordinary gossip will ",
      "not generally converge to the equally weighted network average."
    )
  }
  
  symmetry_error <- max(abs(W - t(W)))
  if (symmetry_error > 1e-6) {
    warning(
      "W is not symmetric. The code can run, but the usual symmetric ",
      "average-consensus theory will not apply directly."
    )
  }
  
  dimnames(W) <- list(site_names, site_names)
  W
}

# Resolve W for decentralized_fedgee(), in priority order:
#   1. Explicit W (validated via resolve_weight_matrix)
#   2. structure = "hub"/"ring"/"visn" (built via build_weight_matrix, then
#      still passed through resolve_weight_matrix for the same validation
#      a hand-built W would get)
#   3. Fallback: look for a pre-built W in .GlobalEnv (backward
#      compatibility with older simulation scripts).
get_weight_matrix <- function(W, structure, data_list,
                              hub = 1, K_neigh = 2,
                              region_id = NULL, hub_sites = NULL,
                              verbose = TRUE) {
  if (!is.null(W)) {
    return(resolve_weight_matrix(W, data_list))
  }
  
  if (!is.null(structure)) {
    structure <- match.arg(structure, c("hub", "ring", "visn", "complete"))
    K <- length(data_list)
    wm <- build_weight_matrix(
      K = K, structure = structure, hub = hub, K_neigh = K_neigh,
      region_id = region_id, hub_sites = hub_sites, verbose = verbose
    )
    return(resolve_weight_matrix(wm$W, data_list))
  }
  
  # Neither W nor structure supplied -- fall back to .GlobalEnv (this call
  # will itself raise an informative error if nothing is found there).
  resolve_weight_matrix(NULL, data_list)
}

# Reshape a list of K (p x p) matrices into a K x (p*p) matrix, one row per
# site, so a single neighbor exchange can be done as one matrix product W %*% rows.
matrix_list_to_rows <- function(matrix_list, p) {
  K <- length(matrix_list)
  out <- matrix(0, nrow = K, ncol = p * p)
  
  for (k in seq_len(K)) {
    if (!is.null(matrix_list[[k]])) {
      out[k, ] <- as.vector(matrix_list[[k]])
    }
  }
  out
}

rows_to_matrix_list <- function(x, p, parameter_names = NULL) {
  lapply(seq_len(nrow(x)), function(k) {
    matrix(
      x[k, ], nrow = p, ncol = p,
      dimnames = list(parameter_names, parameter_names)
    )
  })
}

score_list_to_rows <- function(score_list, p) {
  K <- length(score_list)
  out <- matrix(0, nrow = K, ncol = p)
  
  for (k in seq_len(K)) {
    if (!is.null(score_list[[k]])) {
      out[k, ] <- as.vector(score_list[[k]])
    }
  }
  out
}

###############################################################################
# 1. INITIALIZATION -- CENTRALIZED VERSION (kept for backward compatibility)
#
#    Prep_FedGEE() is the original pooled/centralized initializer. It scans
#    every site's design matrix to flag site-constant covariates, then
#    meta-analyzes every site's local GLM fit into ONE global initial value
#    before any site starts iterating. That pre-loop, all-sites-at-once
#    aggregation is appropriate for the CENTRALIZED / one-step Fed-GEE
#    algorithm (Section 2 of the theory note), where a coordinator legitimately
#    sees every site's summary statistics.
#
#    decentralized_fedgee() no longer calls this function. It is kept,
#    UNMODIFIED, only so existing centralized-Fed-GEE code that depends on it
#    keeps working. See Prep_DecFedGEE() below for the decentralized
#    replacement.
###############################################################################
Prep_FedGEE <- function(data_list,
                        main_formula,
                        family_obj,
                        corstr,
                        id_col) {
  
  N_sites <- length(data_list)
  y_name  <- all.vars(main_formula)[1]
  
  # ------------------------------------------------------------------
  # 1. Identify site-constant columns in the model matrix
  #    A column is "site-constant" if it has zero variance at ANY site.
  #    (E.g., 'teaching' is 1 everywhere at site A, 0 everywhere at site B)
  # ------------------------------------------------------------------
  full_X_example <- model.matrix(main_formula, data = data_list[[1]])
  all_colnames   <- colnames(full_X_example)
  p_full         <- length(all_colnames)
  
  # For each site, check which columns have zero variance
  site_constant_flags <- matrix(FALSE, nrow = N_sites, ncol = p_full)
  colnames(site_constant_flags) <- all_colnames
  
  for (i in seq_len(N_sites)) {
    X_i <- model.matrix(main_formula, data = data_list[[i]])
    col_var <- apply(X_i, 2, var)
    site_constant_flags[i, ] <- (col_var < 1e-15)
  }
  
  # A column is a true "site-level covariate" if it has zero variance at EVERY site.
  # (E.g., 'teaching' = 1 for all patients at hospital A, = 0 at hospital B)
  #
  # Columns that are incidentally constant at only SOME sites (e.g., 'female' at
  # an all-female hospital) are kept in the local formula. The local GEE at that
  # site may fail, but tryCatch handles it and the site is excluded from meta-analysis.
  site_level_constant <- apply(site_constant_flags, 2, all)
  # (Intercept) is always constant -- that's expected, not a problem
  site_level_constant["(Intercept)"] <- FALSE
  
  constant_cols <- names(which(site_level_constant))
  varying_cols  <- setdiff(all_colnames, c("(Intercept)", constant_cols))
  
  if (length(constant_cols) > 0) {
    message("FedGEE: Detected site-level covariates (constant within every site): ",
            paste(constant_cols, collapse = ", "),
            "\n  These will be initialized at 0 and estimated via between-site variation.")
  }
  
  # ------------------------------------------------------------------
  # 2. Build reduced formula for local GEE fits (excluding site-constant terms)
  #
  # Built from TERM LABELS, never design-matrix column names. model.matrix()
  # expands a factor `f` into f2, f3, ... which are not variables in the data,
  # so a formula built from column names fails with "object 'f2' not found"
  # and every local fit silently returns NULL. attr(,"assign") maps each design
  # column to its generating term; a term is dropped only when ALL of its
  # columns are site-level constant.
  # ------------------------------------------------------------------
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
    # All covariates are site-constant -- local fits are intercept-only
    reduced_formula <- as.formula(paste(y_name, "~ 1"))
  }
  
  # ------------------------------------------------------------------
  # 3. Fit local GLMs on reduced formula for initial beta values
  #    We use GLM (not GEE) here because:
  #    - We only need rough initial values; the iterative protocol refines them
  #    - GLM model-based SEs are always invertible (no sandwich singularity)
  #    - With few patients per site, GEE sandwich can be singular
  #
  #    Sites where the GLM has aliased (NA) coefficients are excluded.
  #    This happens when a factor level is absent or has zero outcome variation.
  # ------------------------------------------------------------------
  glm_local <- data_list %>%
    map(function(df) {
      fit <- tryCatch(
        suppressWarnings(glm(formula = reduced_formula, data = df, family = family_obj)),
        error = function(e) NULL
      )
      # Reject fits that didn't converge or have aliased coefficients
      if (is.null(fit)) return(NULL)
      if (!fit$converged) return(NULL)
      if (any(is.na(coef(fit)))) return(NULL)
      return(fit)
    })
  
  valid_idx <- !map_lgl(glm_local, is.null)
  glm_valid <- glm_local[valid_idx]
  
  n_failed <- sum(!valid_idx)
  if (n_failed > 0) {
    message(sprintf("FedGEE Prep: %d / %d sites had GLM issues (aliased coefficients, non-convergence, or errors). These are excluded from initialization but will participate in the iterative protocol.",
                    n_failed, N_sites))
  }
  
  if (length(glm_valid) == 0) {
    # Fallback: initialize all coefficients at 0
    message("FedGEE Prep: All local GLMs failed. Initializing all coefficients at 0.")
    initial_values <- matrix(0, nrow = p_full, ncol = 1)
    rownames(initial_values) <- all_colnames
    initial_se <- rep(NA_real_, p_full)
    names(initial_se) <- all_colnames
    
    # Still need alpha
    alpha <- map(data_list, ~ 0)
    alpha_full <- as.list(alpha)
    
    return(list(
      initial_values  = initial_values,
      initial_se      = initial_se,
      alpha_list      = alpha_full,
      valid_sites     = integer(0),
      constant_cols   = constant_cols,
      reduced_formula = reduced_formula
    ))
  }
  
  # ------------------------------------------------------------------
  # 4. Meta-analyze reduced coefficients using model-based (naive) SEs
  # ------------------------------------------------------------------
  beta_hat_reduced <- map(glm_valid, ~ as.matrix(coef(.x)))
  var_beta_reduced <- map(glm_valid, ~ as.matrix(vcov(.x)))
  inv_var          <- map(var_beta_reduced, ~ tryCatch(solve(.x), error = function(e) NULL))
  
  ok       <- !map_lgl(inv_var, is.null)
  beta_hat_reduced <- beta_hat_reduced[ok]
  inv_var          <- inv_var[ok]
  
  if (length(inv_var) == 0) {
    message("FedGEE Prep: All variance inversions failed. Initializing from first valid GLM.")
    # Use first valid GLM's coefficients as initial values (no meta-analysis)
    beta_reduced <- as.matrix(coef(glm_valid[[1]]))
    se_meta_reduced <- sqrt(diag(vcov(glm_valid[[1]])))
  } else {
    inv_var_beta <- map2(inv_var, beta_hat_reduced, ~ .x %*% .y)
    den_FE       <- Reduce(`+`, inv_var)
    num_FE       <- Reduce(`+`, inv_var_beta)
    beta_reduced <- solve(den_FE, num_FE)
    
    var_meta_reduced <- solve(den_FE)
    se_meta_reduced  <- sqrt(diag(var_meta_reduced))
  }
  
  # ------------------------------------------------------------------
  # 4b. Extract alpha from local GEE fits (only if corstr != "independence")
  #     This is a separate step because alpha extraction needs GEE,
  #     but initial beta does not.
  # ------------------------------------------------------------------
  if (corstr != "independence") {
    gee_local <- data_list %>%
      map(function(df) {
        df[[".id_var"]] <- df[[id_col]]
        tryCatch(
          geeglm(formula = reduced_formula, data = df, family = family_obj,
                 id = .id_var, corstr = corstr),
          error = function(e) NULL
        )
      })
    
    alpha <- map(gee_local, function(m) {
      if (is.null(m)) return(0)
      a <- m$geese$alpha
      if (length(a) == 0) return(0)
      return(a)
    })
  } else {
    alpha <- map(data_list, ~ 0)
  }
  
  # ------------------------------------------------------------------
  # 5. Assemble full initial beta vector
  #    Reduced coefficients go in their correct positions.
  #    Site-constant coefficients initialized at 0.
  # ------------------------------------------------------------------
  reduced_names <- rownames(beta_reduced)
  if (is.null(reduced_names)) reduced_names <- names(coef(glm_valid[[1]]))
  
  initial_values <- matrix(0, nrow = p_full, ncol = 1)
  rownames(initial_values) <- all_colnames
  
  initial_se <- rep(NA_real_, p_full)
  names(initial_se) <- all_colnames
  
  for (nm in reduced_names) {
    if (nm %in% all_colnames) {
      idx_full    <- which(all_colnames == nm)
      idx_reduced <- which(reduced_names == nm)
      initial_values[idx_full, 1] <- beta_reduced[idx_reduced, 1]
      initial_se[idx_full]        <- se_meta_reduced[idx_reduced]
    }
  }
  
  # ------------------------------------------------------------------
  # 6. Package alpha list (already computed in step 4b)
  #    Ensure every site has an alpha entry (use 0 for failed sites)
  # ------------------------------------------------------------------
  alpha_full <- vector("list", N_sites)
  for (i in seq_len(N_sites)) {
    alpha_full[[i]] <- if (i <= length(alpha)) alpha[[i]] else 0
  }
  
  return(list(
    initial_values  = initial_values,
    initial_se      = initial_se,
    alpha_list      = alpha_full,
    valid_sites     = which(valid_idx),
    constant_cols   = constant_cols,
    reduced_formula = reduced_formula
  ))
}

###############################################################################
# 1b. INITIALIZATION -- DECENTRALIZED CONSENSUS META-GLM
#
# At site i:
#   1. Fit a local GLM using reduced_formula.
#   2. Form the model-based precision matrix P_i = V_i^{-1} and the
#      precision-weighted coefficient vector q_i = P_i beta_i.
#   3. Run L_init rounds of neighbor-only average consensus on P_i and q_i.
#   4. Recover the site-specific initializer beta_i^(0) from
#         beta_i^(0) = P_tilde_i^{-1} q_tilde_i.
#
# A site with a failed/noninvertible local GLM contributes zeros to the
# initialization summaries, but it can still receive a valid initializer from
# its neighbors and participates normally in the subsequent Dec-Fed-GEE fit.
# If no usable initialization information reaches a site, that site falls back
# to the zero vector. Site-level covariates are excluded from reduced_formula,
# initialized at zero, and estimated later using the full main_formula.
###############################################################################
Prep_DecFedGEE <- function(
    data_list,
    main_formula,
    family_obj,
    corstr,
    id_col,
    W,
    L_init = 50L,
    site_constant_cols = character(0),
    verbose = TRUE
) {

  # ==========================================================================
  # 1. Validate inputs and set up site/model information
  # ==========================================================================

  if (!is.list(data_list) || length(data_list) == 0L) {
    stop("data_list must be a nonempty list with one data frame per site.")
  }

  N_sites <- length(data_list)

  site_names <- names(data_list)
  if (is.null(site_names)) {
    site_names <- as.character(seq_len(N_sites))
  }

  if (anyDuplicated(site_names)) {
    stop("Site names in data_list must be unique.")
  }

  if (length(L_init) != 1L ||
      is.na(L_init) ||
      L_init < 0L ||
      L_init != as.integer(L_init)) {
    stop("L_init must be one nonnegative integer.")
  }

  if (!inherits(main_formula, "formula")) {
    stop("main_formula must be an R formula.")
  }

  if (!is.character(id_col) || length(id_col) != 1L) {
    stop("id_col must be the name of one patient/cluster ID column.")
  }

  # Validate and, when names are available, reorder W to match data_list.
  W <- resolve_weight_matrix(W, data_list)

  response_name <- all.vars(main_formula)[1L]

  for (i in seq_len(N_sites)) {

    if (!is.data.frame(data_list[[i]])) {
      stop("Each element of data_list must be a data frame.")
    }

    if (!(response_name %in% names(data_list[[i]]))) {
      stop(
        "Response variable '",
        response_name,
        "' is missing from site ",
        site_names[i],
        "."
      )
    }

    if (!(id_col %in% names(data_list[[i]]))) {
      stop(
        "ID variable '",
        id_col,
        "' is missing from site ",
        site_names[i],
        "."
      )
    }
  }


  # ==========================================================================
  # 2. Construct the reduced formula for local GLM/GEE initialization
  #
  # main_formula defines the final decentralized GEE model.
  #
  # Site-level covariates are constant within an individual site and therefore
  # cannot be estimated from a local model. Any formula term involving a
  # variable named in site_constant_cols is removed from reduced_formula.
  #
  # These coefficients remain initialized at zero and are estimated later using
  # the full decentralized GEE score and bread.
  # ==========================================================================

  full_terms <- terms(
    main_formula,
    data = data_list[[1]]
  )

  term_labels <- attr(full_terms, "term.labels")
  include_intercept <- attr(full_terms, "intercept") == 1L

  predictor_variables <- all.vars(
    delete.response(full_terms)
  )

  unknown_constant_variables <- setdiff(
    site_constant_cols,
    predictor_variables
  )

  if (length(unknown_constant_variables) > 0L) {
    warning(
      "The following site_constant_cols are not present in main_formula: ",
      paste(unknown_constant_variables, collapse = ", ")
    )
  }

  # Remove every formula term that involves at least one site-level variable.
  term_uses_site_constant <- vapply(
    term_labels,
    function(term_label) {

      variables_in_term <- all.vars(
        as.formula(
          paste("~", term_label),
          env = environment(main_formula)
        )
      )

      any(variables_in_term %in% site_constant_cols)
    },
    logical(1)
  )

  constant_terms <- term_labels[term_uses_site_constant]
  varying_terms  <- term_labels[!term_uses_site_constant]

  response_expression <- paste(
    deparse(main_formula[[2L]], width.cutoff = 500L),
    collapse = ""
  )

  reduced_formula <- reformulate(
    termlabels = varying_terms,
    response = response_expression,
    intercept = include_intercept,
    env = environment(main_formula)
  )

  if (length(varying_terms) == 0L && !include_intercept) {
    stop(
      "After removing site-level terms, reduced_formula has no locally ",
      "estimable coefficient. Include an intercept or at least one ",
      "within-site covariate."
    )
  }


  # ==========================================================================
  # 3. Determine coefficient names for the full and reduced models
  # ==========================================================================

  full_X_example <- model.matrix(
    main_formula,
    data = data_list[[1]]
  )

  reduced_X_example <- model.matrix(
    reduced_formula,
    data = data_list[[1]]
  )

  full_parameter_names    <- colnames(full_X_example)
  reduced_parameter_names <- colnames(reduced_X_example)

  p_full    <- length(full_parameter_names)
  p_reduced <- length(reduced_parameter_names)

  if (!all(reduced_parameter_names %in% full_parameter_names)) {
    stop(
      "The reduced model matrix contains coefficient names that cannot be ",
      "matched to the full model matrix. Check factor levels and contrasts."
    )
  }

  # Map full model-matrix columns back to formula terms. This also handles
  # factor variables and interactions involving site-level covariates.
  full_assign <- attr(full_X_example, "assign")

  constant_term_indices <- match(
    constant_terms,
    term_labels
  )

  constant_term_indices <- constant_term_indices[
    !is.na(constant_term_indices)
  ]

  constant_cols <- full_parameter_names[
    full_assign %in% constant_term_indices
  ]

  if (length(constant_cols) > 0L && verbose) {
    message(
      "Prep_DecFedGEE: site-level coefficient(s) initialized at zero: ",
      paste(constant_cols, collapse = ", ")
    )
  }


  # ==========================================================================
  # 4. Fit local GLMs and construct precision-weighted local summaries
  #
  # At site i:
  #   beta_i = local GLM estimate
  #   V_i    = model-based covariance estimate
  #   P_i    = V_i^{-1}
  #   q_i    = P_i beta_i
  #
  # Sites with failed or noninvertible local GLMs contribute P_i = 0 and q_i = 0
  # to initialization consensus. They still participate in later GEE iterations.
  # ==========================================================================

  local_beta_by_site <- matrix(
    NA_real_,
    nrow = N_sites,
    ncol = p_reduced,
    dimnames = list(site_names, reduced_parameter_names)
  )

  local_precision_list <- lapply(
    seq_len(N_sites),
    function(i) {
      matrix(
        0,
        nrow = p_reduced,
        ncol = p_reduced,
        dimnames = list(
          reduced_parameter_names,
          reduced_parameter_names
        )
      )
    }
  )

  local_weighted_beta <- matrix(
    0,
    nrow = N_sites,
    ncol = p_reduced,
    dimnames = list(site_names, reduced_parameter_names)
  )

  site_status <- rep(
    "local_glm_excluded_from_consensus",
    N_sites
  )
  names(site_status) <- site_names

  valid_local_glm <- rep(FALSE, N_sites)

  for (i in seq_len(N_sites)) {

    fit_i <- tryCatch(
      suppressWarnings(
        glm(
          formula = reduced_formula,
          data = data_list[[i]],
          family = family_obj
        )
      ),
      error = function(e) NULL
    )

    if (is.null(fit_i) || !isTRUE(fit_i$converged)) {
      next
    }

    beta_i <- coef(fit_i)

    # Reject aliased, missing, or nonfinite local coefficients.
    if (length(beta_i) != p_reduced ||
        is.null(names(beta_i)) ||
        !all(reduced_parameter_names %in% names(beta_i)) ||
        any(!is.finite(beta_i))) {
      next
    }

    beta_i <- beta_i[reduced_parameter_names]

    variance_i <- tryCatch(
      as.matrix(vcov(fit_i)),
      error = function(e) NULL
    )

    if (is.null(variance_i) ||
        !all(dim(variance_i) == c(p_reduced, p_reduced)) ||
        is.null(rownames(variance_i)) ||
        is.null(colnames(variance_i)) ||
        !all(reduced_parameter_names %in% rownames(variance_i)) ||
        !all(reduced_parameter_names %in% colnames(variance_i))) {
      next
    }

    variance_i <- variance_i[
      reduced_parameter_names,
      reduced_parameter_names,
      drop = FALSE
    ]

    if (any(!is.finite(variance_i))) {
      next
    }

    precision_i <- safe_solve(
      variance_i,
      ridge = 0
    )

    if (is.null(precision_i) ||
        any(!is.finite(precision_i))) {
      next
    }

    precision_i <- (
      precision_i + t(precision_i)
    ) / 2

    weighted_beta_i <- precision_i %*% matrix(
      beta_i,
      ncol = 1L,
      dimnames = list(reduced_parameter_names, NULL)
    )

    if (any(!is.finite(weighted_beta_i))) {
      next
    }

    local_beta_by_site[i, ] <- beta_i
    local_precision_list[[i]] <- precision_i
    local_weighted_beta[i, ] <- as.vector(weighted_beta_i)

    valid_local_glm[i] <- TRUE
    site_status[i] <- "local_glm_contributed_to_consensus"
  }

  n_valid_local_glm <- sum(valid_local_glm)

  if (verbose) {
    message(
      sprintf(
        paste0(
          "Prep_DecFedGEE: %d / %d local GLMs contributed to ",
          "the consensus meta-GLM initializer."
        ),
        n_valid_local_glm,
        N_sites
      )
    )
  }


  # ==========================================================================
  # 5. Run consensus on precision matrices and weighted coefficient vectors
  #
  # Average consensus gives each site approximations to:
  #   P_bar = (1/K) sum_i P_i
  #   q_bar = (1/K) sum_i q_i
  #
  # The factor 1/K cancels in beta_i^(0) = P_bar_i^{-1} q_bar_i.
  # ==========================================================================

  precision_rows <- matrix_list_to_rows(
    local_precision_list,
    p = p_reduced
  )

  consensus_precision_rows <- consensus_rounds(
    values = precision_rows,
    W = W,
    L = L_init,
    row_names = site_names
  )

  consensus_weighted_beta <- consensus_rounds(
    values = local_weighted_beta,
    W = W,
    L = L_init,
    row_names = site_names,
    col_names = reduced_parameter_names
  )

  consensus_precision_list <- rows_to_matrix_list(
    consensus_precision_rows,
    p = p_reduced,
    parameter_names = reduced_parameter_names
  )


  # ==========================================================================
  # 6. Recover each site's consensus meta-GLM initial coefficient vector
  #
  # A small scale-adjusted ridge is used only when the consensus precision
  # matrix cannot be inverted directly. If no usable initialization information
  # reaches a site, that site falls back to zero initialization.
  # ==========================================================================

  initial_values <- matrix(
    0,
    nrow = N_sites,
    ncol = p_full,
    dimnames = list(site_names, full_parameter_names)
  )

  reduced_positions_in_full <- match(
    reduced_parameter_names,
    full_parameter_names
  )

  consensus_fallback <- rep(FALSE, N_sites)

  for (i in seq_len(N_sites)) {

    precision_i <- consensus_precision_list[[i]]
    precision_i <- (precision_i + t(precision_i)) / 2

    weighted_beta_i <- matrix(
      consensus_weighted_beta[i, ],
      ncol = 1L,
      dimnames = list(reduced_parameter_names, NULL)
    )

    has_consensus_information <- (
      all(is.finite(precision_i)) &&
      all(is.finite(weighted_beta_i)) &&
      max(abs(precision_i)) > sqrt(.Machine$double.eps)
    )

    if (!has_consensus_information) {
      consensus_fallback[i] <- TRUE
      site_status[i] <- "consensus_meta_glm_fallback_zero"
      next
    }

    initial_beta_i <- safe_solve(
      precision_i,
      weighted_beta_i,
      ridge = 0
    )

    # If direct inversion fails, use a small scale-adjusted ridge.
    if (is.null(initial_beta_i)) {

      precision_scale <- max(
        1,
        max(abs(diag(precision_i)))
      )

      initialization_ridge <- 1e-8 * precision_scale

      initial_beta_i <- safe_solve(
        precision_i,
        weighted_beta_i,
        ridge = initialization_ridge
      )
    }

    if (is.null(initial_beta_i) ||
        any(!is.finite(initial_beta_i))) {
      consensus_fallback[i] <- TRUE
      site_status[i] <- "consensus_meta_glm_fallback_zero"
      next
    }

    initial_values[
      i,
      reduced_positions_in_full
    ] <- as.vector(initial_beta_i)
  }

  n_consensus_fallback <- sum(consensus_fallback)

  if (verbose) {
    message(
      sprintf(
        paste0(
          "Prep_DecFedGEE: consensus meta-GLM initialization completed; ",
          "L_init = %d, zero fallback used at %d / %d sites."
        ),
        as.integer(L_init),
        n_consensus_fallback,
        N_sites
      )
    )
  }


  # ==========================================================================
  # 7. Estimate local working-correlation parameters
  #
  # alpha_i is estimated using only site i's data and reduced_formula.
  # If local GEE estimation fails, alpha_i falls back to zero.
  # ==========================================================================

  alpha_list <- vector(
    mode = "list",
    length = N_sites
  )
  names(alpha_list) <- site_names

  for (i in seq_len(N_sites)) {

    if (corstr == "independence") {
      alpha_list[[i]] <- 0
      next
    }

    site_data_i <- data_list[[i]]
    site_data_i[[".id_var"]] <- site_data_i[[id_col]]

    gee_i <- tryCatch(
      suppressWarnings(
        geeglm(
          formula = reduced_formula,
          data = site_data_i,
          family = family_obj,
          id = .id_var,
          corstr = corstr
        )
      ),
      error = function(e) NULL
    )

    alpha_i_is_valid <- (
      !is.null(gee_i) &&
      length(gee_i$geese$alpha) > 0L &&
      all(is.finite(gee_i$geese$alpha))
    )

    if (alpha_i_is_valid) {
      alpha_list[[i]] <- gee_i$geese$alpha
    } else {
      alpha_list[[i]] <- 0
    }
  }


  # ==========================================================================
  # 8. Return initialization results
  # ==========================================================================

  list(
    initial_values_by_site = initial_values,
    alpha_list = alpha_list,
    constant_cols = constant_cols,
    reduced_formula = reduced_formula,
    site_status = site_status,
    L_init = as.integer(L_init)
  )
}

###############################################################################
# 2. LOCAL GEE COMPUTATION AT ONE SITE
###############################################################################

get_site_stats <- function(site_data,
                           beta_global,
                           alpha,
                           main_formula,
                           family_obj,
                           id_col,
                           corstr,
                           sandwich_level = "site",
                           md_correction  = FALSE,
                           B_global       = NULL) {
  
  beta_curr <- as.matrix(beta_global)
  y_name    <- all.vars(main_formula)[1]
  
  # Full design matrix and response for the site
  Full_X <- model.matrix(main_formula, data = site_data)
  Full_y <- site_data[[y_name]]
  p      <- ncol(Full_X)
  
  # Identify unique clusters (patients)
  cluster_IDs <- unique(site_data[[id_col]])
  
  # Accumulators
  Bread_clusters <- vector("list", length(cluster_IDs))
  Score_clusters <- vector("list", length(cluster_IDs))
  
  idx <- 0L
  for (j in cluster_IDs) {
    idx <- idx + 1L
    
    row_idx <- which(site_data[[id_col]] == j)
    X_j     <- Full_X[row_idx, , drop = FALSE]
    y_j     <- Full_y[row_idx]
    n_j     <- length(y_j)
    
    # Linear predictor and link quantities
    eta_j      <- as.vector(X_j %*% beta_curr)
    mu_j       <- family_obj$linkinv(eta_j)
    dmu_deta_j <- family_obj$mu.eta(eta_j)
    var_mu_j   <- pmax(family_obj$variance(mu_j), 1e-12)
    
    # D_j = diag(dmu/deta) %*% X
    D_j <- dmu_deta_j * X_j  # n_j x p, element-wise column scaling
    
    # Working covariance: V = A^{1/2} R A^{1/2}
    A_half_j <- sqrt(var_mu_j)
    
    if (corstr == "independence" || n_j == 1) {
      # V_inv is diagonal — avoid full matrix ops
      V_inv_diag <- 1 / var_mu_j
      V_inv_r    <- V_inv_diag * (y_j - mu_j)
      V_inv_D    <- V_inv_diag * D_j
    } else {
      R_j <- get_Ri(corstr, alpha, n_j)
      # V = diag(A_half) %*% R %*% diag(A_half)
      V_j <- (A_half_j %o% A_half_j) * R_j
      V_inv_j <- tryCatch(solve(V_j), error = function(e) NULL)
      if (is.null(V_inv_j)) next
      V_inv_r <- V_inv_j %*% (y_j - mu_j)
      V_inv_D <- V_inv_j %*% D_j
    }
    
    # Bread_j = D' V^{-1} D
    B_j <- crossprod(D_j, V_inv_D)
    
    # Score_j = D' V^{-1} r
    S_j <- crossprod(D_j, V_inv_r)
    
    Bread_clusters[[idx]] <- B_j
    Score_clusters[[idx]] <- S_j
  }
  
  # Remove NULLs from skipped clusters
  keep <- !sapply(Bread_clusters, is.null)
  Bread_clusters <- Bread_clusters[keep]
  Score_clusters <- Score_clusters[keep]
  
  if (length(Bread_clusters) == 0) return(NULL)
  
  # ------ Aggregate bread and score ------
  B_site <- Reduce(`+`, Bread_clusters)
  S_site <- Reduce(`+`, Score_clusters)
  
  # ------ Compute meat based on sandwich_level ------
  if (sandwich_level == "site") {
    # Hospital-level sandwich: one outer product for entire site
    # Accounts for ALL within-hospital correlation (between + within patient)
    if (md_correction && !is.null(B_global)) {
      # Score-space leverage correction via the SIMILARITY TRANSFORM.
      #   G_i = B^{-1/2} B_i B^{-1/2}   (symmetric, same eigenvalues as
      #                                  B_i B^{-1}, and sum_i G_i = I)
      #   A_i = B^{1/2} (I - G_i)^{-1/2} B^{-1/2}
      # The naive leverage B_i B^{-1} is NOT symmetric, so calling
      # eigen(..., symmetric = TRUE) on I - H silently discards the upper
      # triangle and returns the wrong eigenvalues.
      A_i <- .dec_leverage_op(B_site, B_global, power = 0.5)
      M_site <- if (is.null(A_i)) tcrossprod(S_site) else tcrossprod(A_i %*% S_site)
    } else {
      M_site <- tcrossprod(S_site)
    }
    
  } else {
    # Patient-level sandwich: sum of per-patient outer products
    # Only accounts for within-patient correlation; assumes patients independent
    if (md_correction && !is.null(B_global)) {
      M_list <- lapply(seq_along(Score_clusters), function(i) {
        S_j <- Score_clusters[[i]]
        A_j <- .dec_leverage_op(Bread_clusters[[i]], B_global, power = 0.5)
        if (is.null(A_j)) tcrossprod(S_j) else tcrossprod(A_j %*% S_j)
      })
    } else {
      M_list <- lapply(Score_clusters, function(S_j) tcrossprod(S_j))
    }
    M_site <- Reduce(`+`, M_list)
  }
  
  return(list(
    Bread     = B_site,
    Score     = S_site,
    Meat      = M_site,
    n_clusters = sum(keep)
  ))
}

###############################################################################
# 3. CONSENSUS HELPERS
###############################################################################

# Apply L rounds of neighbor-only average consensus to a K x q object.
# Each row corresponds to one site. This works for parameters, scores, and
# vectorized matrices such as Bread and Meat.
consensus_rounds <- function(values, W, L, row_names = NULL, col_names = NULL) {
  if (length(L) != 1L || is.na(L) || L < 0 || L != as.integer(L)) {
    stop("Consensus rounds L must be one nonnegative integer.")
  }
  
  out <- as.matrix(values)
  if (nrow(out) != nrow(W)) {
    stop("The consensus object must have one row per site.")
  }
  
  if (L > 0L) {
    for (ell in seq_len(as.integer(L))) {
      out <- W %*% out
    }
  }
  
  if (!is.null(row_names)) rownames(out) <- row_names
  if (!is.null(col_names)) colnames(out) <- col_names
  out
}

# Allow a constant, vector, or function-valued Fisher-scoring step size.
resolve_step_size <- function(step_size, t) {
  eta_t <- if (is.function(step_size)) {
    step_size(t)
  } else if (length(step_size) == 1L) {
    step_size
  } else {
    if (t > length(step_size)) {
      stop("step_size is a vector but has fewer than n_iter entries.")
    }
    step_size[t]
  }
  
  if (length(eta_t) != 1L || !is.finite(eta_t) || eta_t <= 0) {
    stop("Each Fisher-scoring step size must be one positive finite number.")
  }
  as.numeric(eta_t)
}

###############################################################################
# 4. DECENTRALIZED TRAINING AND INFERENCE
#
# Estimation iteration:
#   1. L_beta rounds of parameter consensus.
#   2. Local GEE Score and Bread evaluated at the consensus parameter.
#   3. L_S rounds of Score consensus and L_B rounds of Bread consensus.
#   4. One Fisher-scoring update.
#   5. Check update, consensus, and Newton-direction residuals.
#
# Final inference:
#   1. Run a final L_beta-round parameter consensus.
#   2. Recompute local Bread and Meat at the fixed final estimates.
#   3. Apply L_B rounds to Bread and L_S rounds to Meat.
#   4. Construct the sandwich covariance from the consensus averages.
###############################################################################

train_DecentralizedFedGEE <- function(data_list,
                                      initial_beta_by_site,
                                      main_formula,
                                      alpha,
                                      family_obj,
                                      corstr,
                                      id_col,
                                      W,
                                      L_beta = 1L,
                                      L_S = 1L,
                                      L_B = 1L,
                                      n_iter = 100L,
                                      tol = 1e-8,
                                      tol_update = tol,
                                      tol_consensus = tol,
                                      tol_score = tol,
                                      step_size = 0.5,
                                      sandwich_level = "patient",
                                      md_correction = FALSE,
                                      ridge = 1e-5,
                                      verbose = TRUE) {
  stopifnot(sandwich_level %in% c("site", "patient"))
  
  K <- length(data_list)
  beta_curr <- as.matrix(initial_beta_by_site)
  parameter_names <- colnames(beta_curr)
  p <- ncol(beta_curr)
  site_names <- rownames(beta_curr)
  
  if (nrow(beta_curr) != K) {
    stop("initial_beta_by_site must have one row per site.")
  }
  if (length(alpha) != K) {
    stop("alpha must contain one entry per site.")
  }
  if (length(n_iter) != 1L || is.na(n_iter) || n_iter < 1L ||
      n_iter != as.integer(n_iter)) {
    stop("n_iter must be one positive integer.")
  }
  
  # Validate all consensus-round inputs before starting.
  invisible(consensus_rounds(matrix(0, K, 1), W, L_beta))
  invisible(consensus_rounds(matrix(0, K, 1), W, L_S))
  invisible(consensus_rounds(matrix(0, K, 1), W, L_B))
  
  tolerances <- c(
    update = tol_update,
    consensus = tol_consensus,
    score = tol_score
  )
  if (any(!is.finite(tolerances)) || any(tolerances < 0)) {
    stop("All convergence tolerances must be nonnegative finite numbers.")
  }
  
  if (verbose) {
    cat("Starting decentralized FedGEE\n")
    cat("  Family:", family_obj$family, "| Link:", family_obj$link, "\n")
    cat("  Working correlation:", corstr, "\n")
    cat("  Sandwich level:", sandwich_level, "\n")
    cat("  Consensus rounds: L_beta =", L_beta,
        "| L_S =", L_S, "| L_B =", L_B, "\n")
    cat("  Tolerances: update =", format(tol_update, scientific = TRUE),
        "| consensus =", format(tol_consensus, scientific = TRUE),
        "| score =", format(tol_score, scientific = TRUE), "\n")
    cat("  Max iterations:", n_iter, "\n")
  }
  
  beta_history <- list(beta_curr)
  diagnostic_history <- vector("list", n_iter)
  converged <- FALSE
  iterations_completed <- 0L
  
  for (t in seq_len(n_iter)) {
    eta_t <- resolve_step_size(step_size, t)
    
    # 1. Parameter consensus.
    beta_tilde <- consensus_rounds(
      beta_curr, W, L_beta,
      row_names = site_names,
      col_names = parameter_names
    )
    
    # 2. Local GEE quantities evaluated at the consensus parameter.
    site_outputs <- vector("list", K)
    for (i in seq_len(K)) {
      site_outputs[[i]] <- get_site_stats(
        site_data = data_list[[i]],
        beta_global = beta_tilde[i, ],
        alpha = alpha[[i]],
        main_formula = main_formula,
        family_obj = family_obj,
        id_col = id_col,
        corstr = corstr,
        sandwich_level = sandwich_level,
        md_correction = FALSE,
        B_global = NULL
      )
    }
    
    failed_stats <- which(vapply(site_outputs, is.null, logical(1)))
    if (length(failed_stats) > 0L) {
      if (verbose) {
        cat("D-FedGEE stopped: local GEE statistics failed at site(s):",
            paste(failed_stats, collapse = ", "), "\n")
      }
      return(NULL)
    }
    
    B_rows <- matrix_list_to_rows(lapply(site_outputs, `[[`, "Bread"), p)
    S_rows <- score_list_to_rows(lapply(site_outputs, `[[`, "Score"), p)
    
    # 3. Separate consensus procedures for Score and Bread.
    S_tilde_rows <- consensus_rounds(
      S_rows, W, L_S,
      row_names = site_names,
      col_names = parameter_names
    )
    B_tilde_rows <- consensus_rounds(
      B_rows, W, L_B,
      row_names = site_names
    )
    B_tilde_list <- rows_to_matrix_list(
      B_tilde_rows, p, parameter_names
    )
    
    # 4. Fisher-scoring update at every site.
    beta_next <- matrix(
      NA_real_, nrow = K, ncol = p,
      dimnames = list(site_names, parameter_names)
    )
    newton_direction <- matrix(
      NA_real_, nrow = K, ncol = p,
      dimnames = list(site_names, parameter_names)
    )
    
    failed_solve <- integer(0)
    for (i in seq_len(K)) {
      S_tilde_i <- matrix(
        S_tilde_rows[i, ], ncol = 1,
        dimnames = list(parameter_names, NULL)
      )
      direction_i <- safe_solve(
        B_tilde_list[[i]], S_tilde_i, ridge = ridge
      )
      if (is.null(direction_i)) {
        failed_solve <- c(failed_solve, i)
        next
      }
      
      newton_direction[i, ] <- as.vector(direction_i)
      beta_next[i, ] <- beta_tilde[i, ] + eta_t * as.vector(direction_i)
    }
    
    if (length(failed_solve) > 0L) {
      if (verbose) {
        cat("D-FedGEE stopped: consensus Bread was singular at site(s):",
            paste(failed_solve, collapse = ", "), "\n")
      }
      return(NULL)
    }
    
    if (!all(is.finite(beta_next)) || !all(is.finite(newton_direction))) {
      if (verbose) cat("D-FedGEE stopped because non-finite values were produced.\n")
      return(NULL)
    }
    
    # 5. Three convergence diagnostics.
    update_error <- sqrt(rowSums((beta_next - beta_curr)^2))
    
    neighbor_average <- W %*% beta_next
    dimnames(neighbor_average) <- list(site_names, parameter_names)
    consensus_error <- sqrt(rowSums((beta_next - neighbor_average)^2))
    
    score_error <- sqrt(rowSums(newton_direction^2))
    
    max_update <- max(update_error)
    max_consensus <- max(consensus_error)
    max_score <- max(score_error)
    
    diagnostic_history[[t]] <- data.frame(
      iteration = t,
      step_size = eta_t,
      max_update_error = max_update,
      max_consensus_error = max_consensus,
      max_score_error = max_score
    )
    
    beta_curr <- beta_next
    beta_history[[t + 1L]] <- beta_curr
    iterations_completed <- t
    
    if (verbose) {
      cat(sprintf(
        "  Iteration %d | update %.3e | consensus %.3e | score %.3e\n",
        t, max_update, max_consensus, max_score
      ))
    }
    
    if (max_update > 1e4 || max_score > 1e4) {
      if (verbose) cat("D-FedGEE divergence detected.\n")
      return(NULL)
    }
    
    if (max_update < tol_update &&
        max_consensus < tol_consensus &&
        max_score < tol_score) {
      converged <- TRUE
      if (verbose) cat("  D-FedGEE converged.\n")
      break
    }
  }
  
  diagnostic_history <- do.call(
    rbind, diagnostic_history[seq_len(iterations_completed)]
  )
  
  # Final parameter consensus. No additional Fisher-scoring update is made.
  beta_hat <- consensus_rounds(
    beta_curr, W, L_beta,
    row_names = site_names,
    col_names = parameter_names
  )
  
  # Recompute local Bread and patient/site-level Meat at the final fixed beta.
  local_final <- vector("list", K)
  for (i in seq_len(K)) {
    local_final[[i]] <- get_site_stats(
      site_data = data_list[[i]],
      beta_global = beta_hat[i, ],
      alpha = alpha[[i]],
      main_formula = main_formula,
      family_obj = family_obj,
      id_col = id_col,
      corstr = corstr,
      sandwich_level = sandwich_level,
      md_correction = FALSE,
      B_global = NULL
    )
  }
  
  failed_final <- which(vapply(local_final, is.null, logical(1)))
  if (length(failed_final) > 0L) {
    if (verbose) {
      cat("D-FedGEE inference failed at site(s):",
          paste(failed_final, collapse = ", "), "\n")
    }
    return(NULL)
  }
  
  # Consensus Bread is an approximation to K^{-1} sum_i B_i.
  B_final_rows <- matrix_list_to_rows(lapply(local_final, `[[`, "Bread"), p)
  B_bar_rows <- consensus_rounds(
    B_final_rows, W, L_B,
    row_names = site_names
  )
  B_bar_list <- rows_to_matrix_list(B_bar_rows, p, parameter_names)
  
  # Optional correction is recomputed using the approximate global Bread sum.
  # get_site_stats() is intentionally reused without modification.
  if (md_correction) {
    
    corrected_meat <- vector("list", K)
    
    for (i in seq_len(K)) {
      
      # B_bar_list[[i]] is the consensus approximation to
      # (1 / K) * sum_k B_k(beta_hat).
      # Therefore, multiply by K to recover the global Bread.
      B_global_i <- K * B_bar_list[[i]]
      B_global_i <- (B_global_i + t(B_global_i)) / 2
      
      corrected_i <- get_site_stats(
        site_data = data_list[[i]],
        beta_global = beta_hat[i, ],
        alpha = alpha[[i]],
        main_formula = main_formula,
        family_obj = family_obj,
        id_col = id_col,
        corstr = corstr,
        sandwich_level = sandwich_level,
        md_correction = TRUE,
        B_global = B_global_i
      )
      
      if (
        is.null(corrected_i) ||
        is.null(corrected_i$Meat) ||
        any(!is.finite(corrected_i$Meat))
      ) {
        if (verbose) {
          cat(
            "D-FedGEE leverage-adjusted Meat failed at site",
            i, "\n"
          )
        }
        return(NULL)
      }
      
      # get_site_stats() has already applied
      # (I - H)^(-1/2) to the score contribution.
      corrected_meat[[i]] <-
        (corrected_i$Meat + t(corrected_i$Meat)) / 2
    }
    
    # Replace the uncorrected local meats.
    for (i in seq_len(K)) {
      local_final[[i]]$Meat <- corrected_meat[[i]]
    }
  }
  
  # Meat is score-derived, so L_S controls its final consensus rounds.
  # M_bar approximates K^{-1} sum_i M_i.
  M_final_rows <- matrix_list_to_rows(lapply(local_final, `[[`, "Meat"), p)
  M_bar_rows <- consensus_rounds(
    M_final_rows, W, L_S,
    row_names = site_names
  )
  M_bar_list <- rows_to_matrix_list(M_bar_rows, p, parameter_names)
  
  # With average-consensus outputs,
  # Var(beta_hat) = K^{-1} B_bar^{-1} M_bar B_bar^{-T}.
  vcov_by_site <- vector("list", K)
  se_by_site <- matrix(
    NA_real_, nrow = K, ncol = p,
    dimnames = list(site_names, parameter_names)
  )
  
  for (i in seq_len(K)) {
    B_inv_i <- safe_solve(B_bar_list[[i]], ridge = ridge)
    if (is.null(B_inv_i)) next
    
    vcov_i <- (1 / K) * B_inv_i %*% M_bar_list[[i]] %*% t(B_inv_i)
    vcov_i <- (vcov_i + t(vcov_i)) / 2
    dimnames(vcov_i) <- list(parameter_names, parameter_names)
    
    vcov_by_site[[i]] <- vcov_i
    se_by_site[i, ] <- sqrt(pmax(diag(vcov_i), 0))
  }
  names(vcov_by_site) <- site_names
  
  # --------------------------------------------------------------------
  # Network-average decentralized estimator (Section 3.2 of the theory
  # note):
  #   beta_dec = (1/K) sum_i beta_hat_i
  # This is the estimator Theorems 5-6 establish consistency and
  # asymptotic normality for -- a single length-p vector, not a K x p
  # matrix of per-site rows. Because W is doubly stochastic, this equals
  # (1/K) sum_i beta_curr_i as well (the final consensus step only
  # redistributes mass across sites; it does not change the network
  # average). Per-site rows (coefficients_by_site) are retained only as
  # a diagnostic of residual cross-site disagreement, not as the
  # reported point estimate.
  #
  # vcov_dec is reported as the corresponding average of the per-site
  # sandwich covariances, giving a single p x p matrix (and se_dec its
  # square-root diagonal) that a site could reconstruct locally via one
  # additional round of neighbor-to-neighbor consensus, exactly as noted
  # in the theory note for beta_dec itself.
  # --------------------------------------------------------------------
  beta_dec <- setNames(colMeans(beta_hat), parameter_names)
  
  valid_vcov <- vcov_by_site[!vapply(vcov_by_site, is.null, logical(1))]
  vcov_dec <- if (length(valid_vcov) > 0) {
    v <- Reduce(`+`, valid_vcov) / length(valid_vcov)
    dimnames(v) <- list(parameter_names, parameter_names)
    v
  } else {
    matrix(NA_real_, p, p, dimnames = list(parameter_names, parameter_names))
  }
  se_dec <- setNames(sqrt(pmax(diag(vcov_dec), 0)), parameter_names)
  
  n_clusters_local <- vapply(local_final, function(z) z$n_clusters, numeric(1))
  n_valid_sites <- sum(n_clusters_local > 0)
  n_patients <- sum(n_clusters_local)
  
  if (sandwich_level == "site") {
    total_clusters <- n_valid_sites
    df_residual <- n_valid_sites - p
  } else {
    total_clusters <- n_patients
    df_residual <- n_patients - p
  }
  
  out <- list(
    coefficients = beta_dec,
    vcov = vcov_dec,
    se = se_dec,
    coefficients_by_site = beta_hat,
    coefficients_before_final_consensus = beta_curr,
    se_by_site = se_by_site,
    vcov_by_site = vcov_by_site,
    Bread_by_site = B_bar_list,
    Meat_by_site = M_bar_list,
    local_Bread_at_final = lapply(local_final, `[[`, "Bread"),
    local_Meat_at_final = lapply(local_final, `[[`, "Meat"),
    df_residual = df_residual,
    total_clusters = total_clusters,
    n_sites = K,
    n_valid_sites = n_valid_sites,
    n_patients = n_patients,
    iterations = iterations_completed,
    converged = converged,
    history = beta_history,
    diagnostic_history = diagnostic_history,
    final_diagnostics = if (nrow(diagnostic_history) > 0L) {
      diagnostic_history[nrow(diagnostic_history), , drop = FALSE]
    } else {
      NULL
    },
    sandwich_level = sandwich_level,
    md_correction = md_correction,
    corstr = corstr,
    n_iter = n_iter,
    tol = tol,
    tol_update = tol_update,
    tol_consensus = tol_consensus,
    tol_score = tol_score,
    step_size = step_size,
    ridge = ridge,
    L_beta = L_beta,
    L_S = L_S,
    L_B = L_B,
    W = W
  )
  
  class(out) <- "DecentralizedFedGEE"
  out
}

###############################################################################
# 5. PRINT METHOD
###############################################################################

print.DecentralizedFedGEE <- function(x, ...) {
  cat("Decentralized Federated GEE\n")
  cat("  Working correlation:", x$corstr, "\n")
  cat("  Sandwich level:", x$sandwich_level, "\n")
  cat("  Mancl-DeRouen:", x$md_correction, "\n")
  if (!is.null(x$L_init)) {
    cat("  Initialization: consensus meta-GLM | L_init =", x$L_init, "\n")
  }
  cat("  Consensus rounds: L_beta =", x$L_beta,
      "| L_S =", x$L_S, "| L_B =", x$L_B, "\n")
  cat("  Sites:", x$n_sites, "| Patients:", x$n_patients, "\n")
  cat("  Max iterations:", x$n_iter, "| Converged:", x$converged,
      "in", x$iterations, "iterations\n")
  
  if (!is.null(x$final_diagnostics)) {
    d <- x$final_diagnostics
    cat(sprintf(
      "  Final diagnostics: update %.3e | consensus %.3e | score %.3e\n",
      d$max_update_error, d$max_consensus_error, d$max_score_error
    ))
  }
  
  cat("\n  Decentralized estimator (network average, beta_dec):\n")
  print(round(x$coefficients, 4))
  
  cat("\n  Standard errors (from the averaged decentralized sandwich covariance):\n")
  print(round(x$se, 4))
  
  max_disagreement <- max(sqrt(rowSums(
    (x$coefficients_by_site - matrix(x$coefficients,
                                     nrow = nrow(x$coefficients_by_site),
                                     ncol = ncol(x$coefficients_by_site),
                                     byrow = TRUE))^2
  )))
  cat(sprintf(
    "\n  Max residual cross-site disagreement ||beta_i - beta_dec||: %.3e\n",
    max_disagreement
  ))
  cat("  (Per-site rows are available in $coefficients_by_site and $se_by_site\n",
      "   as diagnostics; beta_dec above is the reported point estimate.)\n", sep = "")
  
  invisible(x)
}

###############################################################################
# 6. CONVENIENCE WRAPPER
###############################################################################

decentralized_fedgee <- function(data_list,
                                 main_formula,
                                 family_obj = binomial(link = "logit"),
                                 corstr = "independence",
                                 id_col = "pat_id",
                                 W = NULL,
                                 structure = NULL,
                                 hub = 1,
                                 K_neigh = 2,
                                 region_id = NULL,
                                 hub_sites = NULL,
                                 site_constant_cols = character(0),
                                 L_init = NULL,
                                 L_beta = 1L,
                                 L_S = 1L,
                                 L_B = 1L,
                                 sandwich_level = "patient",
                                 md_correction = FALSE,
                                 n_iter = 50L,
                                 tol = 1e-8,
                                 tol_update = tol,
                                 tol_consensus = tol,
                                 tol_score = tol,
                                 step_size = 0.5,
                                 ridge = 1e-5,
                                 verbose = TRUE) {
  W <- get_weight_matrix(
    W = W,
    structure = structure,
    data_list = data_list,
    hub = hub,
    K_neigh = K_neigh,
    region_id = region_id,
    hub_sites = hub_sites,
    verbose = verbose
  )
  
  if (is.null(L_init)) L_init <- L_beta
  
  prep <- Prep_DecFedGEE(
    data_list = data_list,
    main_formula = main_formula,
    family_obj = family_obj,
    corstr = corstr,
    id_col = id_col,
    W = W,
    L_init = L_init,
    site_constant_cols = site_constant_cols,
    verbose = verbose
  )
  
  fit <- train_DecentralizedFedGEE(
    data_list = data_list,
    initial_beta_by_site = prep$initial_values_by_site,
    main_formula = main_formula,
    alpha = prep$alpha_list,
    family_obj = family_obj,
    corstr = corstr,
    id_col = id_col,
    W = W,
    L_beta = L_beta,
    L_S = L_S,
    L_B = L_B,
    n_iter = n_iter,
    tol = tol,
    tol_update = tol_update,
    tol_consensus = tol_consensus,
    tol_score = tol_score,
    step_size = step_size,
    sandwich_level = sandwich_level,
    md_correction = md_correction,
    ridge = ridge,
    verbose = verbose
  )
  
  if (is.null(fit)) return(NULL)
  
  # Retain initialization diagnostics. These fields do not alter the estimator
  # returned by the decentralized training routine.
  fit$L_init <- prep$L_init
  fit$initial_beta_by_site <- prep$initial_values_by_site
  fit$initial_site_status <- prep$site_status
  
  fit
}
