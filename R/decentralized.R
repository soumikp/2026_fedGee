###############################################################################
# Decentralized (serverless) Fed-GEE
#
# No server. Sites exchange summaries only with graph neighbours through a
# doubly stochastic gossip matrix W; L rounds of consensus shrink disagreement
# at rate rho^L, rho = second-largest |eigenvalue| of W.
#
# Estimation:
#   1. Warm start: each site fits a local GLM on the locally estimable
#      formula, forms P_i = V_i^{-1} and q_i = P_i beta_i, and runs L_init
#      consensus rounds on both; beta_i^(0) = P_tilde_i^{-1} q_tilde_i.
#   2. Each iteration: L_beta rounds on beta, local bread and score at the
#      consensus beta, L_S rounds on score and L_B on bread, one damped
#      Fisher-scoring step.
# Inference:
#   Local bread and meat at the final beta, L_B / L_S consensus rounds, and a
#   sandwich built from the consensus averages. Site i uses its own estimate
#   of the global bread, K * Bbar_i, so the leverage identity sum_i G_i = I
#   holds only to O(rho^L). Under-provisioned L makes the KC/MD correction
#   over-correct (conservative), not under-correct.
#
# Covariates constant within every site must be named in site_constant_cols:
# detecting them automatically would need all-sites visibility, which the
# decentralized protocol is designed to avoid.
###############################################################################

###############################################################################
# 0. GOSSIP WEIGHT MATRICES
###############################################################################

#' Metropolis-Hastings gossip weights
#'
#' Builds a symmetric, doubly stochastic weight matrix from an adjacency
#' matrix using the self-inclusive degree \eqn{d_i = } (neighbours + 1):
#' \eqn{w_{ij} = 1/\max(d_i, d_j)} on edges and \eqn{w_{ii} = 1 - \sum_{j\ne
#' i} w_{ij}}.
#'
#' @param A Symmetric 0/1 adjacency matrix with ones on the diagonal.
#' @return A \eqn{K \times K} weight matrix.
#' @export
mh_weights <- function(A) {
  A <- as.matrix(A)
  stopifnot(nrow(A) == ncol(A), isSymmetric(A))
  if (!all(diag(A) == 1)) stop("A must have self-loops: set diag(A) <- 1.")
  deg <- rowSums(A != 0)
  W <- (A != 0) / outer(deg, deg, pmax)
  diag(W) <- 0
  diag(W) <- 1 - rowSums(W)
  W
}

.adj_hub <- function(K, hub = 1) {
  A <- diag(1, K)
  A[hub, ] <- 1
  A[, hub] <- 1
  A
}

.adj_ring <- function(K, K_neigh = 2) {
  if (K_neigh %% 2 != 0 || K_neigh >= K) stop("K_neigh must be even and < K.")
  A <- diag(1, K)
  for (i in seq_len(K)) {
    for (d in seq_len(K_neigh / 2)) {
      A[i, ((i - 1 + d) %% K) + 1] <- 1
      A[i, ((i - 1 - d) %% K) + 1] <- 1
    }
  }
  A
}

# Sites in the same region are fully connected; regional hubs are fully
# connected to each other; non-hub sites in different regions share no edge.
.adj_visn <- function(region_id, hub_sites) {
  A <- diag(1, length(region_id))
  for (r in unique(region_id)) {
    idx <- which(region_id == r)
    A[idx, idx] <- 1
  }
  A[hub_sites, hub_sites] <- 1
  A
}

#' Build a gossip weight matrix for a network topology
#'
#' @param K Number of sites.
#' @param structure \code{"hub"} (star), \code{"ring"}, \code{"visn"}
#'   (regions fully connected internally, hubs connected to each other) or
#'   \code{"complete"}.
#' @param hub Hub index for \code{"hub"}.
#' @param K_neigh Even number of neighbours per site for \code{"ring"}.
#' @param region_id Length-\code{K} region labels for \code{"visn"}.
#' @param hub_sites Indices of each region's hub for \code{"visn"}.
#' @return A list with \code{W} (weights), \code{A} (adjacency),
#'   \code{structure} and \code{rho}, the second-largest absolute eigenvalue
#'   of \code{W}; consensus error decays like \eqn{\rho^L}.
#' @examples
#' build_weight_matrix(10, "ring")$rho
#' build_weight_matrix(10, "hub")$rho
#' @export
build_weight_matrix <- function(K,
                                structure = c("hub", "ring", "visn", "complete"),
                                hub = 1,
                                K_neigh = 2,
                                region_id = NULL,
                                hub_sites = NULL) {
  structure <- match.arg(structure)
  A <- switch(structure,
    hub = .adj_hub(K, hub),
    ring = .adj_ring(K, K_neigh),
    complete = matrix(1, K, K),
    visn = {
      if (is.null(region_id) || is.null(hub_sites)) {
        stop("structure = \"visn\" needs `region_id` (length K) and `hub_sites`.")
      }
      stopifnot(length(region_id) == K)
      .adj_visn(region_id, hub_sites)
    }
  )
  W <- mh_weights(A)
  ev <- sort(abs(eigen(W, symmetric = TRUE, only.values = TRUE)$values),
    decreasing = TRUE
  )
  list(W = W, A = A, structure = structure, rho = if (K > 1) ev[2] else 0)
}

# Validate W against data_list, reordering by site name when both are named.
.dec_resolve_W <- function(W, structure, data_list, hub, K_neigh, region_id,
                           hub_sites, tolerance = 1e-8) {
  K <- length(data_list)
  if (is.null(W)) {
    if (is.null(structure)) stop("Supply either `W` or `structure`.")
    W <- build_weight_matrix(K, structure, hub, K_neigh, region_id, hub_sites)$W
  }
  W <- as.matrix(W)
  site_names <- names(data_list)
  if (!all(dim(W) == c(K, K))) {
    stop(sprintf("W must be %d x %d to match data_list.", K, K))
  }
  if (!is.null(rownames(W)) && setequal(rownames(W), site_names) &&
    setequal(colnames(W), site_names)) {
    W <- W[site_names, site_names, drop = FALSE]
  }
  if (any(!is.finite(W)) || any(W < -tolerance)) {
    stop("W must be finite and nonnegative.")
  }
  if (max(abs(rowSums(W) - 1)) > tolerance) stop("Rows of W must sum to 1.")
  if (max(abs(colSums(W) - 1)) > 1e-6) {
    warning("W is not doubly stochastic; gossip will not converge to the plain network average.")
  }
  if (max(abs(W - t(W))) > 1e-6) {
    warning("W is not symmetric; the usual average-consensus theory does not apply directly.")
  }
  dimnames(W) <- list(site_names, site_names)
  W
}

# L rounds of neighbour averaging on a K x q object (one row per site).
.dec_consensus <- function(values, W, L) {
  if (length(L) != 1L || is.na(L) || L < 0 || L != as.integer(L)) {
    stop("Consensus rounds must be one nonnegative integer.")
  }
  out <- as.matrix(values)
  for (ell in seq_len(L)) out <- W %*% out
  dimnames(out) <- dimnames(as.matrix(values))
  out
}

.dec_rows <- function(mats) do.call(rbind, lapply(mats, as.vector))
.dec_unrows <- function(rows, p, nm) {
  lapply(seq_len(nrow(rows)), function(k) matrix(rows[k, ], p, p, dimnames = list(nm, nm)))
}

.dec_step <- function(step_size, t) {
  eta <- if (is.function(step_size)) {
    step_size(t)
  } else if (length(step_size) == 1L) {
    step_size
  } else {
    if (t > length(step_size)) stop("step_size vector is shorter than n_iter.")
    step_size[t]
  }
  if (length(eta) != 1L || !is.finite(eta) || eta <= 0) {
    stop("Each Fisher-scoring step size must be one positive finite number.")
  }
  as.numeric(eta)
}

###############################################################################
# 1. CONSENSUS META-GLM WARM START
###############################################################################
.dec_prep <- function(data_list, main_formula, family_obj, corstr, id_col, W,
                      L_init, site_constant_cols, verbose) {
  K <- length(data_list)
  site_names <- names(data_list)

  tt <- terms(main_formula, data = data_list[[1]])
  term_labels <- attr(tt, "term.labels")
  unknown <- setdiff(site_constant_cols, all.vars(delete.response(tt)))
  if (length(unknown)) {
    warning("site_constant_cols not in main_formula: ", paste(unknown, collapse = ", "))
  }
  uses_const <- vapply(term_labels, function(lab) {
    any(all.vars(stats::as.formula(paste("~", lab))) %in% site_constant_cols)
  }, logical(1))
  reduced <- reformulate(
    termlabels = if (any(!uses_const)) term_labels[!uses_const] else "1",
    response = main_formula[[2L]],
    intercept = attr(tt, "intercept") == 1L,
    env = environment(main_formula)
  )

  full_names <- colnames(model.matrix(main_formula, data = data_list[[1]]))
  red_names <- colnames(model.matrix(reduced, data = data_list[[1]]))
  pr <- length(red_names)
  constant_cols <- setdiff(full_names, red_names)
  if (length(constant_cols) && verbose) {
    message("Site-level coefficient(s) initialized at zero: ", paste(constant_cols, collapse = ", "))
  }

  P_list <- rep(list(matrix(0, pr, pr)), K)
  q_rows <- matrix(0, K, pr, dimnames = list(site_names, red_names))
  status <- setNames(rep("local_glm_excluded", K), site_names)

  for (i in seq_len(K)) {
    fit <- tryCatch(suppressWarnings(glm(reduced, data = data_list[[i]], family = family_obj)),
      error = function(e) NULL
    )
    if (is.null(fit) || !isTRUE(fit$converged)) next
    b <- coef(fit)
    if (!all(red_names %in% names(b)) || any(!is.finite(b[red_names]))) next
    V <- tryCatch(as.matrix(vcov(fit))[red_names, red_names], error = function(e) NULL)
    P <- if (is.null(V)) NULL else .safe_solve(V)
    if (is.null(P) || any(!is.finite(P))) next
    P_list[[i]] <- P
    q_rows[i, ] <- as.vector(P %*% b[red_names])
    status[i] <- "local_glm_used"
  }

  P_tilde <- .dec_unrows(.dec_consensus(.dec_rows(P_list), W, L_init), pr, red_names)
  q_tilde <- .dec_consensus(q_rows, W, L_init)

  init <- matrix(0, K, length(full_names), dimnames = list(site_names, full_names))
  for (i in seq_len(K)) {
    P <- P_tilde[[i]]
    if (!all(is.finite(P)) || max(abs(P)) <= sqrt(.Machine$double.eps)) {
      status[i] <- "fallback_zero"
      next
    }
    b0 <- .safe_solve(P, q_tilde[i, ])
    if (is.null(b0)) b0 <- .safe_solve(P, q_tilde[i, ], ridge = 1e-8 * max(1, abs(diag(P))))
    if (is.null(b0) || any(!is.finite(b0))) {
      status[i] <- "fallback_zero"
      next
    }
    init[i, red_names] <- as.vector(b0)
  }

  alpha <- lapply(data_list, .local_alpha,
    reduced_formula = reduced,
    family_obj = family_obj, id_col = id_col, corstr = corstr
  )
  .warn_clamped(alpha)

  list(
    initial_beta_by_site = init, alpha = alpha, constant_cols = constant_cols,
    reduced_formula = reduced, site_status = status
  )
}

###############################################################################
# 2. LOCAL BREAD / SCORE / MEAT AT ONE SITE
#
# With B_global (the site's own estimate K * Bbar_i) and power > 0, the site
# score (site level) or each patient score (patient level) is multiplied by
# the similarity-transform leverage operator before forming the meat.
###############################################################################
.dec_site_stats <- function(site_data, beta, alpha, main_formula, family_obj,
                            id_col, corstr, sandwich_level, power = 0,
                            B_global = NULL) {
  st <- .cluster_stats(site_data, beta, alpha, main_formula, family_obj,
    id_col, corstr,
    keep_clusters = sandwich_level == "patient"
  )
  if (is.null(st)) {
    return(NULL)
  }

  op <- function(Bu) diag(nrow(Bu))
  if (power > 0 && !is.null(B_global)) {
    B_half <- .mat_pow_safe(B_global, 0.5, tol = 1e-10)$mat
    B_neghalf <- .mat_pow_safe(B_global, -0.5, tol = 1e-10)$mat
    if (is.null(B_half) || is.null(B_neghalf)) {
      return(NULL)
    }
    op <- function(Bu) .leverage_op(Bu, B_half, B_neghalf, power, tol = 1e-6)$mat
  }

  meat <- if (sandwich_level == "site") {
    tcrossprod(op(st$B_site) %*% st$S_site)
  } else {
    Reduce(`+`, Map(function(Bj, Sj) tcrossprod(op(Bj) %*% Sj), st$B_clusters, st$S_clusters))
  }
  list(Bread = st$B_site, Score = st$S_site, Meat = meat, n_clusters = st$n_clusters)
}

###############################################################################
# 3. MAIN ENTRY POINT
###############################################################################

#' Decentralized (serverless) federated GEE
#'
#' Fits the Fed-GEE model with no central server. Sites exchange bread,
#' score and coefficient summaries only with neighbours through gossip
#' averaging on the weight matrix \code{W}. With a complete graph (or enough
#' consensus rounds) the fit reproduces \code{\link{fedgee}}.
#'
#' The number of rounds needed scales like
#' \eqn{L \gtrsim \log K / (2(1 - \rho))}, with \eqn{\rho} returned by
#' \code{\link{build_weight_matrix}}. Too few rounds harm the variance before
#' the point estimate: check \code{$se_by_site} and
#' \code{$final_diagnostics}. Convergence in \eqn{L} need not be monotone on
#' badly mixing graphs.
#'
#' @inheritParams fedgee
#' @param W A \eqn{K \times K} doubly stochastic gossip matrix whose rows and
#'   columns follow \code{data_list}. If \code{NULL}, built from
#'   \code{structure}.
#' @param structure,hub,K_neigh,region_id,hub_sites Topology passed to
#'   \code{\link{build_weight_matrix}} when \code{W} is \code{NULL}.
#' @param site_constant_cols Names of covariates constant within every site
#'   (initialized at zero, estimated from between-site variation).
#' @param L_init,L_beta,L_S,L_B Consensus rounds for the warm start,
#'   coefficients, scores/meat and breads. \code{L_init} defaults to
#'   \code{L_beta}.
#' @param sandwich_level \code{"patient"} (default) or \code{"site"}.
#' @param correction Leverage correction of the meat: \code{"none"}
#'   (default), \code{"KC"} or \code{"MD"}. Each site uses its own consensus
#'   estimate of the global bread.
#' @param md_correction Deprecated. In earlier code \code{TRUE} applied the
#'   KC (power 1/2) correction despite its name; it now maps to
#'   \code{correction = "KC"} with a warning.
#' @param step_size Fisher-scoring step: a number, a vector (one per
#'   iteration) or a function of the iteration.
#' @param tol_update,tol_consensus,tol_score Convergence tolerances on the
#'   update, neighbour disagreement and Newton direction.
#' @param ridge Ridge added before inverting consensus breads.
#' @return An object of class \code{DecentralizedFedGEE} with the network
#'   average \code{coefficients}, averaged \code{vcov} and \code{se},
#'   per-site rows (\code{coefficients_by_site}, \code{se_by_site}),
#'   iteration diagnostics and \code{df_residual} (\eqn{K - p} at site level,
#'   \eqn{N - p} at patient level).
#' @examples
#' data(ChickWeight)
#' cw <- ChickWeight
#' cw$site <- as.integer(cw$Chick) %% 8
#' dl <- split(cw, cw$site)
#' dfit <- decentralized_fedgee(dl, weight ~ Time,
#'   family_obj = gaussian(), id_col = "Chick", structure = "ring",
#'   L_beta = 20, L_S = 20, L_B = 20, verbose = FALSE
#' )
#' dfit
#' @export
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
                                 sandwich_level = c("patient", "site"),
                                 correction = c("none", "KC", "MD"),
                                 md_correction = NULL,
                                 n_iter = 50L,
                                 tol = 1e-8,
                                 tol_update = tol,
                                 tol_consensus = tol,
                                 tol_score = tol,
                                 step_size = 0.5,
                                 ridge = 1e-5,
                                 verbose = TRUE) {
  call <- match.call()
  sandwich_level <- match.arg(sandwich_level)
  correction <- match.arg(correction)
  if (!is.null(md_correction)) {
    warning(
      "`md_correction` is deprecated. It always applied the KC (power 1/2) ",
      "correction; use correction = \"KC\" or \"MD\" instead."
    )
    if (isTRUE(md_correction)) correction <- "KC"
  }
  corstr <- match.arg(corstr, c("independence", "exchangeable", "ar1", "unstructured"))
  if (is.function(family_obj)) family_obj <- family_obj()
  if (!is.list(data_list) || length(data_list) < 2L) {
    stop("data_list must be a list of at least two site data frames.")
  }
  if (is.null(names(data_list))) names(data_list) <- paste0("site", seq_along(data_list))
  if (anyDuplicated(names(data_list))) stop("Site names in data_list must be unique.")
  if (length(n_iter) != 1L || n_iter < 1L || n_iter != as.integer(n_iter)) {
    stop("n_iter must be one positive integer.")
  }
  tols <- c(tol_update, tol_consensus, tol_score)
  if (any(!is.finite(tols)) || any(tols < 0)) stop("Tolerances must be nonnegative and finite.")
  for (d in data_list) {
    if (!id_col %in% names(d)) stop("id_col '", id_col, "' is missing from a site.")
  }

  W <- .dec_resolve_W(W, structure, data_list, hub, K_neigh, region_id, hub_sites)
  if (is.null(L_init)) L_init <- L_beta
  for (L in list(L_init, L_beta, L_S, L_B)) .dec_consensus(matrix(0, nrow(W), 1), W, L)

  prep <- .dec_prep(
    data_list, main_formula, family_obj, corstr, id_col, W, L_init,
    site_constant_cols, verbose
  )

  K <- length(data_list)
  site_names <- names(data_list)
  beta <- prep$initial_beta_by_site
  par_names <- colnames(beta)
  p <- ncol(beta)

  if (verbose) {
    cat(sprintf(
      "Decentralized FedGEE | %d sites | %s | corstr = %s | L_beta/L_S/L_B = %d/%d/%d\n",
      K, family_obj$family, corstr, L_beta, L_S, L_B
    ))
  }

  local_stats <- function(beta_rows, power = 0, B_list = NULL) {
    lapply(seq_len(K), function(i) {
      .dec_site_stats(data_list[[i]], beta_rows[i, ], prep$alpha[[i]], main_formula,
        family_obj, id_col, corstr, sandwich_level,
        power = power, B_global = if (is.null(B_list)) NULL else B_list[[i]]
      )
    })
  }
  fail <- function(msg, idx) {
    if (verbose) cat("Decentralized FedGEE stopped:", msg, paste(idx, collapse = ", "), "\n")
    NULL
  }

  history <- list(beta)
  diag_hist <- vector("list", n_iter)
  converged <- FALSE
  iters <- 0L

  for (t in seq_len(n_iter)) {
    eta <- .dec_step(step_size, t)
    beta_tilde <- .dec_consensus(beta, W, L_beta)
    out <- local_stats(beta_tilde)
    bad <- which(vapply(out, is.null, logical(1)))
    if (length(bad)) {
      return(fail("local GEE statistics failed at site(s)", bad))
    }

    S_tilde <- .dec_consensus(do.call(rbind, lapply(out, function(o) as.vector(o$Score))), W, L_S)
    B_tilde <- .dec_unrows(.dec_consensus(.dec_rows(lapply(out, `[[`, "Bread")), W, L_B), p, par_names)

    direction <- matrix(NA_real_, K, p, dimnames = list(site_names, par_names))
    for (i in seq_len(K)) {
      d_i <- .safe_solve(B_tilde[[i]], S_tilde[i, ], ridge = ridge)
      if (!is.null(d_i)) direction[i, ] <- as.vector(d_i)
    }
    bad <- which(!stats::complete.cases(direction))
    if (length(bad)) {
      return(fail("consensus bread was singular at site(s)", bad))
    }
    beta_next <- beta_tilde + eta * direction
    dimnames(beta_next) <- list(site_names, par_names)

    upd <- max(sqrt(rowSums((beta_next - beta)^2)))
    cons <- max(sqrt(rowSums((beta_next - W %*% beta_next)^2)))
    score <- max(sqrt(rowSums(direction^2)))
    diag_hist[[t]] <- data.frame(
      iteration = t, step_size = eta, max_update_error = upd,
      max_consensus_error = cons, max_score_error = score
    )
    beta <- beta_next
    history[[t + 1L]] <- beta
    iters <- t
    if (verbose) {
      cat(sprintf("  iter %2d | update %.2e | consensus %.2e | score %.2e\n", t, upd, cons, score))
    }
    if (!is.finite(upd) || upd > 1e4 || score > 1e4) {
      return(fail("divergence at iteration", t))
    }
    if (upd < tol_update && cons < tol_consensus && score < tol_score) {
      converged <- TRUE
      break
    }
  }
  diag_hist <- do.call(rbind, diag_hist[seq_len(iters)])
  if (!converged) warning("Decentralized FedGEE did not converge in ", n_iter, " iterations.")

  # ---- Inference at the final consensus beta ----
  beta_hat <- .dec_consensus(beta, W, L_beta)
  final <- local_stats(beta_hat)
  bad <- which(vapply(final, is.null, logical(1)))
  if (length(bad)) {
    return(fail("inference failed at site(s)", bad))
  }

  B_bar <- .dec_unrows(.dec_consensus(.dec_rows(lapply(final, `[[`, "Bread")), W, L_B), p, par_names)
  power <- .correction_power(correction)
  if (power > 0) {
    # Each site's consensus Bbar_i approximates (1/K) sum_k B_k.
    B_glob <- lapply(B_bar, function(Bb) K * (Bb + t(Bb)) / 2)
    final_c <- local_stats(beta_hat, power = power, B_list = B_glob)
    bad <- which(vapply(final_c, function(o) is.null(o) || any(!is.finite(o$Meat)), logical(1)))
    if (length(bad)) {
      return(fail("leverage-corrected meat failed at site(s)", bad))
    }
    for (i in seq_len(K)) final[[i]]$Meat <- final_c[[i]]$Meat
  }
  M_bar <- .dec_unrows(.dec_consensus(.dec_rows(lapply(final, `[[`, "Meat")), W, L_S), p, par_names)

  # With average-consensus outputs, Var(beta_hat) = K^{-1} Bbar^{-1} Mbar Bbar^{-1}.
  vcov_by_site <- vector("list", K)
  se_by_site <- matrix(NA_real_, K, p, dimnames = list(site_names, par_names))
  for (i in seq_len(K)) {
    Bi <- .safe_solve(B_bar[[i]], ridge = ridge)
    if (is.null(Bi)) next
    V <- Bi %*% M_bar[[i]] %*% Bi / K
    V <- (V + t(V)) / 2
    dimnames(V) <- list(par_names, par_names)
    vcov_by_site[[i]] <- V
    se_by_site[i, ] <- sqrt(pmax(diag(V), 0))
  }
  names(vcov_by_site) <- site_names

  # Reported estimator: the network average beta_dec = (1/K) sum_i beta_i.
  # W doubly stochastic => unchanged by the final consensus step. vcov_dec is
  # the matching average of per-site sandwiches (one more consensus round).
  ok <- Filter(Negate(is.null), vcov_by_site)
  vcov_dec <- if (length(ok)) {
    Reduce(`+`, ok) / length(ok)
  } else {
    matrix(NA_real_, p, p, dimnames = list(par_names, par_names))
  }
  n_cl <- vapply(final, `[[`, numeric(1), "n_clusters")
  n_valid <- sum(n_cl > 0)

  structure(
    list(
      coefficients = setNames(colMeans(beta_hat), par_names),
      vcov = vcov_dec,
      se = setNames(sqrt(pmax(diag(vcov_dec), 0)), par_names),
      coefficients_by_site = beta_hat,
      se_by_site = se_by_site,
      vcov_by_site = vcov_by_site,
      Bread_by_site = B_bar,
      Meat_by_site = M_bar,
      df_residual = if (sandwich_level == "site") n_valid - p else sum(n_cl) - p,
      n_sites = K,
      n_valid_sites = n_valid,
      n_patients = sum(n_cl),
      iterations = iters,
      converged = converged,
      history = history,
      diagnostic_history = diag_hist,
      final_diagnostics = if (NROW(diag_hist)) diag_hist[nrow(diag_hist), , drop = FALSE] else NULL,
      initial_beta_by_site = prep$initial_beta_by_site,
      initial_site_status = prep$site_status,
      alpha = prep$alpha,
      constant_cols = prep$constant_cols,
      sandwich_level = sandwich_level,
      correction = correction,
      corstr = corstr,
      family = family_obj,
      formula = main_formula,
      W = W,
      L_init = as.integer(L_init), L_beta = L_beta, L_S = L_S, L_B = L_B,
      step_size = step_size, ridge = ridge, n_iter = n_iter,
      call = call
    ),
    class = "DecentralizedFedGEE"
  )
}

###############################################################################
# 4. METHODS
###############################################################################

#' Print a decentralized FedGEE fit
#'
#' @param x A \code{DecentralizedFedGEE} fit.
#' @param ... Unused.
#' @export
print.DecentralizedFedGEE <- function(x, ...) {
  cat("Decentralized Federated GEE\n")
  cat("  Family / link       :", x$family$family, "/", x$family$link, "\n")
  cat("  Working correlation :", x$corstr, "\n")
  cat("  Sites               :", x$n_sites, "| Patients:", x$n_patients, "\n")
  cat("  Sandwich level      :", x$sandwich_level, "\n")
  cat("  SS correction       :", x$correction, "\n")
  cat(
    "  Consensus rounds    : L_init", x$L_init, "| L_beta", x$L_beta,
    "| L_S", x$L_S, "| L_B", x$L_B, "\n"
  )
  cat("  Converged           :", x$converged, "in", x$iterations, "iterations\n")

  dfv <- if (x$df_residual > 0) x$df_residual else Inf
  cat("  Reference dist.     :", if (is.finite(dfv)) sprintf("t(%d)", dfv) else "normal", "\n\n")
  stat <- x$coefficients / x$se
  mult <- qt(0.975, dfv)
  tab <- data.frame(
    Estimate = x$coefficients, SE = x$se,
    `2.5%` = x$coefficients - mult * x$se, `97.5%` = x$coefficients + mult * x$se,
    statistic = stat, p.value = 2 * pt(-abs(stat), dfv), check.names = FALSE
  )
  tab[] <- lapply(tab, signif, 4)
  print(tab)

  gap <- max(sqrt(rowSums(sweep(x$coefficients_by_site, 2, x$coefficients)^2)))
  cat(sprintf("\n  Max site disagreement ||beta_i - beta_dec||: %.2e\n", gap))
  invisible(x)
}

#' @export
coef.DecentralizedFedGEE <- function(object, ...) object$coefficients

#' @export
vcov.DecentralizedFedGEE <- function(object, ...) object$vcov
