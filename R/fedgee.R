###############################################################################
# Fed-GEE: Federated Generalized Estimating Equations (centralized server)
#
# ONE fit computes every small-sample variance at once. Point estimates are
# identical across corrections; only the sandwich meat and the reference
# distribution change, so `correction` and `df` merely select which variant
# print()/summary()/vcov()/confint() report. Any other variant can be shown
# later without refitting, e.g. summary(fit, correction = "MD", df = "K-1").
#
# Site-level corrections (score space; sites transmit only p x p summaries):
#   "none" -- uncorrected sandwich, one outer product per site
#   "KC"   -- Kauermann-Carroll, (I - G_i)^{-1/2} via the similarity transform
#   "MD"   -- Mancl-DeRouen,     (I - G_i)^{-1}
#   "FG"   -- Fay-Graubard (2001), diagonal site-level leverage correction
#
# Reference distributions (site level):
#   "bm"  -- Bell-McCaffrey df in score space (design-only, federated)
#   "K-1", "K-p" -- fixed cluster counts;  "z" -- normal
# Patient-level sandwich (sandwich_level = "patient") supports none/KC/MD with
# df = N - p (or z).
###############################################################################

utils::globalVariables(".id_var")

.CORRECTIONS <- c("none", "KC", "MD", "FG")
.DF_METHODS <- c("bm", "K-1", "K-p", "z")
.FG_BOUND <- 0.75

###############################################################################
# 1. PREPARATION: site-constant columns, reduced formula, local alpha
#
# Covariates constant within EVERY site (e.g. a hospital-level indicator)
# cannot be estimated locally; they are dropped from the local formula used
# for alpha and recovered by the server from between-site variation. The
# coordinator sees every site's design, so this is detected automatically.
# Newton-Raphson starts from zero: meta-analytic starts from local GLMs are
# unstable at high p (sparse local binary covariates) and zero converges in
# a handful of iterations across all regimes tested.
###############################################################################
.prep_fedgee <- function(data_list, main_formula, family_obj, corstr, id_col) {
  X1 <- model.matrix(main_formula, data = data_list[[1]])
  colnames_full <- colnames(X1)

  constant <- Reduce(`&`, lapply(data_list, function(d) {
    X <- model.matrix(main_formula, data = d)
    if (!identical(colnames(X), colnames_full)) {
      stop(
        "Design matrices differ across sites. Make factor levels identical ",
        "at every site (e.g. set factor(x, levels = ...) before splitting)."
      )
    }
    apply(X, 2, function(x) all(abs(x - x[1]) < 1e-12))
  }))
  constant["(Intercept)" == colnames_full] <- FALSE
  constant_cols <- colnames_full[constant]

  reduced <- .reduced_formula(main_formula, X1, constant_cols, data_list[[1]])
  alpha <- lapply(data_list, .local_alpha,
    reduced_formula = reduced,
    family_obj = family_obj, id_col = id_col, corstr = corstr
  )
  .warn_clamped(alpha)

  list(
    initial_beta = matrix(0, length(colnames_full), 1,
      dimnames = list(colnames_full, NULL)
    ),
    alpha = alpha,
    constant_cols = constant_cols,
    reduced_formula = reduced
  )
}

###############################################################################
# 2. BELL-McCAFFREY DF IN SCORE SPACE
#
# Working model: S_i(b) = U_i - B_i (b - b0), U_i independent with
# Var(U_i) = tau B_i. For a variant with raw-score operator T_i, the variance of
# coefficient j is sum_i (a_i' S_i)^2 with a_i = T_i' B^{-1} e_j. Under the
# working model this is tau times a quadratic form in N(0, I) draws with K x K
# kernel
#     C_ik = 1{i = k} a_i' B_i a_i - (B_i a_i)' B^{-1} (B_k a_k),
# so the Satterthwaite df is tr(C)^2 / tr(C^2). Design-only: needs {B_i}.
# Balanced sites give exactly K - 1; matches clubSandwich CR2 at W = I.
###############################################################################
.bm_df <- function(T_list, B_sites, B_inv) {
  K <- length(B_sites)
  p <- nrow(B_inv)
  vapply(seq_len(p), function(j) {
    a <- lapply(T_list, function(Ti) crossprod(Ti, B_inv[, j]))
    Ba <- vapply(seq_len(K), function(i) drop(B_sites[[i]] %*% a[[i]]), numeric(p))
    Ba <- matrix(Ba, nrow = p)
    d <- vapply(seq_len(K), function(i) sum(a[[i]] * Ba[, i]), numeric(1))
    C <- diag(d, nrow = K) - crossprod(Ba, B_inv %*% Ba)
    trC2 <- sum(C * C)
    if (!is.finite(trC2) || trC2 < 1e-300) NA_real_ else sum(diag(C))^2 / trC2
  }, numeric(1))
}

# Fay-Graubard's own df (dtildeH) via saws. About 1 s per coefficient at
# K = 50, p = 20, so it is opt-in.
.fg_own_df <- function(beta_hat, S_sites, B_sites) {
  K <- length(S_sites)
  p <- length(beta_hat)
  u <- do.call(rbind, lapply(S_sites, as.numeric))
  om <- array(0, c(K, p, p))
  for (i in seq_len(K)) om[i, , ] <- B_sites[[i]]
  vapply(seq_len(p), function(j) {
    tm <- matrix(0, 1, p)
    tm[1, j] <- 1
    out <- tryCatch(
      suppressWarnings(saws::saws(
        list(coefficients = beta_hat, u = u, omega = om),
        test = tm, method = "d5", bound = .FG_BOUND
      )),
      error = function(e) NULL
    )
    if (is.null(out)) NA_real_ else as.numeric(out$df)
  }, numeric(1))
}

###############################################################################
# 3. ALL VARIANCE VARIANTS AT THE CONVERGED ESTIMATE
###############################################################################
.fedgee_variances <- function(sites, beta_hat, par_names, sandwich_level,
                              fg_df = FALSE) {
  K <- length(sites)
  p <- length(beta_hat)
  B_sites <- lapply(sites, `[[`, "B_site")
  S_sites <- lapply(sites, `[[`, "S_site")
  B <- Reduce(`+`, B_sites)
  B_inv <- .safe_solve(B)
  if (is.null(B_inv)) stop("The aggregated bread matrix is singular.")
  B_half <- .mat_pow_safe(B, 0.5, tol = 1e-12)$mat
  B_neghalf <- .mat_pow_safe(B, -0.5, tol = 1e-12)$mat
  n_patients <- sum(vapply(sites, `[[`, numeric(1), "n_clusters"))

  sandwich <- function(scores) {
    V <- B_inv %*% Reduce(`+`, lapply(scores, tcrossprod)) %*% B_inv
    V <- (V + t(V)) / 2
    dimnames(V) <- list(par_names, par_names)
    V
  }
  variant <- function(V, df_bm = rep(NA_real_, p), df_fg = rep(NA_real_, p)) {
    list(
      vcov = V, se = setNames(sqrt(pmax(diag(V), 0)), par_names),
      df_bm = setNames(df_bm, par_names), df_fg = setNames(df_fg, par_names)
    )
  }

  # --- site level: operators T_i on the raw site score ----------------------
  lev <- function(power) {
    lapply(B_sites, .leverage_op,
      B_half = B_half, B_neghalf = B_neghalf,
      power = power, tol = 0.05
    )
  }
  kc_ops <- lev(0.5)
  md_ops <- lev(1)
  q_diag <- lapply(B_sites, function(Bi) diag(Bi %*% B_inv))
  T_list <- list(
    none = rep(list(diag(p)), K),
    KC = lapply(kc_ops, `[[`, "mat"),
    MD = lapply(md_ops, `[[`, "mat"),
    FG = lapply(q_diag, function(q) diag((1 - pmin(.FG_BOUND, q))^(-0.5), nrow = p))
  )
  site <- lapply(names(T_list), function(nm) {
    Ts <- T_list[[nm]]
    scores <- Map(function(Ti, Si) Ti %*% Si, Ts, S_sites)
    variant(
      sandwich(scores),
      df_bm = .bm_df(Ts, B_sites, B_inv),
      df_fg = if (nm == "FG" && fg_df) .fg_own_df(beta_hat, S_sites, B_sites) else rep(NA_real_, p)
    )
  })
  names(site) <- names(T_list)

  # --- patient level (only when requested; needs per-cluster pieces) --------
  patient <- NULL
  if (sandwich_level == "patient") {
    B_cl <- unlist(lapply(sites, `[[`, "B_clusters"), recursive = FALSE)
    S_cl <- unlist(lapply(sites, `[[`, "S_clusters"), recursive = FALSE)
    patient <- lapply(c(none = 0, KC = 0.5, MD = 1), function(power) {
      scores <- if (power == 0) {
        S_cl
      } else {
        Map(function(Bj, Sj) {
          .leverage_op(Bj, B_half, B_neghalf, power = power, tol = 0.05)$mat %*% Sj
        }, B_cl, S_cl)
      }
      variant(sandwich(scores))
    })
  }

  # --- diagnostics ----------------------------------------------------------
  n_c <- vapply(sites, `[[`, numeric(1), "n_clusters")
  M_site <- Reduce(`+`, lapply(S_sites, tcrossprod))
  diagnostics <- list(
    site_size_cv = if (K > 1) stats::sd(n_c) / mean(n_c) else NA_real_,
    max_leverage = max(unlist(q_diag)),
    meat_rank = qr(M_site, tol = 1e-10)$rank,
    kc_clipped = sum(vapply(kc_ops, `[[`, logical(1), "clipped")),
    md_clipped = sum(vapply(md_ops, `[[`, logical(1), "clipped")),
    fg_capped = sum(vapply(q_diag, function(q) any(q > .FG_BOUND), logical(1))),
    kc_min_eig = min(vapply(kc_ops, `[[`, numeric(1), "min_eig"), na.rm = TRUE)
  )

  list(
    site = site, patient = patient, B = B, n_patients = n_patients,
    diagnostics = diagnostics
  )
}

###############################################################################
# 4. SELECT ONE VARIANT + REFERENCE DISTRIBUTION
###############################################################################
.fedgee_select <- function(object, correction, df) {
  correction <- match.arg(correction, .CORRECTIONS)
  df <- match.arg(df, .DF_METHODS)
  p <- length(object$coefficients)
  K <- object$n_sites

  if (object$sandwich_level == "patient") {
    if (correction == "FG") {
      stop("The FG correction is available only with sandwich_level = \"site\".")
    }
    v <- object$variances$patient[[correction]]
    dfv <- if (df == "z") Inf else object$n_patients - p
    label <- if (df == "z") "z" else "N-p"
  } else {
    v <- object$variances$site[[correction]]
    dfv <- switch(df,
      bm = v$df_bm,
      `K-1` = K - 1,
      `K-p` = K - p,
      z = Inf
    )
    label <- df
  }
  dfv <- rep_len(dfv, p)
  dfv[!is.finite(dfv) | dfv <= 0] <- Inf
  list(
    vcov = v$vcov, se = v$se, df = setNames(dfv, names(object$coefficients)),
    correction = correction, df_method = label
  )
}

###############################################################################
# 5. MAIN ENTRY POINT
###############################################################################

#' Federated Generalized Estimating Equations (FedGEE)
#'
#' Fits a GEE across sites that share only \eqn{p \times p} summaries
#' (bread and score) with a central server. The point estimate equals pooled
#' GEE when all sites use a common working correlation (always true for
#' \code{corstr = "independence"}). One fit computes every small-sample
#' variance; \code{correction} and \code{df} only choose which one is
#' reported, and \code{\link{summary.FedGEE}} can show any other without
#' refitting.
#'
#' @section Small-sample corrections (site level):
#' \describe{
#'   \item{\code{"none"}}{Uncorrected sandwich, one outer product per site.}
#'   \item{\code{"KC"}}{Kauermann--Carroll: site score multiplied by
#'     \eqn{B^{1/2}(I - G_i)^{-1/2}B^{-1/2}}, \eqn{G_i = B^{-1/2}B_iB^{-1/2}}.
#'     Unbiased under the working model; the recommended default.}
#'   \item{\code{"MD"}}{Mancl--DeRouen: power \eqn{-1}; conservative.}
#'   \item{\code{"FG"}}{Fay--Graubard: diagonal leverage correction with
#'     bound 0.75.}
#' }
#' @section Degrees of freedom:
#' \code{"bm"} is the Bell--McCaffrey (Satterthwaite) df computed in score
#' space from the site breads only; it adapts to unequal site sizes and is
#' the recommended pairing with KC. \code{"K-1"} and \code{"K-p"} are fixed
#' counts; \code{"z"} uses the normal. The patient-level sandwich always uses
#' \eqn{N - p} (or \code{"z"}).
#'
#' @section Rank limit:
#' The site-level meat is a sum of \eqn{K} outer products of scores that sum
#' to zero, so its rank is \eqn{\min(p, K - 1)}. With \eqn{K \le p} the
#' sandwich is singular and a warning is issued.
#'
#' @param data_list A list of data frames, one per site. Factor columns must
#'   have the same levels at every site.
#' @param main_formula Model formula.
#' @param family_obj A family object, e.g. \code{binomial()}.
#' @param corstr Working correlation: \code{"independence"},
#'   \code{"exchangeable"}, \code{"ar1"} or \code{"unstructured"}. The
#'   correlation parameter is estimated locally at each site.
#' @param id_col Name of the within-site cluster (patient) ID column.
#' @param sandwich_level \code{"site"} (default; one cluster per site) or
#'   \code{"patient"} (patients independent within and across sites).
#' @param correction Reported correction: \code{"KC"} (default),
#'   \code{"MD"}, \code{"FG"} or \code{"none"}.
#' @param df Reported reference distribution: \code{"bm"} (default),
#'   \code{"K-1"}, \code{"K-p"} or \code{"z"}.
#' @param fg_df Logical; also compute Fay--Graubard's own df via the
#'   \pkg{saws} package (slow for large \eqn{p}). Stored in
#'   \code{$variances$site$FG$df_fg}.
#' @param n_iter Maximum Newton--Raphson iterations.
#' @param tol Convergence tolerance on the update norm.
#' @param verbose Logical; print iteration progress.
#' @return An object of class \code{FedGEE}: \code{coefficients},
#'   \code{vcov}, \code{se}, \code{df} (per coefficient) for the selected
#'   variant; \code{variances} holding every variant; \code{diagnostics}
#'   (site-size CV, maximum leverage, meat rank, clipping counts); and fit
#'   metadata.
#' @references
#' Kauermann G, Carroll RJ (2001). JASA 96:1387--1396.
#' Mancl LA, DeRouen TA (2001). Biometrics 57:126--134.
#' Fay MP, Graubard BI (2001). Biometrics 57:1198--1206.
#' Bell RM, McCaffrey DF (2002). Survey Methodology 28:169--181.
#' @examples
#' data(ChickWeight)
#' cw <- ChickWeight
#' cw$site <- as.integer(cw$Chick) %% 12
#' fit <- fedgee(split(cw, cw$site), weight ~ Time + Diet,
#'   family_obj = gaussian(), id_col = "Chick", verbose = FALSE
#' )
#' fit
#' summary(fit, correction = "MD", df = "K-1")
#' @importFrom stats as.formula coef gaussian binomial glm pnorm pt qnorm qt var vcov
#'   model.matrix terms reformulate setNames confint delete.response
#' @importFrom utils globalVariables
#' @export
fedgee <- function(data_list,
                   main_formula,
                   family_obj = binomial(link = "logit"),
                   corstr = "independence",
                   id_col = "pat_id",
                   sandwich_level = c("site", "patient"),
                   correction = c("KC", "MD", "FG", "none"),
                   df = c("bm", "K-1", "K-p", "z"),
                   fg_df = FALSE,
                   n_iter = 50,
                   tol = 1e-8,
                   verbose = TRUE) {
  call <- match.call()
  sandwich_level <- match.arg(sandwich_level)
  correction <- match.arg(correction)
  df <- match.arg(df)
  corstr <- match.arg(corstr, c("independence", "exchangeable", "ar1", "unstructured"))
  if (is.function(family_obj)) family_obj <- family_obj()
  if (!is.list(data_list) || length(data_list) < 2L) {
    stop("data_list must be a list of at least two site data frames.")
  }
  if (fg_df && !requireNamespace("saws", quietly = TRUE)) {
    stop("fg_df = TRUE needs the 'saws' package.")
  }
  for (d in data_list) {
    if (!id_col %in% names(d)) stop("id_col '", id_col, "' is missing from a site.")
  }
  if (is.null(names(data_list))) names(data_list) <- paste0("site", seq_along(data_list))

  prep <- .prep_fedgee(data_list, main_formula, family_obj, corstr, id_col)
  beta <- prep$initial_beta
  par_names <- rownames(beta)
  K <- length(data_list)
  p <- length(beta)

  if (verbose) {
    cat(sprintf(
      "FedGEE | %d sites | %s | corstr = %s\n", K,
      family_obj$family, corstr
    ))
  }

  site_stats <- function(b, keep_clusters = FALSE) {
    Map(function(d, a) {
      .cluster_stats(d, b, a, main_formula, family_obj, id_col, corstr,
        keep_clusters = keep_clusters
      )
    }, data_list, prep$alpha)
  }

  converged <- FALSE
  for (k in seq_len(n_iter)) {
    sites <- Filter(Negate(is.null), site_stats(beta))
    if (length(sites) == 0) stop("No site returned usable statistics.")
    upd <- .safe_solve(
      Reduce(`+`, lapply(sites, `[[`, "B_site")),
      Reduce(`+`, lapply(sites, `[[`, "S_site"))
    )
    if (is.null(upd)) stop("Aggregated bread is singular at iteration ", k, ".")
    beta <- beta + upd
    step <- sqrt(sum(upd^2))
    if (verbose) cat(sprintf("  iter %2d | step %.2e\n", k, step))
    if (!is.finite(step) || step > 1e4) stop("Newton-Raphson diverged at iteration ", k, ".")
    if (step < tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) warning("FedGEE did not converge in ", n_iter, " iterations.")

  # Final pass at the converged estimate: all variance variants.
  sites <- Filter(Negate(is.null), site_stats(beta, keep_clusters = sandwich_level == "patient"))
  K_used <- length(sites)
  beta_hat <- setNames(as.vector(beta), par_names)
  var <- .fedgee_variances(sites, beta_hat, par_names, sandwich_level, fg_df = fg_df)

  if (sandwich_level == "site" && var$diagnostics$meat_rank < p) {
    warning(sprintf(
      paste0(
        "Site-level meat has rank %d < p = %d (rank <= min(p, K - 1) with K = %d). ",
        "The sandwich is singular; reduce p or use more sites."
      ),
      var$diagnostics$meat_rank, p, K_used
    ))
  }

  out <- list(
    coefficients = beta_hat,
    variances = list(site = var$site, patient = var$patient),
    diagnostics = var$diagnostics,
    B_global = var$B,
    alpha = prep$alpha,
    constant_cols = prep$constant_cols,
    n_sites = K_used,
    n_patients = var$n_patients,
    df_residual = if (sandwich_level == "site") K_used - p else var$n_patients - p,
    sandwich_level = sandwich_level,
    corstr = corstr,
    family = family_obj,
    formula = main_formula,
    iterations = k,
    converged = converged,
    call = call
  )
  sel <- .fedgee_select(out, correction, df)
  out$vcov <- sel$vcov
  out$se <- sel$se
  out$df <- sel$df
  out$correction <- sel$correction
  out$df_method <- sel$df_method
  class(out) <- "FedGEE"
  out
}

###############################################################################
# 6. METHODS
###############################################################################

#' Summarize a FedGEE fit
#'
#' Builds the coefficient table for any stored variance variant. Because
#' \code{\link{fedgee}} computes every correction in one fit, switching
#' \code{correction} or \code{df} here needs no refit.
#'
#' @param object A \code{FedGEE} fit.
#' @param correction One of \code{"KC"}, \code{"MD"}, \code{"FG"},
#'   \code{"none"}. Defaults to the fit's choice.
#' @param df One of \code{"bm"}, \code{"K-1"}, \code{"K-p"}, \code{"z"}.
#'   Defaults to the fit's choice.
#' @param level Confidence level.
#' @param x A \code{summary.FedGEE} object.
#' @param ... Unused.
#' @return An object of class \code{summary.FedGEE}.
#' @export
summary.FedGEE <- function(object, correction = object$correction,
                           df = NULL, level = 0.95, ...) {
  if (is.null(df)) df <- if (object$df_method %in% .DF_METHODS) object$df_method else "bm"
  sel <- .fedgee_select(object, correction, df)
  est <- object$coefficients
  stat <- est / sel$se
  mult <- qt(1 - (1 - level) / 2, df = sel$df)
  tab <- data.frame(
    Estimate = est, SE = sel$se, df = sel$df,
    lower = est - mult * sel$se, upper = est + mult * sel$se,
    statistic = stat, p.value = 2 * pt(-abs(stat), df = sel$df),
    check.names = FALSE
  )
  names(tab)[4:5] <- sprintf("%g%%", 100 * c((1 - level) / 2, 1 - (1 - level) / 2))
  structure(
    list(
      coefficients = tab, correction = sel$correction, df_method = sel$df_method,
      object = object
    ),
    class = "summary.FedGEE"
  )
}

#' @rdname summary.FedGEE
#' @export
print.summary.FedGEE <- function(x, ...) {
  o <- x$object
  cat("Federated GEE\n")
  cat("  Family / link       :", o$family$family, "/", o$family$link, "\n")
  cat("  Working correlation :", o$corstr, "\n")
  cat("  Sites               :", o$n_sites, "| Patients:", o$n_patients, "\n")
  cat("  Sandwich level      :", o$sandwich_level, "\n")
  cat("  SS correction       :", x$correction, "\n")
  cat("  Reference dist.     :", switch(x$df_method,
    bm = "t, Bell-McCaffrey df",
    z = "normal",
    paste0("t(", x$df_method, ")")
  ), "\n")
  cat("  Converged           :", o$converged, "in", o$iterations, "iterations\n")
  d <- o$diagnostics
  cat(sprintf(
    "  Site-size CV %.2f | max leverage %.2f | meat rank %d\n\n",
    d$site_size_cv, d$max_leverage, d$meat_rank
  ))
  tab <- x$coefficients
  tab[] <- lapply(tab, function(col) signif(col, 4))
  print(tab)
  invisible(x)
}

#' Print Method for Federated GEE
#'
#' @param x An object of class \code{FedGEE}.
#' @param ... Passed to \code{\link{summary.FedGEE}}.
#' @export
print.FedGEE <- function(x, ...) {
  print(summary(x, ...))
  invisible(x)
}

#' @export
coef.FedGEE <- function(object, ...) object$coefficients

#' Variance-covariance matrix of a FedGEE fit
#'
#' @param object A \code{FedGEE} fit.
#' @param correction Which stored variant; defaults to the fit's choice.
#' @param ... Unused.
#' @export
vcov.FedGEE <- function(object, correction = object$correction, ...) {
  .fedgee_select(object, correction, "z")$vcov
}

#' Confidence intervals for a FedGEE fit
#'
#' @param object A \code{FedGEE} fit.
#' @param parm Coefficient names or indices; default all.
#' @param level Confidence level.
#' @param correction,df As in \code{\link{summary.FedGEE}}.
#' @param ... Unused.
#' @export
confint.FedGEE <- function(object, parm, level = 0.95,
                           correction = object$correction, df = NULL, ...) {
  tab <- summary(object, correction = correction, df = df, level = level)$coefficients
  ci <- as.matrix(tab[, 4:5])
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}
