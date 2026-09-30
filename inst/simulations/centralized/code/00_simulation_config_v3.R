###############################################################################
# FedGEE Simulation Framework — v2
#
# Major changes from v1:
#   1. ONE fedgee_v2() call per replicate produces all variance variants.
#      Long-format output: one row per (estimator, variance_variant, inference).
#   2. New methods added: KC, MD, FG, VD (variance decomposition).
#   3. Bell-McCaffrey Satterthwaite df reported alongside K-p df.
#   4. GLMM fit is now controlled by a grid flag (`fit_glmm`) — only fit when
#      the estimand comparison is the goal (saves ~80% wall-clock for the
#      experiments that don't need it).
#   5. Per-replicate diagnostics stored: KC/MD eigenvalue minima, FG max obs
#      leverage, count of regularization clips, n_sites contributing to VD.
#   6. Parameter grid extended for new experiments (E1-E8). Grid flags:
#      `experiment` (E1-E8), `p_extra` (extra nuisance covariates),
#      `site_size_dist` (uniform/gamma/extreme), `corstr_truth`/`corstr_work`,
#      `fit_glmm`.
###############################################################################

library(geepack)
library(lme4)
library(dplyr)
library(purrr)
library(tibble)
library(Matrix)

# Source the v2 implementation
# Path-independent: the cluster driver sets FEDGEE_SIM_DIR, otherwise fall back
# to the working directory. A bare source("01_fedgee_v2.R") breaks whenever the
# job is launched from anywhere other than the code directory.
if (!exists("fedgee_v2")) {
  source(file.path(Sys.getenv("FEDGEE_SIM_DIR", unset = getwd()),
                   "01_fedgee_v3.R"))
}

###############################################################################
# A. DATA-GENERATING PROCESS (extended)
#
# New arguments:
#   p_extra: number of additional nuisance covariates to include
#            (mix of continuous, binary patient-level, and one site-level binary)
#   site_size_dist: "uniform" (current default), "gamma" (VHA-like), "extreme"
###############################################################################
generate_sim_data <- function(rho_H, rho_P,
                              n_hospitals,
                              patients_per_hosp_range = c(10, 30),
                              visits_per_patient_range = c(2, 6),
                              beta_0 = -1,
                              true_beta_x = 1,
                              include_teaching = FALSE,
                              beta_teach = 0.5,
                              p_extra = 0,
                              site_size_dist = "uniform") {

  stopifnot(rho_P >= rho_H, rho_H >= 0, rho_P < 1)

  U <- rnorm(n_hospitals)

  # Site-level teaching indicator
  teaching <- rep(0L, n_hospitals)
  if (include_teaching) {
    teaching[sample(seq_len(n_hospitals), floor(n_hospitals / 2))] <- 1L
  }

  # ---- Site sizes ----
  if (site_size_dist == "uniform") {
    n_pat_per_hosp <- sample(
      patients_per_hosp_range[1]:patients_per_hosp_range[2],
      n_hospitals, replace = TRUE)
  } else if (site_size_dist %in% c("cv078", "cv100")) {
    # v3: unequal site sizes with the SAME mean (~20 patients/site) as the
    # uniform arm, so imbalance is not confounded with total sample size.
    # Gamma shape tuned so the realised CV after the floor of 5 is ~0.78
    # (VHA-matched) or ~1.00 (checked on 1e6 draws, 2026-09-27).
    cv_par <- if (site_size_dist == "cv078") 0.805 else 1.062
    raw <- rgamma(n_hospitals, shape = 1 / cv_par^2, scale = 20 * cv_par^2)
    n_pat_per_hosp <- pmax(round(raw), 5)
  } else if (site_size_dist == "gamma") {
    # v2 legacy: mean ~100, sd ~100
    raw <- rgamma(n_hospitals, shape = 1, scale = 100)
    n_pat_per_hosp <- pmax(round(raw), 5)
  } else if (site_size_dist == "extreme") {
    # 10% large (1000 each), 90% small (20 each)
    n_large <- max(1, round(0.1 * n_hospitals))
    n_pat_per_hosp <- c(rep(1000, n_large), rep(20, n_hospitals - n_large))
    n_pat_per_hosp <- sample(n_pat_per_hosp)
  } else {
    stop("Unknown site_size_dist: ", site_size_dist)
  }

  # ---- Extra site-level covariate (one binary, exists if p_extra >= 1) ----
  has_site_extra <- p_extra >= 1
  if (has_site_extra) {
    site_extra <- rbinom(n_hospitals, 1, 0.5)
  }

  # ---- Patient table ----
  patients <- data.frame(
    hosp_id = rep(seq_len(n_hospitals), times = n_pat_per_hosp),
    pat_num = unlist(lapply(n_pat_per_hosp, seq_len)))
  patients$pat_id <- paste0(patients$hosp_id, "_", patients$pat_num)

  # ---- Generate visits ----
  visits_list <- lapply(seq_len(nrow(patients)), function(p) {
    h     <- patients$hosp_id[p]
    n_vis <- sample(visits_per_patient_range[1]:visits_per_patient_range[2], 1)
    x     <- rnorm(n_vis)

    eta <- beta_0 + true_beta_x * x
    if (include_teaching) eta <- eta + beta_teach * teaching[h]
    prob_marg <- plogis(eta)

    W   <- rnorm(1)
    eps <- rnorm(n_vis)
    Z   <- sqrt(rho_H) * U[h] +
           sqrt(rho_P - rho_H) * W +
           sqrt(1 - rho_P) * eps

    y <- as.integer(Z < qnorm(prob_marg))

    out <- data.frame(
      hosp_id = h,
      pat_id  = patients$pat_id[p],
      teaching = teaching[h],
      x = x, y = y)

    # Add nuisance covariates: NoiseN (continuous), NoiseB (binary patient)
    if (p_extra >= 1) out$site_extra <- site_extra[h]
    if (p_extra >= 2) {
      n_continuous <- ceiling((p_extra - 1) / 2)
      for (q in seq_len(n_continuous)) {
        out[[paste0("NoiseN", q)]] <- rnorm(n_vis)
      }
      n_binary <- (p_extra - 1) - n_continuous
      if (n_binary > 0) {
        for (q in seq_len(n_binary)) {
          out[[paste0("NoiseB", q)]] <- rbinom(n_vis, 1, 0.3)
        }
      }
    }
    out
  })

  bind_rows(visits_list) |> arrange(hosp_id, pat_id)
}

###############################################################################
# B. BUILD FORMULA based on grid row
###############################################################################
build_formula <- function(include_teaching, p_extra) {
  rhs_terms <- "x"
  if (include_teaching) rhs_terms <- c(rhs_terms, "teaching")
  if (p_extra >= 1) rhs_terms <- c(rhs_terms, "site_extra")
  if (p_extra >= 2) {
    n_continuous <- ceiling((p_extra - 1) / 2)
    n_binary <- (p_extra - 1) - n_continuous
    if (n_continuous > 0) {
      rhs_terms <- c(rhs_terms, paste0("NoiseN", seq_len(n_continuous)))
    }
    if (n_binary > 0) {
      rhs_terms <- c(rhs_terms, paste0("NoiseB", seq_len(n_binary)))
    }
  }
  as.formula(paste("y ~", paste(rhs_terms, collapse = " + ")))
}

###############################################################################
# C. META-ANALYTIC GEE  (per-site GEE + inverse-variance pool)
#
# Each site fits its OWN GEE (patients as clusters, robust/sandwich SE, working
# independence to match the pooled and federated fits), then the site
# coefficients are combined by fixed-effect inverse-variance weighting. This is
# the honest privacy-preserving comparator: only per-site (beta, se) leave a
# site, never patient rows.
#
# Small-site / rare-outcome failures (separation, non-convergence, non-finite
# robust SE) are handled by DROPPING that site from the pool and reporting the
# number of contributing sites via attribute "n_sites". A pooled estimate needs
# >= 2 usable sites.
###############################################################################
fit_meta_gee <- function(data, formula, SEP_THRESH = 10) {
  rhs_terms <- attr(terms(formula), "term.labels")
  resp      <- as.character(formula)[2]
  sites <- split(data, data$hosp_id)
  site_fits <- lapply(sites, function(sd) {
    tryCatch({
      # geeglm needs cluster rows contiguous; sort by patient within the site.
      sd <- sd[order(sd$pat_id), ]
      sd$pat_num_id <- as.numeric(factor(sd$pat_id))
      # Site-level covariates (site_extra, teaching) are CONSTANT within a
      # site, so they are structurally inestimable in a single-site fit --
      # geeglm errors ("rank deficient") rather than aliasing them. Build a
      # per-site formula from only the terms that actually vary here. The
      # patient/visit-level effect x always varies, so it is retained.
      varying <- rhs_terms[vapply(rhs_terms,
                    function(tm) length(unique(sd[[tm]])) > 1, logical(1))]
      if (length(varying) == 0) return(NULL)
      f_site <- as.formula(paste(resp, "~", paste(varying, collapse = " + ")))
      # v3 PRE-SCREEN (2026-09-27): geepack's compiled solver (update_beta) can
      # loop forever on a separated site -- found at K=18, p=5, rare outcome,
      # site with 1 event in 83 rows. The post-fit separation guard below never
      # runs because geeglm never returns. Screen BEFORE geeglm: drop the site
      # if it has < 2 events or non-events, or if a plain glm (which handles
      # separation safely) fails to converge, gives an extreme coefficient, or
      # pins fitted probabilities at 0/1. Same drop-and-count rule as below.
      yv <- sd[[resp]]
      if (min(sum(yv == 1), sum(yv == 0)) < 2) return(NULL)
      g0 <- tryCatch(suppressWarnings(glm(f_site, data = sd, family = binomial)),
                     error = function(e) NULL)
      if (is.null(g0) || !g0$converged) return(NULL)
      c0 <- coef(g0)
      if (any(!is.finite(c0)) || max(abs(c0)) > SEP_THRESH) return(NULL)
      fv <- fitted(g0)
      if (any(fv < 1e-8 | fv > 1 - 1e-8)) return(NULL)
      # Belt-and-suspenders: swallow any stray stdout from a residual
      # rank-deficiency path so 189k tasks don't bloat the logs.
      invisible(utils::capture.output(
        fit <- suppressWarnings(
          geeglm(f_site, id = pat_num_id, data = sd, family = binomial,
                 corstr = "independence"))))
      co <- coef(fit)
      vc <- vcov(fit)                       # geepack vcov = robust sandwich
      se <- suppressWarnings(sqrt(diag(vc)))  # separated sites -> NaN, dropped below
      # Separation guard: a site with (near-)complete separation returns a huge
      # log-OR with a deceptively small robust SE, which then dominates the
      # inverse-variance pool. Any |coef| beyond SEP_THRESH (OR > e^10 ~ 22000)
      # marks the whole site fit as unreliable -> drop it and count it. This is
      # what a careful federated analyst would do; the drop rate is reported so
      # the underlying instability stays visible.
      co_fin <- co[is.finite(co)]
      if (length(co_fin) == 0 || max(abs(co_fin)) > SEP_THRESH) return(NULL)
      keep <- is.finite(co) & is.finite(se) & se > 0
      if (!any(keep)) return(NULL)
      tibble(
        term = names(co)[keep],
        beta = unname(co[keep]),
        se   = unname(se[keep])
      )
    }, error = function(e) NULL)
  })
  n_attempted <- length(site_fits)
  site_fits <- Filter(Negate(is.null), site_fits)
  if (length(site_fits) < 2) return(NULL)

  site_fits <- bind_rows(site_fits)
  site_fits <- site_fits[is.finite(site_fits$se) & site_fits$se > 0 &
                           is.finite(site_fits$beta), ]
  if (nrow(site_fits) < 2) return(NULL)

  pooled <- site_fits |>
    group_by(term) |>
    summarise(
      w = list(1 / se^2), b = list(beta),
      n_sites = n(),
      .groups = "drop"
    ) |>
    mutate(
      beta = map2_dbl(b, w, ~ sum(.y * .x) / sum(.y)),
      se   = map_dbl(w, ~ 1 / sqrt(sum(.x)))
    ) |>
    select(term, beta, se, n_sites)

  attr(pooled, "n_attempted") <- n_attempted
  pooled
}

###############################################################################
# D. SINGLE REPLICATION FUNCTION
#
# Output: long-format data frame, one row per (method, inference) combination.
# Required columns for downstream aggregation:
#   estimator   -- "pooled_gee" | "fed" | "meta_glm" | "glmm"
#   variant     -- variance variant: "patient_z", "site_z", "site_t", "site_satt",
#                  "kc", "md", "fg", "vd", "naive" (for non-fed methods)
#   inference   -- "z" | "t_kp" | "t_satt"
#   coef        -- "x" | "teaching" | etc.
#   beta_hat, se, df, ci_lo, ci_hi, covers
###############################################################################
run_one_rep <- function(params, rep_id) {

  # ---- Unpack ----
  rho_H <- params$rho_H; rho_P <- params$rho_P; K <- params$K
  pat_range <- c(params$pat_lo, params$pat_hi)
  vis_range <- c(params$vis_lo, params$vis_hi)
  incl_teach <- params$include_teaching
  beta_teach <- params$beta_teach
  true_beta_x <- params$true_beta_x
  beta_0 <- if (is.null(params$beta_0)) -1 else params$beta_0
  p_extra <- if (is.null(params$p_extra)) 0 else params$p_extra
  site_size_dist <- if (is.null(params$site_size_dist)) "uniform" else params$site_size_dist
  fit_glmm_flag  <- if (is.null(params$fit_glmm)) FALSE else params$fit_glmm
  experiment_id  <- if (is.null(params$experiment)) "E0" else params$experiment

  form <- build_formula(incl_teach, p_extra)

  # ---- Generate data ----
  d <- generate_sim_data(
    rho_H = rho_H, rho_P = rho_P, n_hospitals = K,
    patients_per_hosp_range = pat_range,
    visits_per_patient_range = vis_range,
    beta_0 = beta_0, true_beta_x = true_beta_x,
    include_teaching = incl_teach, beta_teach = beta_teach,
    p_extra = p_extra, site_size_dist = site_size_dist
  )
  data_list <- split(d, d$hosp_id)
  n_total <- nrow(d); n_patients <- n_distinct(d$pat_id)

  # ---- True parameter map for coverage ----
  true_betas <- c("x" = true_beta_x)
  if (incl_teach) true_betas["teaching"] <- beta_teach

  # ---- Coefficients to report (focus on parameters with known truth) ----
  coefs_to_report <- names(true_betas)

  # ---- Universal metadata (added to every output row) ----
  meta_cols <- list(
    experiment = experiment_id, rep = rep_id,
    rho_H = rho_H, rho_P = rho_P, K = K,
    beta_0 = beta_0, p_extra = p_extra,
    site_size_dist = site_size_dist,
    include_teaching = incl_teach,
    n_total = n_total, n_patients = n_patients
  )

  # ---- Helper to build one output row ----
  emit_row <- function(estimator, variant, inference, coef_name,
                       beta_hat, se, df, extra = list()) {
    truth <- true_betas[[coef_name]]

    if (is.na(beta_hat) || is.na(se) || se <= 0) {
      ci_lo <- NA_real_; ci_hi <- NA_real_; covers <- NA
    } else {
      crit <- if (inference == "z" || is.infinite(df)) 1.96 else qt(0.975, df = max(df, 1))
      ci_lo <- beta_hat - crit * se
      ci_hi <- beta_hat + crit * se
      covers <- (ci_lo <= truth) && (truth <= ci_hi)
    }

    base <- c(meta_cols, list(
      estimator = estimator, variant = variant, inference = inference,
      coef = coef_name, true_value = truth,
      beta_hat = beta_hat, se = se, df = df,
      bias = beta_hat - truth, abs_bias = abs(beta_hat - truth),
      ci_lo = ci_lo, ci_hi = ci_hi,
      ci_width = ci_hi - ci_lo, covers = covers
    ))
    for (nm in names(extra)) base[[nm]] <- extra[[nm]]
    as_tibble(base)
  }

  rows <- list()

  # ============================================================================
  # 1. POOLED GEE (centralized) — one row per coef, z-based
  # ============================================================================
  rows$pooled <- tryCatch({
    d$hosp_num <- as.numeric(factor(d$hosp_id))
    fit <- suppressWarnings(
      geeglm(form, id = hosp_num, data = d, family = binomial,
             corstr = "independence"))
    bind_rows(lapply(coefs_to_report, function(cn) {
      idx <- which(names(coef(fit)) == cn)
      if (length(idx) == 0) return(emit_row("pooled_gee", "site_z", "z", cn, NA, NA, Inf))
      emit_row("pooled_gee", "site_z", "z", cn,
               unname(coef(fit)[idx]), unname(sqrt(vcov(fit)[idx, idx])), Inf)
    }))
  }, error = function(e) {
    bind_rows(lapply(coefs_to_report, function(cn)
      emit_row("pooled_gee", "site_z", "z", cn, NA, NA, Inf)))
  })

  # ============================================================================
  # 2. FED-GEE — single fit, multiple variants
  # ============================================================================
  # v3: FG's own df (saws) is only computed for the reported coefficients.
  old_opt <- options(fedgee.df_coefs = coefs_to_report); on.exit(options(old_opt), add = TRUE)
  fed_fit <- tryCatch(
    fedgee_v2(data_list, form, binomial(), "independence", "pat_id",
              n_iter = 50, tol = 1e-8, verbose = FALSE),
    error = function(e) NULL)

  fed_variants_to_emit <- c("patient_z", "site", "kc", "md", "fg", "hc2obs", "vd")

  for (vname in fed_variants_to_emit) {
    for (cn in coefs_to_report) {

      # v3: every site-clustered variant is evaluated under every reference
      # distribution; FG additionally under its own df (dtildeH via saws).
      inference_modes <- switch(vname,
        patient_z = "z",
        vd        = "t_kp",
        fg        = c("z", "t_k1", "t_kp", "t_nu", "t_bm", "t_fgdf"),
        c("z", "t_k1", "t_kp", "t_nu", "t_bm"))

      if (is.null(fed_fit)) {
        for (im in inference_modes) {
          rows[[paste0("fed_", vname, "_", cn, "_", im)]] <-
            emit_row("fed", vname, im, cn, NA, NA, NA)
        }
        next
      }

      idx <- which(names(fed_fit$coefficients) == cn)
      if (length(idx) == 0) next
      v <- fed_fit$variants[[vname]]

      bh <- unname(fed_fit$coefficients[idx])
      se <- unname(v$se[idx])

      for (im in inference_modes) {
        df_use <- switch(im,
                         z      = Inf,
                         t_k1   = v$df_k1,
                         t_kp   = v$df_kp,
                         t_nu   = v$df_nu[idx],
                         t_bm   = v$df_bm[idx],
                         t_fgdf = v$df_fg[idx])
        if (length(df_use) == 0 || is.na(df_use)) df_use <- NA_real_
        rows[[paste0("fed_", vname, "_", cn, "_", im)]] <-
          emit_row("fed", vname, im, cn, bh, se, df_use,
                   extra = list(converged = fed_fit$converged,
                                iterations = fed_fit$iterations))
      }
    }
  }

  # Append diagnostics as separate rows (variant = "diag", coef = "diag")
  if (!is.null(fed_fit)) {
    # Drop diagnostic fields that duplicate meta_cols
    diag_safe <- fed_fit$diagnostics[setdiff(names(fed_fit$diagnostics),
                                             names(meta_cols))]
    diag_row <- as_tibble(c(meta_cols, list(
      estimator = "fed", variant = "diag", inference = "diag",
      coef = "diag", true_value = NA_real_,
      beta_hat = NA_real_, se = NA_real_, df = NA_real_,
      bias = NA_real_, abs_bias = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_,
      ci_width = NA_real_, covers = NA),
      diag_safe))
    rows[["fed_diag"]] <- diag_row
  }

  # ============================================================================
  # 3. META-ANALYTIC GEE  (per-site GEE + inverse-variance pool)
  # ============================================================================
  rows$meta <- tryCatch({
    m <- fit_meta_gee(d, form)
    n_att <- if (is.null(m)) NA_integer_ else attr(m, "n_attempted")
    if (is.null(m)) {
      bind_rows(lapply(coefs_to_report, function(cn)
        emit_row("meta_gee", "guarded", "z", cn, NA, NA, Inf,
                 extra = list(n_sites = NA_integer_, n_dropped = NA_integer_))))
    } else {
      bind_rows(lapply(coefs_to_report, function(cn) {
        r <- m[m$term == cn, ]
        if (nrow(r) == 0) {
          emit_row("meta_gee", "guarded", "z", cn, NA, NA, Inf,
                   extra = list(n_sites = 0L, n_dropped = n_att))
        } else {
          emit_row("meta_gee", "guarded", "z", cn, r$beta, r$se, Inf,
                   extra = list(n_sites = r$n_sites,
                                n_dropped = n_att - r$n_sites))
        }
      }))
    }
  }, error = function(e) {
    bind_rows(lapply(coefs_to_report, function(cn)
      emit_row("meta_gee", "guarded", "z", cn, NA, NA, Inf,
               extra = list(n_sites = NA_integer_, n_dropped = NA_integer_))))
  })

  # ============================================================================
  # 4. GLMM (CONDITIONAL on grid flag — gates expensive fit)
  # ============================================================================
  if (fit_glmm_flag) {
    rows$glmm <- tryCatch({
      form_glmm_chr <- paste0(paste(deparse(form), collapse = ""),
                              " + (1 | hosp_id) + (1 | pat_id)")
      form_glmm <- as.formula(form_glmm_chr)
      fit <- glmer(form_glmm, data = d, family = binomial, nAGQ = 1,
                   control = glmerControl(optimizer = "bobyqa",
                                          optCtrl = list(maxfun = 1e5)))
      bind_rows(lapply(coefs_to_report, function(cn) {
        fe <- fixef(fit)
        idx <- which(names(fe) == cn)
        if (length(idx) == 0) return(emit_row("glmm", "naive", "z", cn, NA, NA, Inf))
        se_glmm <- sqrt(diag(vcov(fit)))[idx]
        emit_row("glmm", "naive", "z", cn,
                 unname(fe[idx]), unname(se_glmm), Inf,
                 extra = list(
                   sigma2_hosp = VarCorr(fit)$hosp_id[1, 1],
                   sigma2_pat  = VarCorr(fit)$pat_id[1, 1]
                 ))
      }))
    }, error = function(e) {
      bind_rows(lapply(coefs_to_report, function(cn)
        emit_row("glmm", "naive", "z", cn, NA, NA, Inf)))
    })
  }

  bind_rows(rows)
}

###############################################################################
# E. PARAMETER GRID — extended for E1-E8
#
# Each experiment is a separate block of grid rows tagged with `experiment`.
# Build the grid by composing experiment-specific blocks.
###############################################################################

# ---- Default grid row template ----
grid_row_default <- function() {
  data.frame(
    rho_H = 0.05, rho_P = 0.30, K = 50,
    pat_lo = 10, pat_hi = 30, vis_lo = 2, vis_hi = 6,
    beta_0 = -1, true_beta_x = 1,
    include_teaching = FALSE, beta_teach = 0,
    p_extra = 0, site_size_dist = "uniform",
    fit_glmm = FALSE,
    stringsAsFactors = FALSE
  )
}

# ---- Helper: append default columns (those not already in base) ----
.append_defaults <- function(base) {
  defaults <- grid_row_default()
  missing_cols <- setdiff(names(defaults), names(base))
  for (cn in missing_cols) base[[cn]] <- defaults[[cn]]
  base
}

# =============================================================================
# CALIBRATION EXPERIMENT (EC)
#
# Goal: characterize how the fed-GEE estimator and its small-sample variance
# corrections (KC, MD) behave -- across reference distributions (z, t_{K-p},
# t_Satterthwaite) -- as a function of:
#   * site count K            (the cluster count the sandwich/df depend on)
#   * predictor count p       (drives residual df K - p)
#   * outcome prevalence      (rare / moderate / common, via beta_0)
#   * within/between ICC      (rho_H)
# with pooled GEE (oracle) and per-site meta-GEE (privacy-preserving) as
# reference arms.
#
# Prevalence is set through the intercept beta_0 (x ~ N(0,1), true_beta_x = 1):
#   beta_0 = -3  -> marginal prevalence ~ 6%   (rare)
#   beta_0 = -1  -> marginal prevalence ~ 27%  (moderate)
#   beta_0 =  0  -> marginal prevalence ~ 50%  (common)
#
# Feasibility filter: keep only cells with residual df = K - p >= 5, where
# p = p_extra + 2 (intercept + x). This drops e.g. K=10 with p>=10, keeping
# every cell estimable and every t_{K-p} reference well defined.
# =============================================================================
build_grid_EC <- function() {
  base <- expand.grid(
    K       = c(10, 20, 30, 50, 75, 100),
    p_extra = c(0, 3, 8, 18),          # total p = 2, 5, 10, 20
    rho_H   = c(0, 0.05, 0.10),
    beta_0  = c(-3, -1, 0),            # rare / moderate / common
    rho_P   = 0.30,
    stringsAsFactors = FALSE)

  # Residual df = K - (p_extra + 2); require >= 5 for a usable t-reference.
  base <- base[(base$K - (base$p_extra + 2)) >= 5, , drop = FALSE]

  base <- .append_defaults(base)          # beta_0 already set -> preserved
  base$experiment      <- "EC"
  base$true_beta_x     <- 1
  base$include_teaching <- FALSE
  base$site_size_dist  <- "uniform"       # small sites keep the small-K regime real
  base$pat_lo <- 10; base$pat_hi <- 30
  base$vis_lo <- 2;  base$vis_hi <- 6
  base$fit_glmm        <- FALSE
  rownames(base) <- NULL
  base
}

# =============================================================================
# EXPERIMENT EV3 (2026-09-27): true site-level FG + Bell-McCaffrey df, with
# UNEQUAL site sizes in the small-K regime (never simulated in v2/calibration).
#
#   K         in {10, 18, 30, 50}      (18 = the 18-VISN VHA application)
#   p total   in {2, 5, 10, 20}        (keep K - p >= 5)
#   beta_0    in {-3, -1, 0}           (~6% / ~27% / ~50% prevalence)
#   site size in {uniform (CV~0.30), cv078 (VHA-matched), cv100}
#                (all with mean ~20 patients/site)
#   rho_H = 0.05, rho_P = 0.30
#
# Variance variants: site (uncorrected), KC, MD, FG (site-level), hc2obs (v2
# "fg"), each under z / t(K-1) / t(K-p) / t(nu) / t(BM); FG also t(FG df).
# =============================================================================
build_grid_EV3 <- function() {
  base <- expand.grid(
    K              = c(10, 18, 30, 50),
    p_extra        = c(0, 3, 8, 18),         # total p = 2, 5, 10, 20
    beta_0         = c(-3, -1, 0),
    site_size_dist = c("uniform", "cv078", "cv100"),
    stringsAsFactors = FALSE)
  base <- base[(base$K - (base$p_extra + 2)) >= 5, , drop = FALSE]
  base <- .append_defaults(base)
  base$experiment       <- "EV3"
  base$rho_H            <- 0.05
  base$rho_P            <- 0.30
  base$true_beta_x      <- 1
  base$include_teaching <- FALSE
  base$pat_lo <- 10; base$pat_hi <- 30
  base$vis_lo <- 2;  base$vis_hi <- 6
  base$fit_glmm         <- FALSE
  rownames(base) <- NULL
  base
}

# ---- Master grid ----
build_param_grid <- function(experiments = c("EV3")) {
  grids <- list()
  if ("EC"  %in% experiments) grids$EC  <- build_grid_EC()
  if ("EV3" %in% experiments) grids$EV3 <- build_grid_EV3()

  grid <- bind_rows(grids)
  grid$scenario_id <- seq_len(nrow(grid))
  grid
}

###############################################################################
# F. TASK TABLE
###############################################################################
build_task_table <- function(n_reps = 200, experiments = c("EV3")) {
  grid <- build_param_grid(experiments)
  tasks <- grid[rep(seq_len(nrow(grid)), each = n_reps), ]
  tasks$rep_id <- rep(seq_len(n_reps), times = nrow(grid))
  tasks$task_id <- seq_len(nrow(tasks))
  rownames(tasks) <- NULL
  tasks
}
