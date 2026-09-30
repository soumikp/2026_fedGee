###############################################################################
#
# FedGEE Simulation Framework
#
# Design:
#   1. Parameter grid defines ALL scenarios (rho_H, rho_P, K, etc.)
#   2. run_one_rep() is a self-contained function: takes one parameter row +
#      a rep ID, returns a tidy data.frame with one row per METHOD
#   3. Runner script dispatches rows to cluster nodes via parallel/future
#   4. Results saved as one RDS per task (crash-safe)
#   5. Aggregation script combines RDS files and produces tables/figures
#
# Files:
#   00_simulation_config.R   <- this file (grid + single-rep function)
#   01_fedgee.R               <- centralized/pooled FedGEE (source'd; defines
#                                fedgee(), used below by the fed_site /
#                                fed_site_md / fed_patient benchmark methods)
#   01_decentralized_fedgee.R <- decentralized Dec-Fed-GEE (source'd; defines
#                                decentralized_fedgee(), used below by the
#                                dec_hub / dec_ring / dec_visn / dec_complete
#                                methods)
#   02_run_cluster.R          <- cluster dispatch script
#   03_aggregate_results.R    <- post-hoc analysis
#
# run_one_rep() now returns up to 10 rows per replication: the 6 centralized
# benchmark methods (pooled_gee, fed_site, fed_site_md, fed_patient,
# meta_glm, glmm) PLUS 4 decentralized methods, one per network topology
# (dec_hub, dec_ring, dec_visn, dec_complete), all fit to the SAME simulated
# dataset `d`.
#
# NOTE ON DECENTRALIZED REPORTING (updated):
#   decentralized_fedgee() now returns the network-average estimator
#   beta_dec = (1/K) sum_i beta_hat_i directly as $coefficients (a named
#   vector), with $vcov and $se its averaged sandwich covariance / SE --
#   see Section 3.2 of the Dec-Fed-GEE theory note. dec_row() below reads
#   these fields directly instead of manually averaging
#   $coefficients_by_site / $se_by_site itself. The per-site matrices are
#   still returned by decentralized_fedgee() and are used here only to
#   report cross-site disagreement diagnostics (beta_sd_across_sites, etc.),
#   never as an alternative point estimate.
###############################################################################

library(geepack)
library(lme4)
library(dplyr)
library(purrr)
library(tibble)
library(Matrix)
library(FedGEE)

###############################################################################
# A. DATA-GENERATING PROCESS
###############################################################################
generate_sim_data <- function(rho_H, rho_P,
                              n_hospitals,
                              patients_per_hosp_range = c(10, 30),
                              visits_per_patient_range = c(2, 6),
                              beta_0 = -1,           
                              true_beta_x = 1,   
                              include_teaching = FALSE,
                              beta_teach = 0.5) {
  
  stopifnot(rho_P >= rho_H, rho_H >= 0, rho_P < 1)
  
  U <- rnorm(n_hospitals)
  
  teaching <- rep(0L, n_hospitals)
  if (include_teaching) {
    teaching[sample(1:n_hospitals, floor(n_hospitals / 2))] <- 1L
  }
  
  # Variable number of patients per hospital
  n_pat_per_hosp <- sample(
    patients_per_hosp_range[1]:patients_per_hosp_range[2],
    n_hospitals, replace = TRUE
  )
  
  # Build patient table
  patients <- data.frame(
    hosp_id = rep(1:n_hospitals, times = n_pat_per_hosp),
    pat_num = unlist(lapply(n_pat_per_hosp, seq_len))
  )
  patients$pat_id <- paste0(patients$hosp_id, "_", patients$pat_num)
  
  # Generate visits
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
    
    data.frame(hosp_id  = h,
               pat_id   = patients$pat_id[p],
               teaching = teaching[h],
               x = x, y = y)
  })
  
  bind_rows(visits_list) |> arrange(hosp_id, pat_id)
}

###############################################################################
# B. HELPER: Meta-analytic GLM
###############################################################################
fit_meta_glm <- function(data) {
  sites <- split(data, data$hosp_id)
  site_fits <- lapply(sites, function(sd) {
    tryCatch({
      fit <- suppressWarnings(glm(y ~ x, data = sd, family = binomial))
      if (!fit$converged || any(is.na(coef(fit)))) return(NULL)
      data.frame(beta = coef(fit)["x"], se = sqrt(vcov(fit)["x", "x"]))
    }, error = function(e) NULL)
  })
  site_fits <- bind_rows(Filter(Negate(is.null), site_fits))
  site_fits <- site_fits[is.finite(site_fits$se) & site_fits$se > 0 &
                           is.finite(site_fits$beta), ]
  if (nrow(site_fits) < 2) return(list(beta = NA, se = NA))
  w    <- 1 / site_fits$se^2
  list(beta = sum(w * site_fits$beta) / sum(w), se = 1 / sqrt(sum(w)))
}

###############################################################################
# B2. HELPER: auto-partition K sites into VISN regions
#
# decentralized_fedgee(..., structure = "visn") needs a `region_id` (length
# K, which region each site belongs to) and `hub_sites` (indices of each
# region's hub). Since the simulation grid only specifies K, this derives a
# reasonable regional network automatically: sequential blocks of
# `region_size` sites per region (last region gets the remainder), with the
# FIRST site in each region as its hub. Regional hubs are always mutually
# connected (see build_adj_visn()), so the resulting network is connected
# for any K >= 1 and any region_size >= 1.
###############################################################################
make_visn_partition <- function(K, region_size = 5) {
  n_regions <- max(1L, ceiling(K / region_size))
  region_id <- rep(seq_len(n_regions), each = region_size, length.out = K)
  hub_sites <- vapply(seq_len(n_regions), function(r) which(region_id == r)[1], integer(1))
  list(region_id = region_id, hub_sites = hub_sites)
}

###############################################################################
# C. SINGLE REPLICATION FUNCTION
#
# Input:  one row of parameter grid + rep ID
# Output: data.frame with one row per METHOD, columns for all metrics
#
# Methods:
#   1. pooled_gee        — geeglm on full data, hospital-clustered
#   2. fed_site          — FedGEE, site-level sandwich
#   3. fed_site_md       — FedGEE, site-level + Mancl-DeRouen
#   4. fed_patient       — FedGEE, patient-level sandwich
#   5. meta_glm          — site-level GLM + inverse-variance meta-analysis
#   6. glmm              — GLMM with hospital + patient random intercepts
#   7. dec_hub           — Decentralized Dec-Fed-GEE, hub-and-spoke topology
#   8. dec_ring          — Decentralized Dec-Fed-GEE, ring topology
#   9. dec_visn          — Decentralized Dec-Fed-GEE, VISN (regional) topology
#  10. dec_complete       — Decentralized Dec-Fed-GEE, complete (fully
#                          connected) topology. This is the "every site talks
#                          to every other site every round" benchmark: with
#                          W = (1/K) * J, one consensus round already equals
#                          the exact network average, so this topology acts
#                          as a best-case reference for how fast dec_hub /
#                          dec_ring / dec_visn approach it as L_beta, L_S,
#                          L_B grow.
#
# dec_hub / dec_ring / dec_visn / dec_complete report the network-average
# decentralized estimator beta_dec = (1/K) sum_i beta_hat_i as
# decentralized_fedgee()'s $coefficients / $se / $vcov (see Section 3.2 of
# the theory note); this is used directly below, not recomputed. Per-site
# rows ($coefficients_by_site / $se_by_site) are used only to report
# cross-site disagreement as a diagnostic.
###############################################################################
run_one_rep <- function(params,
                        rep_id,
                        dec_n_iter = 100L,
                        dec_L_init = 50L,
                        dec_L_beta = 50L,
                        dec_L_S = 50L,
                        dec_L_B = 50L,
                        dec_tol = 1e-8,
                        dec_tol_update = 1e-8,
                        dec_tol_consensus = 1e-8,
                        dec_tol_score = 1e-8,
                        dec_step_size = 0.5,
                        dec_sandwich_level = "site",
                        dec_md_correction = FALSE,
                        dec_ridge = 1e-5,
                        dec_ring_K_neigh = 2L,
                        dec_hub = 1L,
                        dec_visn_region_size = 5L,
                        dec_verbose = FALSE) {
  
  rho_H <- params$rho_H
  rho_P <- params$rho_P
  K <- params$K
  patient_range <- c(params$pat_lo, params$pat_hi)
  visit_range <- c(params$vis_lo, params$vis_hi)
  include_teaching <- isTRUE(params$include_teaching)
  beta_teach <- params$beta_teach
  true_beta_x <- params$true_beta_x
  beta_0 <- if (is.null(params$beta_0) || !is.finite(params$beta_0)) {
    -1
  } else {
    params$beta_0
  }
  
  main_formula <- if (include_teaching) {
    y ~ x + teaching
  } else {
    y ~ x
  }
  
  data <- generate_sim_data(
    rho_H = rho_H,
    rho_P = rho_P,
    n_hospitals = K,
    patients_per_hosp_range = patient_range,
    visits_per_patient_range = visit_range,
    beta_0 = beta_0,
    true_beta_x = true_beta_x,
    include_teaching = include_teaching,
    beta_teach = beta_teach
  )
  
  data_list <- split(data, data$hosp_id)
  n_total <- nrow(data)
  n_patients <- n_distinct(data$pat_id)
  
  meta_row <- function(method,
                       beta_hat,
                       se_hat,
                       beta_teach_hat = NA_real_,
                       se_teach_hat = NA_real_,
                       extra = list()) {
    beta_hat <- as.numeric(beta_hat)[1L]
    se_hat <- as.numeric(se_hat)[1L]
    beta_teach_hat <- as.numeric(beta_teach_hat)[1L]
    se_teach_hat <- as.numeric(se_teach_hat)[1L]
    
    output <- data.frame(
      rho_H = rho_H,
      rho_P = rho_P,
      K = K,
      pat_lo = patient_range[1L],
      pat_hi = patient_range[2L],
      vis_lo = visit_range[1L],
      vis_hi = visit_range[2L],
      include_teaching = include_teaching,
      beta_0 = beta_0,
      rep = rep_id,
      n_total = n_total,
      n_patients = n_patients,
      method = method,
      beta_hat = beta_hat,
      se_hat = se_hat,
      beta_teach_hat = beta_teach_hat,
      se_teach_hat = se_teach_hat,
      bias = beta_hat - true_beta_x,
      abs_bias = abs(beta_hat - true_beta_x),
      ci_lo_z = beta_hat - 1.96 * se_hat,
      ci_hi_z = beta_hat + 1.96 * se_hat,
      covers_z = (beta_hat - 1.96 * se_hat <= true_beta_x) &
        (true_beta_x <= beta_hat + 1.96 * se_hat),
      ci_width_z = 2 * 1.96 * se_hat,
      stringsAsFactors = FALSE
    )
    
    for (name in names(extra)) {
      output[[name]] <- extra[[name]]
    }
    
    output
  }
  
  error_row <- function(method, error_message, extra = list()) {
    meta_row(
      method = method,
      beta_hat = NA_real_,
      se_hat = NA_real_,
      beta_teach_hat = NA_real_,
      se_teach_hat = NA_real_,
      extra = c(
        list(
          converged = FALSE,
          iterations = NA_integer_,
          error_msg = as.character(error_message)
        ),
        extra
      )
    )
  }
  
  results <- list()
  
  # ---------------------------------------------------------------------------
  # 1. Pooled GEE, treating hospital as the independent clustering unit.
  # ---------------------------------------------------------------------------
  results$pooled <- tryCatch({
    data$hosp_num <- as.integer(factor(data$hosp_id))
    fit <- geeglm(
      formula = main_formula,
      id = hosp_num,
      data = data,
      family = binomial(),
      corstr = "independence"
    )
    
    covariance <- vcov(fit)
    beta_teaching_hat <- if (include_teaching) {
      unname(coef(fit)["teaching"])
    } else {
      NA_real_
    }
    se_teaching_hat <- if (include_teaching) {
      sqrt(covariance["teaching", "teaching"])
    } else {
      NA_real_
    }
    
    meta_row(
      method = "pooled_gee",
      beta_hat = unname(coef(fit)["x"]),
      se_hat = sqrt(covariance["x", "x"]),
      beta_teach_hat = beta_teaching_hat,
      se_teach_hat = se_teaching_hat,
      extra = list(converged = TRUE, iterations = NA_integer_, error_msg = NA_character_)
    )
  }, error = function(e) {
    error_row("pooled_gee", conditionMessage(e))
  })
  
  # ---------------------------------------------------------------------------
  # 2. Centralized FedGEE: site-level sandwich, no small-sample correction.
  # ---------------------------------------------------------------------------
  results$fed_site <- tryCatch({
    fit <- FedGEE::fedgee(
      data_list,
      main_formula,
      binomial(),
      "independence",
      "pat_id",
      sandwich_level = "site",
      correction = "none",
      verbose = FALSE
    )
    
    df_residual <- fit$df_residual
    t_critical <- qt(0.975, df = max(df_residual, 1L))
    beta_x <- unname(fit$coefficients["x"])
    se_x <- unname(fit$se["x"])
    
    meta_row(
      method = "fed_site",
      beta_hat = beta_x,
      se_hat = se_x,
      beta_teach_hat = if (include_teaching) unname(fit$coefficients["teaching"]) else NA_real_,
      se_teach_hat = if (include_teaching) unname(fit$se["teaching"]) else NA_real_,
      extra = list(
        converged = fit$converged,
        iterations = fit$iterations,
        df = df_residual,
        ci_lo_t = beta_x - t_critical * se_x,
        ci_hi_t = beta_x + t_critical * se_x,
        covers_t = (beta_x - t_critical * se_x <= true_beta_x) &
          (true_beta_x <= beta_x + t_critical * se_x),
        ci_width_t = 2 * t_critical * se_x,
        error_msg = NA_character_
      )
    )
  }, error = function(e) {
    error_row("fed_site", conditionMessage(e))
  })
  
  # ---------------------------------------------------------------------------
  # 3. Centralized FedGEE: site-level leverage correction.
  # The existing benchmark method name is retained as fed_site_md, while the
  # package call uses correction = "KC" to match its implemented option.
  # ---------------------------------------------------------------------------
  results$fed_site_md <- tryCatch({
    fit <- FedGEE::fedgee(
      data_list,
      main_formula,
      binomial(),
      "independence",
      "pat_id",
      sandwich_level = "site",
      correction = "KC",
      verbose = FALSE
    )
    
    df_residual <- fit$df_residual
    t_critical <- qt(0.975, df = max(df_residual, 1L))
    beta_x <- unname(fit$coefficients["x"])
    se_x <- unname(fit$se["x"])
    
    meta_row(
      method = "fed_site_md",
      beta_hat = beta_x,
      se_hat = se_x,
      beta_teach_hat = if (include_teaching) unname(fit$coefficients["teaching"]) else NA_real_,
      se_teach_hat = if (include_teaching) unname(fit$se["teaching"]) else NA_real_,
      extra = list(
        converged = fit$converged,
        iterations = fit$iterations,
        df = df_residual,
        ci_lo_t = beta_x - t_critical * se_x,
        ci_hi_t = beta_x + t_critical * se_x,
        covers_t = (beta_x - t_critical * se_x <= true_beta_x) &
          (true_beta_x <= beta_x + t_critical * se_x),
        ci_width_t = 2 * t_critical * se_x,
        error_msg = NA_character_
      )
    )
  }, error = function(e) {
    error_row("fed_site_md", conditionMessage(e))
  })
  
  # ---------------------------------------------------------------------------
  # 4. Centralized FedGEE: patient-level sandwich.
  # ---------------------------------------------------------------------------
  results$fed_patient <- tryCatch({
    fit <- FedGEE::fedgee(
      data_list,
      main_formula,
      binomial(),
      "independence",
      "pat_id",
      sandwich_level = "patient",
      correction = "none",
      verbose = FALSE
    )
    
    meta_row(
      method = "fed_patient",
      beta_hat = unname(fit$coefficients["x"]),
      se_hat = unname(fit$se["x"]),
      beta_teach_hat = if (include_teaching) unname(fit$coefficients["teaching"]) else NA_real_,
      se_teach_hat = if (include_teaching) unname(fit$se["teaching"]) else NA_real_,
      extra = list(
        converged = fit$converged,
        iterations = fit$iterations,
        error_msg = NA_character_
      )
    )
  }, error = function(e) {
    error_row("fed_patient", conditionMessage(e))
  })
  
  # ---------------------------------------------------------------------------
  # 5. Inverse-variance meta-analysis of local GLMs.
  # ---------------------------------------------------------------------------
  results$meta_glm <- tryCatch({
    meta_fit <- fit_meta_glm(data)
    meta_row(
      method = "meta_glm",
      beta_hat = meta_fit$beta,
      se_hat = meta_fit$se,
      extra = list(converged = TRUE, iterations = NA_integer_, error_msg = NA_character_)
    )
  }, error = function(e) {
    error_row("meta_glm", conditionMessage(e))
  })
  
  # ---------------------------------------------------------------------------
  # 6. GLMM with hospital and patient random intercepts.
  # ---------------------------------------------------------------------------
  results$glmm <- tryCatch({
    glmm_formula <- if (include_teaching) {
      y ~ x + teaching + (1 | hosp_id) + (1 | pat_id)
    } else {
      y ~ x + (1 | hosp_id) + (1 | pat_id)
    }
    
    fit <- glmer(
      formula = glmm_formula,
      data = data,
      family = binomial(),
      nAGQ = 1,
      control = glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 1e5)
      )
    )
    
    fixed_covariance <- vcov(fit)
    fixed_se <- sqrt(diag(fixed_covariance))
    
    meta_row(
      method = "glmm",
      beta_hat = unname(fixef(fit)["x"]),
      se_hat = unname(fixed_se["x"]),
      beta_teach_hat = if (include_teaching) unname(fixef(fit)["teaching"]) else NA_real_,
      se_teach_hat = if (include_teaching) unname(fixed_se["teaching"]) else NA_real_,
      extra = list(
        converged = TRUE,
        iterations = NA_integer_,
        sigma2_hosp = as.numeric(VarCorr(fit)$hosp_id[1L, 1L]),
        sigma2_pat = as.numeric(VarCorr(fit)$pat_id[1L, 1L]),
        error_msg = NA_character_
      )
    )
  }, error = function(e) {
    error_row("glmm", conditionMessage(e))
  })
  
  # ---------------------------------------------------------------------------
  # Shared reporting helper for decentralized fits. Every decentralized method
  # uses the consensus meta-GLM warm start, a damped Fisher step, and a small
  # ridge to stabilize difficult scenarios.
  #
  # beta_hat and se_hat below are read directly from decentralized_fedgee()'s
  # $coefficients / $se -- the network-average estimator
  # beta_dec = (1/K) sum_i beta_hat_i and its averaged sandwich covariance
  # (Section 3.2 of the theory note). No averaging is performed in this
  # helper. Per-site rows ($coefficients_by_site, $se_by_site,
  # $coefficients_before_final_consensus) are read ONLY to build diagnostic
  # summaries of cross-site disagreement (beta_sd_across_sites, etc.); they
  # are never used as an alternative point estimate.
  # ---------------------------------------------------------------------------
  dec_row <- function(fit, method_tag, structure_tag) {
    if (is.null(fit)) {
      stop("decentralized_fedgee() returned NULL.", call. = FALSE)
    }
    
    if (!("x" %in% names(fit$coefficients)) ||
        !("x" %in% names(fit$se))) {
      stop("The decentralized fit does not contain coefficient 'x'.", call. = FALSE)
    }
    
    # --- Primary point estimate and SE: beta_dec / se_dec, used as-is. ---
    beta_x <- unname(fit$coefficients["x"])
    se_x <- unname(fit$se["x"])
    
    beta_teaching <- if (include_teaching && "teaching" %in% names(fit$coefficients)) {
      unname(fit$coefficients["teaching"])
    } else {
      NA_real_
    }
    
    se_teaching <- if (include_teaching && "teaching" %in% names(fit$se)) {
      unname(fit$se["teaching"])
    } else {
      NA_real_
    }
    
    # --- Diagnostics only, from per-site rows. ---
    beta_x_by_site <- if ("x" %in% colnames(fit$coefficients_by_site)) {
      fit$coefficients_by_site[, "x"]
    } else {
      rep(NA_real_, fit$n_sites)
    }
    se_x_by_site <- if ("x" %in% colnames(fit$se_by_site)) {
      fit$se_by_site[, "x"]
    } else {
      rep(NA_real_, fit$n_sites)
    }
    
    # decentralized_fedgee() runs one extra L_beta-round parameter consensus
    # AFTER the training loop stops (see 01_decentralized_fedgee.R, "Final
    # parameter consensus" step) before beta_dec is formed.
    # fit$coefficients_before_final_consensus is the per-site beta right
    # before that pass; comparing its network average with beta_dec
    # quantifies how much that last smoothing step is doing.
    beta_x_before_final <- if (!is.null(fit$coefficients_before_final_consensus) &&
                               "x" %in% colnames(fit$coefficients_before_final_consensus)) {
      fit$coefficients_before_final_consensus[, "x"]
    } else {
      rep(NA_real_, length(beta_x_by_site))
    }
    beta_x_before_final_mean <- mean(beta_x_before_final, na.rm = TRUE)
    
    df_residual <- fit$df_residual
    t_critical <- qt(0.975, df = max(df_residual, 1L))
    
    final_diagnostics <- fit$final_diagnostics
    final_update_error <- if (!is.null(final_diagnostics)) {
      final_diagnostics$max_update_error[1L]
    } else {
      NA_real_
    }
    final_consensus_error <- if (!is.null(final_diagnostics)) {
      final_diagnostics$max_consensus_error[1L]
    } else {
      NA_real_
    }
    final_score_error <- if (!is.null(final_diagnostics)) {
      final_diagnostics$max_score_error[1L]
    } else {
      NA_real_
    }
    
    meta_row(
      method = method_tag,
      beta_hat = beta_x,
      se_hat = se_x,
      beta_teach_hat = beta_teaching,
      se_teach_hat = se_teaching,
      extra = list(
        converged = fit$converged,
        iterations = fit$iterations,
        df = df_residual,
        total_clusters = fit$total_clusters,
        ci_lo_t = beta_x - t_critical * se_x,
        ci_hi_t = beta_x + t_critical * se_x,
        covers_t = (beta_x - t_critical * se_x <= true_beta_x) &
          (true_beta_x <= beta_x + t_critical * se_x),
        ci_width_t = 2 * t_critical * se_x,
        structure = structure_tag,
        sandwich_level = fit$sandwich_level,
        md_correction = fit$md_correction,
        L_init = fit$L_init,
        init_intercept_sd = if (!is.null(fit$initial_beta_by_site) &&
                                "(Intercept)" %in% colnames(fit$initial_beta_by_site)) {
          sd(fit$initial_beta_by_site[, "(Intercept)"], na.rm = TRUE)
        } else {
          NA_real_
        },
        init_intercept_range = if (!is.null(fit$initial_beta_by_site) &&
                                   "(Intercept)" %in% colnames(fit$initial_beta_by_site)) {
          diff(range(fit$initial_beta_by_site[, "(Intercept)"], na.rm = TRUE))
        } else {
          NA_real_
        },
        L_beta = fit$L_beta,
        L_S = fit$L_S,
        L_B = fit$L_B,
        tol_update = fit$tol_update,
        tol_consensus = fit$tol_consensus,
        tol_score = fit$tol_score,
        step_size = if (length(fit$step_size) == 1L && is.numeric(fit$step_size)) {
          fit$step_size
        } else {
          NA_real_
        },
        ridge = fit$ridge,
        n_valid_sites = fit$n_valid_sites,
        beta_sd_across_sites = sd(beta_x_by_site, na.rm = TRUE),
        se_sd_across_sites = sd(se_x_by_site, na.rm = TRUE),
        beta_range_across_sites = if (any(is.finite(beta_x_by_site))) {
          diff(range(beta_x_by_site[is.finite(beta_x_by_site)]))
        } else {
          NA_real_
        },
        beta_sd_before_final_consensus = sd(beta_x_before_final, na.rm = TRUE),
        beta_shift_from_final_consensus = beta_x - beta_x_before_final_mean,
        final_update_error = final_update_error,
        final_consensus_error = final_consensus_error,
        final_score_error = final_score_error,
        error_msg = NA_character_
      )
    )
  }
  
  dec_error_extra <- function(structure_tag) {
    list(
      structure = structure_tag,
      sandwich_level = dec_sandwich_level,
      md_correction = dec_md_correction,
      L_init = dec_L_init,
      init_intercept_sd = NA_real_,
      init_intercept_range = NA_real_,
      L_beta = dec_L_beta,
      L_S = dec_L_S,
      L_B = dec_L_B,
      tol_update = dec_tol_update,
      tol_consensus = dec_tol_consensus,
      tol_score = dec_tol_score,
      step_size = if (length(dec_step_size) == 1L && is.numeric(dec_step_size)) {
        dec_step_size
      } else {
        NA_real_
      },
      ridge = dec_ridge,
      total_clusters = NA_integer_,
      beta_sd_before_final_consensus = NA_real_,
      beta_shift_from_final_consensus = NA_real_
    )
  }
  
  
  # Common decentralized arguments. Named arguments are used deliberately so
  # future changes in argument order do not silently change the simulation.
  decentralized_common <- list(
    data_list = data_list,
    main_formula = main_formula,
    family_obj = binomial(link = "logit"),
    corstr = "independence",
    id_col = "pat_id",
    W = NULL,
    site_constant_cols = if (include_teaching) "teaching" else character(0),
    L_init = dec_L_init,
    L_beta = dec_L_beta,
    L_S = dec_L_S,
    L_B = dec_L_B,
    sandwich_level = dec_sandwich_level,
    md_correction = dec_md_correction,
    n_iter = dec_n_iter,
    tol = dec_tol,
    tol_update = dec_tol_update,
    tol_consensus = dec_tol_consensus,
    tol_score = dec_tol_score,
    step_size = dec_step_size,
    ridge = dec_ridge,
    verbose = dec_verbose
  )
  
  # ---------------------------------------------------------------------------
  # 7. Decentralized FedGEE: hub-and-spoke topology.
  # ---------------------------------------------------------------------------
  results$dec_hub <- tryCatch({
    fit <- do.call(
      decentralized_fedgee,
      c(
        decentralized_common,
        list(structure = "hub", hub = dec_hub)
      )
    )
    dec_row(fit, "dec_hub", "hub")
  }, error = function(e) {
    error_row(
      "dec_hub",
      conditionMessage(e),
      extra = dec_error_extra("hub")
    )
  })
  
  # ---------------------------------------------------------------------------
  # 8. Decentralized FedGEE: ring topology.
  # ---------------------------------------------------------------------------
  results$dec_ring <- tryCatch({
    fit <- do.call(
      decentralized_fedgee,
      c(
        decentralized_common,
        list(structure = "ring", K_neigh = dec_ring_K_neigh)
      )
    )
    dec_row(fit, "dec_ring", "ring")
  }, error = function(e) {
    error_row(
      "dec_ring",
      conditionMessage(e),
      extra = dec_error_extra("ring")
    )
  })
  
  # ---------------------------------------------------------------------------
  # 9. Decentralized FedGEE: VISN topology.
  # ---------------------------------------------------------------------------
  results$dec_visn <- tryCatch({
    visn <- make_visn_partition(K, region_size = dec_visn_region_size)
    fit <- do.call(
      decentralized_fedgee,
      c(
        decentralized_common,
        list(
          structure = "visn",
          region_id = visn$region_id,
          hub_sites = visn$hub_sites
        )
      )
    )
    dec_row(fit, "dec_visn", "visn")
  }, error = function(e) {
    error_row(
      "dec_visn",
      conditionMessage(e),
      extra = dec_error_extra("visn")
    )
  })
  
  # ---------------------------------------------------------------------------
  # 10. Decentralized FedGEE: complete (fully connected) topology.
  #
  # W = (1/K) * J (every entry 1/K). With a fully connected network every
  # site is a neighbor of every other site, so this is the fastest-mixing
  # topology available and serves as a best-case reference: with even
  # L_beta = L_S = L_B = 1, one consensus round already reproduces the exact
  # network average. Comparing dec_hub / dec_ring / dec_visn against
  # dec_complete at matched L_* is a natural way to quantify the statistical
  # cost of sparser communication graphs.
  # ---------------------------------------------------------------------------
  results$dec_complete <- tryCatch({
    fit <- do.call(
      decentralized_fedgee,
      c(
        decentralized_common,
        list(structure = "complete")
      )
    )
    dec_row(fit, "dec_complete", "complete")
  }, error = function(e) {
    error_row(
      "dec_complete",
      conditionMessage(e),
      extra = dec_error_extra("complete")
    )
  })
  
  bind_rows(results)
}
###############################################################################
# D. PARAMETER GRID
#
# One row = one unique simulation scenario.
# Each scenario is repeated n_reps times.
###############################################################################
build_param_grid <- function() {
  
  # Core grid: rho_H x rho_P x K
  core <- expand.grid(
    rho_H = c(0.00,0.10,0.20),
    rho_P = c(0.25,0.30),
    K     = c(5, 10, 15, 25, 50, 100),
    stringsAsFactors = FALSE
  )
  # Enforce constraint: rho_P > rho_H + small gap
  core <- core[core$rho_P > core$rho_H + 0.05, ]
  
  # Patient and visit ranges
  core$pat_lo <- 10
  core$pat_hi <- 30
  core$vis_lo <- 2
  core$vis_hi <- 6
  
  # Default intercept and slope
  core$beta_0      <- -1
  core$true_beta_x <- 1
  
  # ---------------------------------------------------------------------------
  # 1. Standard scenarios (No teaching, common outcome)
  # ---------------------------------------------------------------------------
  core_no_teach <- core |> 
    mutate(include_teaching = FALSE, beta_teach = 0)
  
  # ---------------------------------------------------------------------------
  # 2. Teaching scenarios (Site-level covariate)
  # ---------------------------------------------------------------------------
  core_teach <- core |>
    filter(K >= 15) |>  # Need enough sites for teaching to be identified
    mutate(include_teaching = TRUE, beta_teach = 0.5)
  
  # ---------------------------------------------------------------------------
  # 3. Rare outcome scenarios (beta_0 = -3)
  # ---------------------------------------------------------------------------
  core_rare <- core |> 
    mutate(beta_0 = -3, include_teaching = FALSE, beta_teach = 0)
  
  # Combine all scenarios
  grid <- bind_rows(core_no_teach, core_teach, core_rare)
  grid$scenario_id <- seq_len(nrow(grid))
  
  return(grid)
}
###############################################################################
# E. TASK TABLE
#
# Expand grid by n_reps. Each row = one call to run_one_rep().
# This is what gets dispatched to the cluster.
###############################################################################
build_task_table <- function(n_reps = 200) {
  grid <- build_param_grid()
  tasks <- grid[rep(seq_len(nrow(grid)), each = n_reps), ]
  tasks$rep_id <- rep(seq_len(n_reps), times = nrow(grid))
  tasks$task_id <- seq_len(nrow(tasks))
  rownames(tasks) <- NULL
  tasks
}
