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
#   01_function_fedgee_improved.R  <- FedGEE implementation (source'd)
#   02_run_cluster.R         <- cluster dispatch script
#   03_aggregate_results.R   <- post-hoc analysis
#
###############################################################################

library(geepack)
library(lme4)
library(dplyr)
library(purrr)
library(tibble)
library(Matrix)

# Load the FedGEE package
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
###############################################################################
###############################################################################
# C. SINGLE REPLICATION FUNCTION
#
# Input:  one row of parameter grid + rep ID
# Output: data.frame with one row per METHOD, columns for all metrics
###############################################################################
run_one_rep <- function(params, rep_id) {
  
  # Unpack parameters
  rho_H       <- params$rho_H
  rho_P       <- params$rho_P
  K           <- params$K
  pat_range   <- c(params$pat_lo, params$pat_hi)
  vis_range   <- c(params$vis_lo, params$vis_hi)
  incl_teach  <- params$include_teaching
  beta_teach  <- params$beta_teach
  true_beta_x <- params$true_beta_x
  
  # Add beta_0 with fallback to -1 if older grid is used
  beta_0      <- if (is.null(params$beta_0)) -1 else params$beta_0
  
  # Formula depends on whether teaching is included
  if (incl_teach) {
    form <- y ~ x + teaching
  } else {
    form <- y ~ x
  }
  
  # Generate data (now passing beta_0)
  d <- generate_sim_data(
    rho_H = rho_H, rho_P = rho_P, n_hospitals = K,
    patients_per_hosp_range = pat_range,
    visits_per_patient_range = vis_range,
    beta_0 = beta_0,                 
    true_beta_x = true_beta_x,       
    include_teaching = incl_teach,
    beta_teach = beta_teach
  )
  
  data_list <- split(d, d$hosp_id)
  n_total   <- nrow(d)
  n_patients <- n_distinct(d$pat_id)
  
  # ---- Shared metadata (Updated to include beta_0 and teaching) ----
  meta_row <- function(method, beta_hat, se_hat, beta_teach_hat = NA_real_, se_teach_hat = NA_real_, extra = list()) {
    out <- data.frame(
      # Parameters
      rho_H = rho_H, rho_P = rho_P, K = K,
      pat_lo = pat_range[1], pat_hi = pat_range[2],
      vis_lo = vis_range[1], vis_hi = vis_range[2],
      include_teaching = incl_teach,
      beta_0 = beta_0,
      rep = rep_id,
      n_total = n_total, n_patients = n_patients,
      # Method
      method = method,
      # Estimates for X
      beta_hat = beta_hat,
      se_hat = se_hat,
      # Estimates for Teaching
      beta_teach_hat = beta_teach_hat,
      se_teach_hat = se_teach_hat,
      # Derived metrics (for beta_x)
      bias = beta_hat - true_beta_x,
      abs_bias = abs(beta_hat - true_beta_x),
      # z-based CI for X
      ci_lo_z = beta_hat - 1.96 * se_hat,
      ci_hi_z = beta_hat + 1.96 * se_hat,
      covers_z = (beta_hat - 1.96 * se_hat <= true_beta_x) &
        (true_beta_x <= beta_hat + 1.96 * se_hat),
      ci_width_z = 2 * 1.96 * se_hat,
      stringsAsFactors = FALSE
    )
    # Append any extra columns (e.g., df, t-based coverage)
    for (nm in names(extra)) out[[nm]] <- extra[[nm]]
    out
  }
  
  results <- list()
  
  # ---- 1. Pooled GEE ----
  tryCatch({
    d$hosp_num <- as.numeric(factor(d$hosp_id))
    fit <- geeglm(form, id = hosp_num, data = d,
                  family = binomial, corstr = "independence")
    
    b_t  <- if(incl_teach) coef(fit)["teaching"] else NA_real_
    se_t <- if(incl_teach) sqrt(vcov(fit)["teaching", "teaching"]) else NA_real_
    
    results$pooled <- meta_row(
      "pooled_gee",
      coef(fit)["x"],
      sqrt(vcov(fit)["x", "x"]),
      beta_teach_hat = b_t,
      se_teach_hat = se_t,
      extra = list(converged = TRUE, iterations = NA_integer_)
    )
  }, error = function(e) {
    results$pooled <<- meta_row("pooled_gee", NA, NA, NA, NA,
                                extra = list(converged = FALSE, iterations = NA_integer_))
  })
  
  # ---- 2. FedGEE: site-level sandwich ----
  tryCatch({
    fit <- fedgee(data_list, form, binomial(), "independence", "pat_id",
                  sandwich_level = "site", md_correction = FALSE, verbose = FALSE)
    
    b_t  <- if(incl_teach) fit$coefficients["teaching"] else NA_real_
    se_t <- if(incl_teach) fit$se["teaching"] else NA_real_
    
    df_r <- fit$df_residual
    t_crit <- qt(0.975, df = max(df_r, 1))
    results$fed_site <- meta_row(
      "fed_site",
      fit$coefficients["x"],
      fit$se["x"],
      beta_teach_hat = b_t,
      se_teach_hat = se_t,
      extra = list(
        converged = fit$converged,
        iterations = fit$iterations,
        df = df_r,
        ci_lo_t = fit$coefficients["x"] - t_crit * fit$se["x"],
        ci_hi_t = fit$coefficients["x"] + t_crit * fit$se["x"],
        covers_t = (fit$coefficients["x"] - t_crit * fit$se["x"] <= true_beta_x) &
          (true_beta_x <= fit$coefficients["x"] + t_crit * fit$se["x"]),
        ci_width_t = 2 * t_crit * fit$se["x"]
      )
    )
  }, error = function(e) {
    results$fed_site <<- meta_row("fed_site", NA, NA, NA, NA,
                                  extra = list(converged = FALSE, iterations = NA_integer_))
  })
  
  # ---- 3. FedGEE: site-level + MD ----
  tryCatch({
    fit <- fedgee(data_list, form, binomial(), "independence", "pat_id",
                  sandwich_level = "site", md_correction = TRUE, verbose = FALSE)
    
    b_t  <- if(incl_teach) fit$coefficients["teaching"] else NA_real_
    se_t <- if(incl_teach) fit$se["teaching"] else NA_real_
    
    df_r <- fit$df_residual
    t_crit <- qt(0.975, df = max(df_r, 1))
    results$fed_md <- meta_row(
      "fed_site_md",
      fit$coefficients["x"],
      fit$se["x"],
      beta_teach_hat = b_t,
      se_teach_hat = se_t,
      extra = list(
        converged = fit$converged,
        iterations = fit$iterations,
        df = df_r,
        ci_lo_t = fit$coefficients["x"] - t_crit * fit$se["x"],
        ci_hi_t = fit$coefficients["x"] + t_crit * fit$se["x"],
        covers_t = (fit$coefficients["x"] - t_crit * fit$se["x"] <= true_beta_x) &
          (true_beta_x <= fit$coefficients["x"] + t_crit * fit$se["x"]),
        ci_width_t = 2 * t_crit * fit$se["x"]
      )
    )
  }, error = function(e) {
    results$fed_md <<- meta_row("fed_site_md", NA, NA, NA, NA,
                                extra = list(converged = FALSE, iterations = NA_integer_))
  })
  
  # ---- 4. FedGEE: patient-level sandwich ----
  tryCatch({
    fit <- fedgee(data_list, form, binomial(), "independence", "pat_id",
                  sandwich_level = "patient", md_correction = FALSE, verbose = FALSE)
    
    b_t  <- if(incl_teach) fit$coefficients["teaching"] else NA_real_
    se_t <- if(incl_teach) fit$se["teaching"] else NA_real_
    
    results$fed_pat <- meta_row(
      "fed_patient",
      fit$coefficients["x"],
      fit$se["x"],
      beta_teach_hat = b_t,
      se_teach_hat = se_t,
      extra = list(converged = fit$converged, iterations = fit$iterations)
    )
  }, error = function(e) {
    results$fed_pat <<- meta_row("fed_patient", NA, NA, NA, NA,
                                 extra = list(converged = FALSE, iterations = NA_integer_))
  })
  
  # ---- 5. Meta-analytic GLM ----
  tryCatch({
    m <- fit_meta_glm(d)
    results$meta <- meta_row(
      "meta_glm",
      m$beta, m$se,
      beta_teach_hat = NA_real_, # Meta GLM cannot estimate site-constant variables
      se_teach_hat = NA_real_,
      extra = list(converged = TRUE, iterations = NA_integer_)
    )
  }, error = function(e) {
    results$meta <<- meta_row("meta_glm", NA, NA, NA, NA,
                              extra = list(converged = FALSE, iterations = NA_integer_))
  })
  
  # ---- 6. GLMM ----
  tryCatch({
    fit <- glmer(form_glmm <- y ~ x + (1 | hosp_id) + (1 | pat_id),
                 data = d, family = binomial, nAGQ = 1,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
    
    se_glmm <- sqrt(diag(vcov(fit)))["x"]
    
    b_t  <- if(incl_teach) fixef(fit)["teaching"] else NA_real_
    se_t <- if(incl_teach) sqrt(diag(vcov(fit)))["teaching"] else NA_real_
    
    results$glmm <- meta_row(
      "glmm",
      fixef(fit)["x"],
      se_glmm,
      beta_teach_hat = b_t,
      se_teach_hat = se_t,
      extra = list(converged = TRUE, iterations = NA_integer_,
                   sigma2_hosp = VarCorr(fit)$hosp_id[1,1],
                   sigma2_pat  = VarCorr(fit)$pat_id[1,1])
    )
  }, error = function(e) {
    results$glmm <<- meta_row("glmm", NA, NA, NA, NA,
                              extra = list(converged = FALSE, iterations = NA_integer_))
  })
  
  # ---- Combine and return ----
  bind_rows(results)
}
###############################################################################
# D. PARAMETER GRID
#
# One row = one unique simulation scenario.
# Each scenario is repeated n_reps times.
###############################################################################
###############################################################################
# D. PARAMETER GRID
#
# One row = one unique simulation scenario.
# Each scenario is repeated n_reps times.
###############################################################################
build_param_grid <- function() {
  
  # Core grid: rho_H x rho_P x K
  core <- expand.grid(
    rho_H = c(0, 0.05, 0.1, 0.15, 0.20),
    rho_P = c(0.225, 0.25, 0.275, 0.300, 0.325, 0.350, 0.375, 0.40),
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