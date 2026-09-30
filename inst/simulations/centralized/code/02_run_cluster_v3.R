###############################################################################
# 02_run_cluster_v2.R
#
# Modes:
#   count        — Print grid summary
#   <int>        — SLURM array mode: run one chunk of tasks (array task ID)
#   local <N>    — Local parallel mode with N cores
#   local N test — Quick test on a small subset
#
# Environment variables:
#   FEDGEE_SIM_DIR    Directory with the .R scripts (default: cwd)
#   FEDGEE_N_REPS     Reps per scenario (default: 200)
#   FEDGEE_EXPS       Comma-separated experiments to run (default: "E1,E2,E3,E5,E7,E8")
#                     Examples: "E1" (just E1), "E1,E7" (E1 and E7 only)
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(purrr)
})

# ---- Locate scripts ----
SCRIPT_DIR <- Sys.getenv("FEDGEE_SIM_DIR", unset = getwd())
source(file.path(SCRIPT_DIR, "01_fedgee_v3.R"))
source(file.path(SCRIPT_DIR, "00_simulation_config_v3.R"))

# ---- v3: FG's own df needs saws. Fail fast rather than silently return NA. ----
if (!requireNamespace("saws", quietly = TRUE)) {
  stop("Package 'saws' is not installed. Run install.packages('saws') on the cluster first.")
}

# ---- Choose which experiments to run ----
exps_env <- Sys.getenv("FEDGEE_EXPS", unset = "EV3")
EXPS <- strsplit(exps_env, ",")[[1]]
EXPS <- trimws(EXPS)

# Results directory. Defaults to sim_results_v3/ beside the code, but the
# cluster driver sets FEDGEE_OUT_DIR so that code/, logs/ and output/ stay
# separate under the project root.
out_dir <- Sys.getenv("FEDGEE_OUT_DIR", unset = file.path(SCRIPT_DIR, "sim_results_v3"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_reps  <- as.integer(Sys.getenv("FEDGEE_N_REPS", unset = "200"))
tasks   <- build_task_table(n_reps = n_reps, experiments = EXPS)
n_tasks <- nrow(tasks)

###############################################################################
# MODE: COUNT
###############################################################################
if (length(args) >= 1 && args[1] == "count") {

  grid <- build_param_grid(experiments = EXPS)

  cat("============================================\n")
  cat("Fed-GEE Simulation v2 — Grid Summary\n")
  cat("============================================\n")
  cat(sprintf("  Experiments:        %s\n", paste(EXPS, collapse = ", ")))
  cat(sprintf("  Total scenarios:    %d\n", nrow(grid)))
  cat(sprintf("  Reps per scenario:  %d\n", n_reps))
  cat(sprintf("  Total task calls:   %d\n", n_tasks))

  cat("\n  Scenarios by experiment:\n")
  by_exp <- grid |> count(experiment, name = "scenarios")
  for (i in seq_len(nrow(by_exp))) {
    cat(sprintf("    %-6s %d scenarios\n", by_exp$experiment[i], by_exp$scenarios[i]))
  }

  cat("\n  Parameter ranges by experiment:\n")
  for (e in EXPS) {
    g <- grid[grid$experiment == e, ]
    if (nrow(g) == 0) next
    cat(sprintf("    [%s]\n", e))
    cat(sprintf("      K:              %s\n", paste(sort(unique(g$K)), collapse = ", ")))
    cat(sprintf("      p_extra:        %s\n", paste(sort(unique(g$p_extra)), collapse = ", ")))
    cat(sprintf("      rho_H:          %s\n", paste(sort(unique(g$rho_H)), collapse = ", ")))
    cat(sprintf("      rho_P:          %s\n", paste(sort(unique(g$rho_P)), collapse = ", ")))
    cat(sprintf("      site_size_dist: %s\n", paste(sort(unique(g$site_size_dist)), collapse = ", ")))
    cat(sprintf("      fit_glmm:       %s\n", paste(sort(unique(g$fit_glmm)), collapse = ", ")))
  }

  cat("\n  Array-size recommendations (assume ~5-30s/task):\n")
  for (n_arr in c(200, 500, 1000)) {
    chunk <- ceiling(n_tasks / n_arr)
    cat(sprintf("    --array=1-%-4d  %d tasks/job\n", n_arr, chunk))
  }
  cat("\n  Test:  Rscript 02_run_cluster_v2.R 1\n")
  cat("============================================\n")
  quit(save = "no")
}

###############################################################################
# Helper: minimal error-row fallback (keeps schema consistent for bind_rows)
###############################################################################
make_error_row <- function(params, rep_id, msg) {
  tibble(
    experiment = params$experiment, rep = rep_id,
    rho_H = params$rho_H, rho_P = params$rho_P, K = params$K,
    beta_0 = params$beta_0, p_extra = params$p_extra,
    site_size_dist = params$site_size_dist,
    include_teaching = params$include_teaching,
    n_total = NA_integer_, n_patients = NA_integer_,
    estimator = "ERROR", variant = "ERROR", inference = "ERROR",
    coef = "ERROR", true_value = NA_real_,
    beta_hat = NA_real_, se = NA_real_, df = NA_real_,
    bias = NA_real_, abs_bias = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_,
    ci_width = NA_real_, covers = NA,
    error_msg = msg
  )
}

###############################################################################
# v3: hard per-rep timeout via a forked child (ported from the calibration
# gap-fill rerun, 02b_rerun_calib.R). A replicate can hang inside compiled
# geepack code (update_beta), which setTimeLimit cannot interrupt; in the
# calibration run this cost 25 of 499 chunks their full 8h wall. mccollect()
# returns NULL on timeout; we then SIGKILL the child and emit an ERROR row.
# mc.set.seed = FALSE keeps the parent's RNG stream, so draws are identical to
# calling run_one_rep() inline after set.seed().
###############################################################################
rep_timeout <- as.numeric(Sys.getenv("FEDGEE_REP_TIMEOUT", unset = "120"))
run_rep_guarded <- function(params, timeout_s) {
  job <- parallel::mcparallel(
    tryCatch(run_one_rep(params, rep_id = params$rep_id),
             error = function(e) structure(list(msg = conditionMessage(e)),
                                           class = "fedgee_rep_error")),
    mc.set.seed = FALSE)
  res <- parallel::mccollect(job, wait = FALSE, timeout = timeout_s)
  if (is.null(res)) {
    tools::pskill(job$pid, tools::SIGKILL)
    parallel::mccollect(job, wait = FALSE)
    return(make_error_row(params, params$rep_id, sprintf("TIMEOUT_%gs", timeout_s)))
  }
  val <- res[[1]]
  if (inherits(val, "try-error"))
    return(make_error_row(params, params$rep_id, paste("child error:", as.character(val))))
  if (inherits(val, "fedgee_rep_error"))
    return(make_error_row(params, params$rep_id, val$msg))
  val
}

###############################################################################
# MODE: SLURM ARRAY
###############################################################################
if (length(args) >= 1 && args[1] != "local") {

  array_id <- as.integer(args[1])
  n_array  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_MAX",
                                    unset = Sys.getenv("SLURM_ARRAY_TASK_COUNT",
                                                       unset = "500")))

  chunk_size <- ceiling(n_tasks / n_array)
  task_start <- (array_id - 1) * chunk_size + 1
  task_end   <- min(array_id * chunk_size, n_tasks)

  if (task_start > n_tasks) {
    cat(sprintf("Array %d: no tasks. Exiting.\n", array_id))
    quit(save = "no")
  }

  my_tasks <- tasks[task_start:task_end, ]
  cat(sprintf("Array %d/%d: tasks %d-%d (%d tasks)\n",
              array_id, n_array, task_start, task_end, nrow(my_tasks)))

  results <- vector("list", nrow(my_tasks))
  t_start <- proc.time()

  for (i in seq_len(nrow(my_tasks))) {
    params <- my_tasks[i, ]
    set.seed(params$scenario_id * 10000 + params$rep_id)

    t0 <- proc.time()
    results[[i]] <- run_rep_guarded(params, rep_timeout)
    elapsed <- (proc.time() - t0)["elapsed"]

    if (i %% 10 == 0 || i == nrow(my_tasks)) {
      total_elapsed <- (proc.time() - t_start)["elapsed"]
      rate <- i / total_elapsed * 60
      eta  <- (nrow(my_tasks) - i) / max(rate, 0.01)
      cat(sprintf("  [%d/%d] %s K=%3d p_x=%2d  %.1fs | %.1f tasks/min | ETA %.1f min\n",
                  i, nrow(my_tasks), params$experiment, params$K, params$p_extra,
                  elapsed, rate, eta))
    }
  }

  out <- bind_rows(results)
  out_file <- file.path(out_dir, sprintf("results_%04d.rds", array_id))
  saveRDS(out, out_file)

  total_time <- (proc.time() - t_start)["elapsed"]
  n_errors <- sum(out$estimator == "ERROR", na.rm = TRUE)
  cat(sprintf("\nSaved %d rows -> %s (%.1f min, %d errors)\n",
              nrow(out), out_file, total_time / 60, n_errors))
}

###############################################################################
# MODE: LOCAL PARALLEL
###############################################################################
if (length(args) >= 1 && args[1] == "local") {

  n_cores <- if (length(args) >= 2) as.integer(args[2]) else max(parallel::detectCores() - 1, 1)

  if (length(args) >= 3 && args[3] == "test") {
    tasks <- tasks |> filter(scenario_id <= 5)
    cat(sprintf("TEST MODE: %d tasks only\n", nrow(tasks)))
  }

  cat(sprintf("Local: %d tasks on %d cores\n", nrow(tasks), n_cores))

  scenario_groups <- split(tasks, tasks$scenario_id)

  run_batch <- function(batch) {
    results <- vector("list", nrow(batch))
    for (i in seq_len(nrow(batch))) {
      params <- batch[i, ]
      set.seed(params$scenario_id * 10000 + params$rep_id)
      results[[i]] <- tryCatch(
        run_one_rep(params, rep_id = params$rep_id),
        error = function(e) make_error_row(params, params$rep_id, conditionMessage(e))
      )
    }
    bind_rows(results)
  }

  library(parallel)
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(geepack); library(lme4); library(dplyr); library(purrr)
      library(tibble); library(Matrix); library(saws)
    })
  })
  clusterExport(cl, c("SCRIPT_DIR", "make_error_row"))
  clusterEvalQ(cl, {
    source(file.path(SCRIPT_DIR, "01_fedgee_v3.R"))
    source(file.path(SCRIPT_DIR, "00_simulation_config_v3.R"))
  })

  t0 <- proc.time()
  all_results <- parLapplyLB(cl, scenario_groups, run_batch)
  stopCluster(cl)

  out <- bind_rows(all_results)
  out_file <- file.path(out_dir, "results_local.rds")
  saveRDS(out, out_file)
  cat(sprintf("Done: %d rows in %.1f min -> %s\n",
              nrow(out), (proc.time() - t0)["elapsed"] / 60, out_file))
}

if (length(args) == 0) {
  cat("Usage:\n")
  cat("  Rscript 02_run_cluster_v2.R count\n")
  cat("  Rscript 02_run_cluster_v2.R <array_id>\n")
  cat("  Rscript 02_run_cluster_v2.R local <cores>\n")
  cat("  Rscript 02_run_cluster_v2.R local 4 test\n")
  cat("\nEnvironment variables:\n")
  cat("  FEDGEE_SIM_DIR  Script directory (default: cwd)\n")
  cat("  FEDGEE_N_REPS   Reps per scenario (default: 200)\n")
  cat("  FEDGEE_EXPS     Comma-separated experiments (default: E1,E2,E3,E5,E7,E8)\n")
  cat("                  Examples: 'E1', 'E1,E7'\n")
}
