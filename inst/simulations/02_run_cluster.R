###############################################################################
# 02_run_cluster.R
#
# Three modes:
#   count   — Print grid size and recommended array settings, then exit
#   <int>   — SLURM array mode: run one chunk of tasks
#   local N — Local parallel mode with N cores
#
# Usage:
#   Rscript 02_run_cluster.R count          # Check before submitting
#   sbatch run_sim.slurm                    # Submit (calls this with $SLURM_ARRAY_TASK_ID)
#   Rscript 02_run_cluster.R local 8        # Local test with 8 cores
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

library(dplyr)
library(purrr)

library(FedGEE)
# Use system.file to make this portable
source(system.file("simulations/00_simulation_config.R", package = "FedGEE", mustWork = FALSE))
# Fallback for local testing if package is not yet compiled
if (!exists("build_param_grid")) source("00_simulation_config.R")

out_dir <- "sim_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_reps  <- 200
tasks   <- build_task_table(n_reps = n_reps)
n_tasks <- nrow(tasks)
n_scenarios <- n_distinct(tasks$scenario_id)

###############################################################################
# MODE: COUNT
###############################################################################
if (length(args) >= 1 && args[1] == "count") {
  
  grid <- build_param_grid()
  
  cat("============================================\n")
  cat("FedGEE Simulation Grid Summary\n")
  cat("============================================\n")
  cat(sprintf("  Unique scenarios:    %d\n", n_scenarios))
  cat(sprintf("  Reps per scenario:   %d\n", n_reps))
  cat(sprintf("  Methods per rep:     6\n"))
  cat(sprintf("  Total task calls:    %d\n", n_tasks))
  cat(sprintf("  Total result rows:   ~%d\n", n_tasks * 6))
  cat("\n  Parameter ranges:\n")
  cat(sprintf("    rho_H: %s\n", paste(sort(unique(grid$rho_H)), collapse = ", ")))
  cat(sprintf("    rho_P: %s\n", paste(sort(unique(grid$rho_P)), collapse = ", ")))
  cat(sprintf("    K:     %s\n", paste(sort(unique(grid$K)), collapse = ", ")))
  cat(sprintf("    Teaching: %d with, %d without\n",
              sum(grid$include_teaching), sum(!grid$include_teaching)))
  cat("\n  Array size recommendations:\n")
  for (n_arr in c(200, 500, 1000)) {
    chunk <- ceiling(n_tasks / n_arr)
    cat(sprintf("    --array=1-%-4d  ->  %d tasks/job  ->  ~%.0f min/job\n",
                n_arr, chunk, chunk * 0.5))
  }
  cat("\n  Test first: Rscript 02_run_cluster.R 1\n")
  cat("============================================\n")
  quit(save = "no")
}

###############################################################################
# MODE: SLURM ARRAY
###############################################################################
if (length(args) >= 1 && args[1] != "local") {
  
  array_id <- as.integer(args[1])
  n_array  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_MAX",
                                    unset = Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1000")))
  
  chunk_size <- ceiling(n_tasks / n_array)
  task_start <- (array_id - 1) * chunk_size + 1
  task_end   <- min(array_id * chunk_size, n_tasks)
  
  if (task_start > n_tasks) {
    cat(sprintf("Array %d: no tasks (total %d, %d slots). Exiting.\n",
                array_id, n_tasks, n_array))
    quit(save = "no")
  }
  
  my_tasks <- tasks[task_start:task_end, ]
  cat(sprintf("Array %d/%d: tasks %d-%d (%d tasks)\n",
              array_id, n_array, task_start, task_end, nrow(my_tasks)))
  
  results  <- vector("list", nrow(my_tasks))
  t_start  <- proc.time()
  
  for (i in seq_len(nrow(my_tasks))) {
    params <- my_tasks[i, ]
    set.seed(params$scenario_id * 10000 + params$rep_id)
    
    t0 <- proc.time()
    results[[i]] <- tryCatch(
      run_one_rep(params, rep_id = params$rep_id),
      error = function(e) {
        data.frame(rho_H = params$rho_H, rho_P = params$rho_P, K = params$K,
                   pat_lo = params$pat_lo, pat_hi = params$pat_hi,
                   vis_lo = params$vis_lo, vis_hi = params$vis_hi,
                   include_teaching = params$include_teaching,
                   rep = params$rep_id, method = "ERROR",
                   beta_hat = NA, se_hat = NA,
                   error_msg = conditionMessage(e),
                   stringsAsFactors = FALSE)
      }
    )
    elapsed <- (proc.time() - t0)["elapsed"]
    
    if (i %% 10 == 0 || i == nrow(my_tasks)) {
      total_elapsed <- (proc.time() - t_start)["elapsed"]
      rate <- i / total_elapsed * 60
      eta  <- (nrow(my_tasks) - i) / max(rate, 0.01)
      cat(sprintf("  [%d/%d] %.1fs | %.1f tasks/min | ETA %.1f min\n",
                  i, nrow(my_tasks), elapsed, rate, eta))
    }
  }
  
  out <- bind_rows(results)
  out_file <- file.path(out_dir, sprintf("results_%04d.rds", array_id))
  saveRDS(out, out_file)
  
  total_time <- (proc.time() - t_start)["elapsed"]
  n_errors   <- sum(out$method == "ERROR", na.rm = TRUE)
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
        error = function(e) {
          data.frame(rho_H = params$rho_H, rho_P = params$rho_P, K = params$K,
                     rep = params$rep_id, method = "ERROR",
                     beta_hat = NA, se_hat = NA, stringsAsFactors = FALSE)
        }
      )
    }
    bind_rows(results)
  }
  
  library(parallel)
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, {
    library(geepack); library(lme4); library(dplyr); library(purrr)
    library(tibble); library(Matrix)
    library(FedGEE)
    source(system.file("simulations/00_simulation_config.R", package = "FedGEE", mustWork = FALSE))
    if (!exists("build_param_grid")) source("00_simulation_config.R")
  })
  
  t0 <- proc.time()
  all_results <- parLapplyLB(cl, scenario_groups, run_batch)
  stopCluster(cl)
  
  out <- bind_rows(all_results)
  out_file <- file.path(out_dir, "results_all.rds")
  saveRDS(out, out_file)
  cat(sprintf("Done: %d rows in %.1f min -> %s\n",
              nrow(out), (proc.time() - t0)["elapsed"] / 60, out_file))
}

###############################################################################
# NO ARGS
###############################################################################
if (length(args) == 0) {
  cat("Usage:\n")
  cat("  Rscript 02_run_cluster.R count          # Grid info\n")
  cat("  Rscript 02_run_cluster.R <array_id>     # SLURM task\n")
  cat("  Rscript 02_run_cluster.R local <cores>   # Local parallel\n")
  cat("  Rscript 02_run_cluster.R local 4 test    # Quick test\n")
}