###############################################################################
# collect_results.R
#
# Combines every results_XXXX.rds in output/ into ONE tidy file, runs
# completeness and sanity checks, and writes a provenance record so the
# manuscript can cite a single run.
#
# Produces, in output/:
#   fedgee_v3_all.rds        combined long-format results
#   fedgee_v3_all.csv.gz     same, portable
#   RUN_PROVENANCE.txt       what was run, when, with which code
#
# Usage (from the project root):
#   Rscript code/collect_results.R
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

PROJECT_ROOT <- Sys.getenv("FEDGEE_PROJECT_ROOT", unset = getwd())
out_dir  <- file.path(PROJECT_ROOT, "output")
code_dir <- file.path(PROJECT_ROOT, "code")

files <- list.files(out_dir, pattern = "^results_[0-9]+\\.rds$", full.names = TRUE)
cat("Found", length(files), "result files in", out_dir, "\n")
if (length(files) == 0) stop("No result files. Has the array finished?")

d <- bind_rows(lapply(files, readRDS))
cat("Combined rows:", nrow(d), "\n\n")

###############################################################################
# Completeness
###############################################################################
cat("=== COMPLETENESS ===\n")
if (all(c("experiment", "rep") %in% names(d))) {
  by_exp <- d |>
    group_by(experiment) |>
    summarise(rows = n(), reps = n_distinct(rep), .groups = "drop")
  print(as.data.frame(by_exp))

  # Flag scenarios with fewer reps than the maximum observed
  key <- intersect(c("experiment", "K", "rho_H", "p_extra", "beta_0",
                     "site_size_dist"),
                   names(d))
  per_scen <- d |>
    group_by(across(all_of(key))) |>
    summarise(reps = n_distinct(rep), .groups = "drop")
  target <- max(per_scen$reps)
  short  <- per_scen |> filter(reps < target)
  cat("\n  Target reps per scenario:", target, "\n")
  if (nrow(short) > 0) {
    cat("  INCOMPLETE scenarios:", nrow(short), "\n")
    print(as.data.frame(head(short, 20)))
  } else {
    cat("  All scenarios complete.\n")
  }
}

###############################################################################
# Sanity checks
###############################################################################
cat("\n=== SANITY CHECKS ===\n")
chk <- function(lbl, ok, extra = "") {
  cat(sprintf("  [%s] %-42s %s\n", if (ok) "PASS" else "WARN", lbl, extra))
}

if ("converged" %in% names(d)) {
  cr <- mean(d$converged, na.rm = TRUE)
  chk("convergence rate > 0.95", isTRUE(cr > 0.95), sprintf("%.4f", cr))
}
if ("se" %in% names(d)) {
  fin <- mean(is.finite(d$se[d$variant != "diag"]), na.rm = TRUE)
  chk("finite SEs (non-diag rows)", isTRUE(fin > 0.95), sprintf("%.4f", fin))
}
if (all(c("variant", "beta_hat", "rep") %in% names(d))) {
  # beta must NOT depend on the variance correction
  spread <- d |>
    filter(estimator == "fed",
           !variant %in% c("diag", "patient_z"), !is.na(beta_hat)) |>
    group_by(across(any_of(c("experiment", "K", "rho_H", "p_extra", "beta_0",
                             "site_size_dist", "rep", "coef")))) |>
    summarise(rng = diff(range(beta_hat)), .groups = "drop")
  mx <- max(spread$rng, na.rm = TRUE)
  chk("beta invariant across corrections", isTRUE(mx < 1e-8),
      sprintf("max spread %.2e", mx))
}
if ("kc_clipped" %in% names(d)) {
  kc <- sum(d$kc_clipped, na.rm = TRUE)
  md <- sum(d$md_clipped, na.rm = TRUE)
  cat(sprintf("  [INFO] eigenvalue clipping: kc=%d md=%d\n", kc, md))
}

###############################################################################
# Write combined output + provenance
###############################################################################
saveRDS(d, file.path(out_dir, "fedgee_v3_all.rds"))
utils::write.csv(d, gzfile(file.path(out_dir, "fedgee_v3_all.csv.gz")),
                 row.names = FALSE)

md5 <- function(f) {
  if (!file.exists(f)) return(NA_character_)
  unname(tools::md5sum(f))
}
prov <- c(
  "Fed-GEE v2 simulation -- run provenance",
  "=======================================",
  paste("Collected      :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Project root   :", PROJECT_ROOT),
  paste("Result files   :", length(files)),
  paste("Combined rows  :", nrow(d)),
  paste("R version      :", R.version.string),
  "",
  "Code checksums (md5):",
  paste("  01_fedgee_v3.R           ", md5(file.path(code_dir, "01_fedgee_v3.R"))),
  paste("  00_simulation_config_v3.R   ", md5(file.path(code_dir, "00_simulation_config_v3.R"))),
  paste("  02_run_cluster_v3.R      ", md5(file.path(code_dir, "02_run_cluster_v3.R"))),
  "",
  "Package versions:",
  paste("  ", sapply(c("geepack", "lme4", "dplyr", "purrr", "tibble", "Matrix"),
                     function(p) tryCatch(paste0(p, " ", as.character(packageVersion(p))),
                                          error = function(e) paste0(p, " NOT INSTALLED"))))
)
writeLines(prov, file.path(out_dir, "RUN_PROVENANCE.txt"))

cat("\nWrote:\n")
cat("  ", file.path(out_dir, "fedgee_v3_all.rds"), "\n")
cat("  ", file.path(out_dir, "fedgee_v3_all.csv.gz"), "\n")
cat("  ", file.path(out_dir, "RUN_PROVENANCE.txt"), "\n")
