###############################################################################
# summarize_v3.R -- tables for the v3 run (or the local pilot).
#
# Usage (from sim_v3/):
#   Rscript code/summarize_v3.R output/fedgee_v3_all.rds      # full run
#   Rscript code/summarize_v3.R output/pilot_smallK.rds       # local pilot
#
# Reports, for the visit-level coefficient x:
#   1. Coverage of nominal 95% CIs: variance (site/KC/MD/FG) x reference
#      distribution (z, t_K-1, t_K-p, t_nu, t_BM, t_FGdf), by K-p band and
#      site-size imbalance.
#   2. SE/SD (mean model SE over empirical SD of beta-hat) by variance.
#   3. Median degrees of freedom by reference and imbalance.
#   4. Stability: eigenvalue clipping, FG bound hits, failures.
###############################################################################
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
f <- if (length(args) >= 1) args[1] else "output/fedgee_v3_all.rds"
d <- readRDS(f)
cat("File:", f, "| rows:", nrow(d), "| errors:", sum(d$estimator == "ERROR", na.rm = TRUE), "\n")

keys <- c("K", "p_extra", "beta_0", "site_size_dist")
fed <- d |>
  filter(estimator == "fed", coef == "x", variant %in% c("site", "kc", "md", "fg")) |>
  mutate(p = p_extra + 2,
         band = ifelse(K - p < 20, "K-p 5-19", "K-p >=20"),
         sizes = factor(site_size_dist, c("uniform", "cv078", "cv100")),
         variant = factor(variant, c("site", "kc", "md", "fg")),
         inference = factor(inference, c("z", "t_k1", "t_kp", "t_nu", "t_bm", "t_fgdf")))

# ---- 1. coverage: scenario means first, then averaged (equal scenario weight)
scen <- fed |>
  group_by(across(all_of(keys)), band, sizes, variant, inference) |>
  summarise(cov = mean(covers, na.rm = TRUE), n = sum(!is.na(covers)),
            df_med = median(df, na.rm = TRUE), .groups = "drop")

cov_tab <- scen |>
  group_by(band, sizes, variant, inference) |>
  summarise(cov = mean(cov), .groups = "drop") |>
  tidyr::pivot_wider(names_from = inference, values_from = cov) |>
  arrange(band, sizes, variant)
cat("\n=== 1. Coverage of nominal 95% intervals (mean over scenarios) ===\n")
print(as.data.frame(mutate(cov_tab, across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

# ---- 2. SE/SD by variance (z rows carry one row per rep per variance)
sesd <- fed |>
  filter(inference == "z") |>
  group_by(across(all_of(keys)), band, sizes, variant) |>
  summarise(se_sd = mean(se, na.rm = TRUE) / sd(beta_hat, na.rm = TRUE), .groups = "drop") |>
  group_by(band, sizes, variant) |>
  summarise(se_sd = round(mean(se_sd), 3), .groups = "drop") |>
  tidyr::pivot_wider(names_from = variant, values_from = se_sd)
cat("\n=== 2. SE/SD (1 = unbiased SE) ===\n")
print(as.data.frame(sesd), row.names = FALSE)

# ---- 3. typical df by reference distribution (KC and FG)
dfs <- scen |>
  filter(inference %in% c("t_k1", "t_kp", "t_nu", "t_bm", "t_fgdf"), variant %in% c("kc", "fg")) |>
  group_by(K, sizes, variant, inference) |>
  summarise(df = round(median(df_med), 1), .groups = "drop") |>
  tidyr::pivot_wider(names_from = inference, values_from = df)
cat("\n=== 3. Median df (over scenarios with that K) ===\n")
print(as.data.frame(dfs), row.names = FALSE)

# ---- 4. stability
dg <- d |> filter(variant == "diag")
if (nrow(dg) > 0) {
  cat("\n=== 4. Stability (per replicate) ===\n")
  print(as.data.frame(dg |>
    group_by(K, site_size_dist) |>
    summarise(reps = n(),
              realised_cv = round(mean(site_size_cv, na.rm = TRUE), 2),
              max_lev = round(mean(max_leverage, na.rm = TRUE), 3),
              kc_clip_pct = round(100 * mean(kc_clipped > 0, na.rm = TRUE), 2),
              fg_cap_pct = round(100 * mean(fg_capped > 0, na.rm = TRUE), 2),
              .groups = "drop")), row.names = FALSE)
}
na_df <- fed |> filter(inference %in% c("t_bm", "t_fgdf")) |>
  summarise(bm_na = mean(is.na(df[inference == "t_bm"])), fg_na = mean(is.na(df[inference == "t_fgdf"])))
cat(sprintf("\nMissing df: BM %.4f%%, FG %.4f%%\n", 100 * na_df$bm_na, 100 * na_df$fg_na))
