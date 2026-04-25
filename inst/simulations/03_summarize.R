###############################################################################
# 03_aggregate_results.R
#
# Combines all saved RDS files and produces tables/figures for each simulation.
# The parameter grid is rich enough that all 6 simulations are VIEWS
# of the same results — just different slices and groupings.
###############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(kableExtra)

out_dir <- "/ihome/spurkayastha/soumik/2026_fedGee/simulations/sim_results"
fig_dir <- "/ihome/spurkayastha/soumik/2026_fedGee/simulations/sim_figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load and combine all results ----
rds_files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
cat("Found", length(rds_files), "result files\n")

df <- bind_rows(lapply(rds_files, readRDS))
cat("Total rows:", nrow(df), "\n")
cat("Unique methods:", paste(unique(df$method), collapse = ", "), "\n")

# Remove error rows
df <- df |> filter(method != "ERROR")

# ---- 2. Summary function ----
summarise_metrics <- function(data, truth = 1) {
  data |>
    filter(is.finite(beta_hat), is.finite(se_hat), se_hat > 0) |>
    summarise(
      n_reps     = n(),
      n_converged = sum(converged, na.rm = TRUE),
      # Bias
      mean_est   = mean(beta_hat),
      mean_bias  = mean(beta_hat - truth),
      MAB        = mean(abs(beta_hat - truth)),
      median_bias = median(beta_hat - truth),
      # Variance
      emp_sd     = sd(beta_hat),
      mean_se    = mean(se_hat),
      se_sd_ratio = mean(se_hat) / sd(beta_hat),
      # MSE
      rmse       = sqrt(mean((beta_hat - truth)^2)),
      # Coverage (z-based)
      cov_z      = mean(covers_z, na.rm = TRUE),
      # Coverage (t-based, if available)
      cov_t      = if ("covers_t" %in% names(data)) mean(covers_t, na.rm = TRUE) else NA_real_,
      # CI width
      mean_width_z = mean(ci_width_z, na.rm = TRUE),
      mean_width_t = if ("ci_width_t" %in% names(data)) mean(ci_width_t, na.rm = TRUE) else NA_real_,
      # Power (reject H0: beta = 0 at 5%)
      power_z    = mean(abs(beta_hat / se_hat) > 1.96, na.rm = TRUE),
      .groups    = "drop"
    )
}

###############################################################################
# SIMULATION 1: Federated GEE = Pooled GEE, Meta-analysis degrades
#
# Slice: include_teaching = FALSE, multiple (rho_H, rho_P) settings
# Compare: pooled_gee, fed_site, meta_glm
###############################################################################
cat("\n=== SIMULATION 1: Fed vs Pooled vs Meta ===\n")

sim1 <- df |>
  filter(!include_teaching,
         method %in% c("fed_site", "fed_site_md", "meta_glm"),
         K == 5) |>
  group_by(rho_H, rho_P, method) |>
  summarise_metrics()

sim1_wide <- sim1 |>
  select(rho_H, rho_P, method, MAB, cov_z) |>
  pivot_wider(names_from = method, values_from = c(MAB, cov_z))

write.csv(sim1_wide, file.path(fig_dir, "sim1_table.csv"), row.names = FALSE)

# Figure
p_sim1 <- sim1 |>
  ggplot(aes(x = rho_P, y = cov_z, color = method)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  facet_wrap(~ paste0("rho_H = ", rho_H), nrow = 1) +
  scale_y_continuous(limits = c(0.7, 1), labels = scales::percent) +
  labs(title = "Sim 1: Coverage by Method",
       x = expression(rho[P]), y = "Coverage (z-based)", color = "Method") +
  theme_bw() + theme(legend.position = "bottom")

# sim1_k5 <- df |>
#   filter(!include_teaching, K == 5, rho_P == 0.6) |>
#   group_by(rho_H, method) |>
#   summarise(
#     cov_z = mean(covers_z, na.rm = TRUE),
#     cov_t = mean(covers_t, na.rm = TRUE),
#     mean_se = mean(se_hat, na.rm = TRUE),
#     emp_sd  = sd(beta_hat, na.rm = TRUE),
#     .groups = "drop"
#   ) |>
#   filter(method %in% c("fed_site", "fed_site_md"))


ggsave(file.path(fig_dir, "sim1_coverage.pdf"), p_sim1, width = 10, height = 6)

###############################################################################
# SIMULATION 2: GEE vs GLMM estimand divergence
#
# Slice: Total correlation grows (rho_P - rho_H ~ constant, but rho_P increases)
# Compare: fed_site vs glmm
###############################################################################
cat("\n=== SIMULATION 2: GEE vs GLMM ===\n")

# Find settings where rho_P - rho_H is approximately 0.3 (or closest)
sim2 <- df |>
  filter(!include_teaching,
         method %in% c("fed_site", "glmm"),
         K == 25,
         abs((rho_P - rho_H) - 0.3) < 0.1) |>
  group_by(rho_H, rho_P, method) |>
  summarise_metrics()

p_sim2 <- sim2 |>
  ggplot(aes(x = rho_H, y = mean_est, color = method, fill = method)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean_est - 1.96 * emp_sd / sqrt(n_reps),
                  ymax = mean_est + 1.96 * emp_sd / sqrt(n_reps)),
              alpha = 0.15, color = NA) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Sim 2: GEE vs GLMM (total correlation grows)",
       x = expression(rho[H]), y = "Mean Estimate", color = "Model") +
  theme_bw() + theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "sim2_estimand.pdf"), p_sim2, width = 7, height = 5)

###############################################################################
# SIMULATION 3: Site-level covariates
#
# Slice: include_teaching = TRUE
# Compare: fed_site vs pooled_gee (check teaching coefficient recovery)
# NOTE: need to also save beta_teach from run_one_rep — if not available,
#       this section can be run separately with a targeted simulation
###############################################################################
cat("\n=== SIMULATION 3: Site-level covariates ===\n")

sim3 <- df |>
  filter(include_teaching,
         method %in% c("pooled_gee", "fed_site"),
         K %in% c(25, 50)) |>
  group_by(rho_H, rho_P, K, method) |>
  summarise_metrics()

write.csv(sim3, file.path(fig_dir, "sim3_table.csv"), row.names = FALSE)
cat("Sim 3: See sim3_table.csv (teaching coefficient recovery from beta_x results)\n")
cat("Note: For teaching-specific coefficient, extend run_one_rep to also save beta_teach.\n")

###############################################################################
# SIMULATION 4: Sandwich level comparison
#
# Slice: Compare fed_site, fed_site_md, fed_patient across rho_H
# Key: patient-level sandwich undercovers when rho_H > 0
###############################################################################
cat("\n=== SIMULATION 4: Sandwich level comparison ===\n")

sim4 <- df |>
  filter(!include_teaching,
         method %in% c("fed_site", "fed_site_md", "fed_patient"),
         K == 25,
         rho_P == 0.6) |>
  group_by(rho_H, method) |>
  summarise_metrics()

p_sim4 <- sim4 |>
  ggplot(aes(x = rho_H, y = cov_z, color = method)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_y_continuous(limits = c(0.7, 1), labels = scales::percent) +
  labs(title = "Sim 4: Coverage by Sandwich Variant (K = 25)",
       x = expression(rho[H]), y = "Coverage", color = "Sandwich") +
  theme_bw() + theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "sim4_sandwich.pdf"), p_sim4, width = 7, height = 5)

# SE/SD ratio table
sim4_tbl <- sim4 |>
  select(rho_H, method, MAB, mean_se, emp_sd, se_sd_ratio, cov_z)
write.csv(sim4_tbl, file.path(fig_dir, "sim4_table.csv"), row.names = FALSE)

###############################################################################
# SIMULATION 5: Number of clusters + small-sample corrections
#
# Slice: Vary K, compare z-based, t-based, MD+t for fed_site / fed_site_md
###############################################################################
cat("\n=== SIMULATION 5: Cluster count and corrections ===\n")

sim5_z <- df |>
  filter(!include_teaching, method == "fed_site", rho_H == 0.2, rho_P == 0.6) |>
  group_by(K) |>
  summarise_metrics() |>
  mutate(correction = "z-based") |>
  select(K, correction, cov = cov_z, mean_width = mean_width_z, se_sd_ratio)

sim5_t <- df |>
  filter(!include_teaching, method == "fed_site", rho_H == 0.2, rho_P == 0.6) |>
  group_by(K) |>
  summarise(cov = mean(covers_t, na.rm = TRUE),
            mean_width = mean(ci_width_t, na.rm = TRUE),
            se_sd_ratio = mean(se_hat, na.rm = TRUE) / sd(beta_hat, na.rm = TRUE),
            .groups = "drop") |>
  mutate(correction = "t-based")

sim5_md <- df |>
  filter(!include_teaching, method == "fed_site_md", rho_H == 0.2, rho_P == 0.6) |>
  group_by(K) |>
  summarise(cov = mean(covers_t, na.rm = TRUE),
            mean_width = mean(ci_width_t, na.rm = TRUE),
            se_sd_ratio = mean(se_hat, na.rm = TRUE) / sd(beta_hat, na.rm = TRUE),
            .groups = "drop") |>
  mutate(correction = "MD + t")

sim5 <- bind_rows(sim5_z, sim5_t, sim5_md)

p_sim5 <- sim5 |>
  ggplot(aes(x = K, y = cov, color = correction)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_y_continuous(limits = c(0.7, 1), labels = scales::percent) +
  labs(title = "Sim 5: Coverage vs Number of Hospitals",
       x = "K (hospitals)", y = "Coverage", color = "Correction") +
  theme_bw() + theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "sim5_clusters.pdf"), p_sim5, width = 7, height = 5)

###############################################################################
# SIMULATION 6: Working correlation choice
#
# This requires exchangeable runs — if those were included in the grid,
# slice here. Otherwise, flag for a separate run.
###############################################################################
cat("\n=== SIMULATION 6: Working correlation choice ===\n")
cat("Note: Simulation 6 requires exchangeable working correlation runs.\n")
cat("If not in the grid, run separately with corstr = 'exchangeable'.\n")

###############################################################################
# BONUS: Comprehensive dashboard figure
###############################################################################
cat("\n=== GENERATING DASHBOARD ===\n")

# Heatmap: coverage by (rho_H, K) for fed_site at fixed rho_P
heatmap_data <- df |>
  filter(!include_teaching, method == "fed_site", rho_P == 0.6) |>
  group_by(rho_H, K) |>
  summarise(coverage = mean(covers_z, na.rm = TRUE), .groups = "drop")

p_heat <- heatmap_data |>
  ggplot(aes(x = factor(K), y = factor(rho_H), fill = coverage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", coverage * 100)), size = 3) +
  scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                       midpoint = 0.95, limits = c(0.7, 1),
                       labels = scales::percent) +
  labs(title = "Coverage Heatmap: FedGEE (site-level sandwich)",
       subtitle = expression(rho[P] * " = 0.6"),
       x = "K (hospitals)", y = expression(rho[H]), fill = "Coverage") +
  theme_minimal()

ggsave(file.path(fig_dir, "dashboard_heatmap.pdf"), p_heat, width = 8, height = 5)

# SE/SD ratio heatmap
ratio_data <- df |>
  filter(!include_teaching, method == "fed_site", rho_P == 0.6) |>
  group_by(rho_H, K) |>
  summarise(ratio = mean(se_hat, na.rm = TRUE) / sd(beta_hat),
            .groups = "drop")

p_ratio <- ratio_data |>
  ggplot(aes(x = factor(K), y = factor(rho_H), fill = ratio)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", ratio)), size = 3) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                       midpoint = 1, limits = c(0.5, 1.5)) +
  labs(title = "SE/SD Ratio: FedGEE (site-level sandwich)",
       subtitle = expression(rho[P] * " = 0.6, ratio < 1 = anti-conservative"),
       x = "K (hospitals)", y = expression(rho[H]), fill = "SE/SD") +
  theme_minimal()

ggsave(file.path(fig_dir, "dashboard_se_ratio.pdf"), p_ratio, width = 8, height = 5)

cat("\nAll figures saved to:", fig_dir, "\n")
cat("Done.\n")