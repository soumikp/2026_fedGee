# Load dependencies directly since package might not be installed in system path for this script
library(geepack)
library(dplyr)
library(ggplot2)
library(purrr)

source("R/fedgee.R")

data(ChickWeight)
df <- ChickWeight
df$site <- paste0("Site_", as.numeric(df$Chick) %% 3 + 1)
data_list <- split(df, df$site)

# 1. FedGEE Model
fed_model <- fedgee(
  data_list      = data_list,
  main_formula   = weight ~ Time + Diet,
  family_obj     = gaussian(link = "identity"),
  corstr         = "exchangeable",     
  id_col         = "Chick",            
  sandwich_level = "site",             
  correction     = "KC",               
  verbose        = FALSE
)

# 2. Pooled Model
pooled_model <- geeglm(
  weight ~ Time + Diet, 
  data = df[order(df$Chick), ], 
  id = Chick, 
  family = gaussian(link = "identity"), 
  corstr = "exchangeable"
)

# Extract FedGEE coefficients and standard errors
fed_est <- fed_model$coefficients
fed_se <- fed_model$se
fed_df <- data.frame(
  term = names(fed_est),
  estimate = fed_est,
  std.error = fed_se,
  model = "FedGEE"
)

# Extract Pooled GEE coefficients and standard errors
pooled_sum <- summary(pooled_model)$coefficients
pooled_df <- data.frame(
  term = rownames(pooled_sum),
  estimate = pooled_sum$Estimate,
  std.error = pooled_sum$Std.err,
  model = "Pooled-GEE"
)

# Combine and calculate 95% CI
plot_data <- bind_rows(fed_df, pooled_df) |>
  mutate(
    conf.low = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error
  )

# Plot
p <- ggplot(plot_data, aes(x = estimate, y = term, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), 
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(
    title = "Comparison of FedGEE vs Pooled GEE",
    x = "Estimate (95% CI)",
    y = "Coefficient",
    color = "Model"
  )

ggsave("fedgee_comparison.png", p, width = 7, height = 5)
