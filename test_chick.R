library(geepack)
library(dplyr)
library(purrr)
library(tibble)
library(Matrix)

source("R/fedgee.R")

data(ChickWeight)
df <- ChickWeight
df$site <- paste0("Site_", as.numeric(df$Chick) %% 3 + 1)
data_list <- split(df, df$site)

fed_model <- fedgee(
  data_list      = data_list,
  main_formula   = weight ~ Time + Diet,
  family_obj     = gaussian(link = "identity"),
  corstr         = "exchangeable",
  id_col         = "Chick",
  sandwich_level = "site",
  correction     = "KC",
  verbose        = TRUE
)
