# FedGEE

**FedGEE** (Federated Generalized Estimating Equations) provides a robust framework for fitting generalized estimating equations across distributed datasets (e.g., multi-center health networks) without pooling patient-level data. It includes small-sample variance corrections (Kauermann-Carroll and Mancl-DeRouen) which are applied directly in the score space to preserve privacy and statistical validity when the number of contributing sites is small.

## Installation

You can install the development version of FedGEE from source using the `devtools` package:

```r
# install.packages("devtools")
# If installed locally:
devtools::install("path/to/FedGEE")

# Or if available on GitHub:
# devtools::install_github("username/FedGEE")
```

## Features
- **Privacy-Preserving**: Only summary-level matrices (e.g., Bread and Score matrices) are aggregated at the server.
- **Small-Sample Corrections**: Includes `KC` (Kauermann-Carroll) and `MD` (Mancl-DeRouen) corrections.
- **Cluster/Patient-Level Adjustments**: Allows sandwich variance estimation at the site-level or patient-level.

## Toy Example

The following example demonstrates how to deploy FedGEE across a list of site-specific data frames. This mimics the setting where `fedgee()` is provided with datasets from multiple distinct hospitals or servers.

```r
library(FedGEE)
set.seed(42)

# 1. Simulate data for 3 distinct sites
simulate_site_data <- function(site_id, n_patients_per_site = 50, records_per_patient = 4) {
  N <- n_patients_per_site * records_per_patient
  
  # Patient IDs
  pat_id <- rep(paste0("Site", site_id, "_P", 1:n_patients_per_site), each = records_per_patient)
  
  # Patient and Observation Covariates
  age <- rep(rnorm(n_patients_per_site, mean = 50, sd = 10), each = records_per_patient)
  treatment <- rbinom(N, 1, prob = 0.5)
  
  # Simulated outcome (Binary)
  eta <- -1.0 + 0.05 * age + 0.8 * treatment
  prob <- exp(eta) / (1 + exp(eta))
  outcome <- rbinom(N, 1, prob)
  
  data.frame(
    pat_id = pat_id, 
    age = age, 
    treatment = treatment, 
    outcome = outcome
  )
}

# Create a list where each element represents data sitting at a distinct site
data_list <- list(
  site1 = simulate_site_data(1),
  site2 = simulate_site_data(2),
  site3 = simulate_site_data(3)
)

# 2. Fit the Federated GEE Model
fed_model <- fedgee(
  data_list      = data_list,
  main_formula   = outcome ~ age + treatment,
  family_obj     = binomial(link = "logit"),
  corstr         = "exchangeable",     # Working correlation structure
  id_col         = "pat_id",           # Column indicating clusters (patients)
  sandwich_level = "site",             # Level for sandwich variance
  correction     = "KC",               # Kauermann-Carroll small-sample correction
  verbose        = TRUE
)

# 3. View the results
print(fed_model)
```
