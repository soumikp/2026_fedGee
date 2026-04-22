# FedGEE

**FedGEE** (Federated Generalized Estimating Equations) provides a robust framework for fitting generalized estimating equations across distributed datasets (e.g., multi-center health networks) without pooling patient-level data. It includes small-sample variance corrections (Kauermann-Carroll and Mancl-DeRouen) which are applied directly in the score space to preserve privacy and statistical validity when the number of contributing sites is small.

## Installation

You can install the development version of FedGEE from source using the `devtools` package:

```r
# install from source available on GitHub:
# devtools::install_github("soumikp/2026_fedGee")
```

## Features
- **Privacy-Preserving**: Only summary-level matrices (e.g., Bread and Score matrices) are aggregated at the server.
- **Small-Sample Corrections**: Includes `KC` (Kauermann-Carroll) and `MD` (Mancl-DeRouen) corrections.
- **Cluster/Patient-Level Adjustments**: Allows sandwich variance estimation at the site-level or patient-level.

## Toy Example

Here is a small example using a known built-in R dataset (`ChickWeight`), simulating a scenario where longitudinal patient records are distributed across multiple "sites". 

```r
library(FedGEE)

# 1. Prepare data using the built-in 'ChickWeight' dataset
# We'll split the data into 3 arbitrary "sites" based on the Chick ID modulo 3
data(ChickWeight)
df <- ChickWeight

# Create a mock 'site' variable
df$site <- paste0("Site_", as.numeric(df$Chick) %% 3 + 1)

# Split the dataframe into a list of dataframes by site
data_list <- split(df, df$site)

# 2. Fit the Federated GEE Model
# We are modeling weight ~ Time + Diet
fed_model <- fedgee(
  data_list      = data_list,
  main_formula   = weight ~ Time + Diet,
  family_obj     = gaussian(link = "identity"),
  corstr         = "exchangeable",     # Working correlation structure
  id_col         = "Chick",            # Column indicating clusters
  sandwich_level = "site",             # Level for sandwich variance
  correction     = "KC",               # Kauermann-Carroll small-sample correction
  verbose        = TRUE
)

# 3. View the results
print(fed_model)
```
