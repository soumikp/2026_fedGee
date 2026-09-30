# Decentralized simulation (Web Appendix I; manuscript Section 5)

Produces Web Tables S1–S2 of the supplement and the spectral quantities quoted in
Section 5 of the main text.

## Scope

This simulation was run in August 2026, before the correction and reference distribution
of the manuscript were finalized. Every interval uses the **uncorrected** site sandwich
with a $t_{K-p}$ reference. It shows how consensus error affects estimation and
inference. It does not evaluate the KC correction or the Bell–McCaffrey degrees of
freedom under consensus; the manuscript states this.

## Design

Latent-Gaussian binary outcomes as in the centralized study, with site-level correlation
0, 0.1, or 0.2, prevalence 0.07 or 0.30, $K \in \{5, 10, 15, 25, 50, 100\}$, and some
scenarios with a binary site-level covariate ($p = 3$). 48 scenarios with 200–400
replicates each. Four networks with Metropolis–Hastings weights: complete; regional
(groups of five sites around a hub, hubs fully connected); ring (each site linked to two
neighbors on each side); single hub. $L = 50$ consensus rounds per stage.
Each replicate sets `set.seed(scenario_id * 10000 + rep_id)`.

## Files

| File | Purpose |
|---|---|
| `code/01_decentralized_fedgee.R` | Decentralized estimator; `build_weight_matrix()`, `mh_weights()` |
| `code/00_simulation_config.R` | Data generation, one replicate, grid. Loads `library(FedGEE)` for the centralized comparator |
| `code/02_run_cluster.R`, `code/02_run_cluster.slurm` | Driver and SLURM submission |
| `code/make_dec_tables.R` | Web Tables S1–S2 (per scenario, no averaging) |
| `results/dec_v4_summary.rds` | $\rho$ by network and $K$; per-scenario coverage, SE/SD, bias/SD, replicate counts |
| `output/sim_results_v4/` | Raw output, 297 files, 160,000 rows (not tracked) |

## Rebuild the tables

```bash
cd inst/simulations/decentralized
Rscript code/make_dec_tables.R
```

The summary is recomputed from `output/sim_results_v4/` when present, otherwise read from
`results/dec_v4_summary.rds`. The spectral quantities are always recomputed from the
network definitions.

## Section 5 numbers for 18 VISNs

```r
source("code/01_decentralized_fedgee.R")
rho <- function(W) sort(abs(eigen(W, symmetric = TRUE, only.values = TRUE)$values), decreasing = TRUE)[2]
rho(build_weight_matrix(18, structure = "hub", hub = 1)$W)       # 0.944 -> 26 rounds for rho^L < 1/sqrt(18)
rho(build_weight_matrix(18, structure = "ring", K_neigh = 2)$W)  # 0.960 -> 36 rounds
```

## Known limitations for an audit

- `02_run_cluster.R` sources files by absolute paths on the original cluster account.
  Edit the two `source()` lines before rerunning.
- The centralized comparator came from the `FedGEE` package installed on the cluster in
  August 2026, a version before 0.2.0. That version was not recorded, so a rerun with the
  current package may differ in the centralized column. The decentralized columns depend
  only on `code/01_decentralized_fedgee.R`.
