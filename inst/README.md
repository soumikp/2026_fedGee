# Reproduction materials for the *Statistics in Medicine* manuscript

**Fed-GEE: Federated Generalized Estimating Equations for Multisite Health Data, with an
Application to Veterans' Health.** Feier Chang and Soumik Purkayastha.

The repository root is the R package `FedGEE`. This `inst/` directory holds what is
needed to reproduce the manuscript: the LaTeX sources, the simulation code and results,
and the VHA application code and coefficient-level results. The patient-level VHA data
are not included (see [Data](#data)).

## Quick start

From the repository root, with R and a TeX distribution installed:

```bash
bash inst/reproduce.sh
```

This takes about 15 seconds. It rebuilds every generated table and figure from the
tracked results, compiles `manuscript/main.pdf` and `manuscript/supplement.pdf`, and
runs `check_main_tables.R`. That check confirms that the numbers typed into Tables 1–3
of `main.tex` equal the script output. The script stops with an error if any number
differs.

## Directory layout

```
inst/
├── README.md                 this file
├── reproduce.sh              rebuild all tables, figures, and PDFs
├── check_main_tables.R       check main-text Tables 1-3 against script output
├── manuscript/               submission: main.tex, supplement.tex, references.bib, PDFs
│   ├── tables/               generated LaTeX tables (do not edit by hand)
│   └── figures/              generated figures
├── simulations/
│   ├── centralized/          simulation study of Section 6 (Tables 1-2, Web Tables S3-S4)
│   └── decentralized/        decentralized simulation (Web Appendix I, Web Tables S1-S2)
└── application/              VHA application (Section 7, Table 3, Figure 1, Web Table S5)
```

Each subdirectory has its own README with details. `inst/_archive/` holds earlier work.
It is kept locally and is not part of the repository.

## Map from manuscript to code

| Manuscript item | Produced by | Input (tracked) | Output |
|---|---|---|---|
| Table 1, Table 2 | `simulations/centralized/code/make_statmed_tables_v2.R` | `simulations/centralized/results/statmed_scenario_table.rds` | rows printed to `results/main_table_rows.txt`, typed into `main.tex` |
| Web Tables S3, S4 | same script | same | `manuscript/tables/supp_sim_tables.tex` |
| Section 6 text (ranges, counts) | same script | same | printed to `results/main_table_rows.txt` |
| Web Tables S1, S2 | `simulations/decentralized/code/make_dec_tables.R` | `simulations/decentralized/results/dec_v4_summary.rds` | `manuscript/tables/supp_dec_tables.tex` |
| Section 5 ($\rho$ for 18 VISNs) | `simulations/decentralized/code/01_decentralized_fedgee.R` (`build_weight_matrix`) | none (network only) | see `simulations/decentralized/README.md` |
| Table 3 | `application/visn_analysis.R` (run inside VA) | VHA data (not shared) | `application/results/*.csv`, typed into `main.tex` |
| Figure 1 | `application/visn_forest.R` | `application/results/` | `manuscript/figures/visn_forest_paper.pdf` |
| Web Table S5 | `application/make_supp_va_table.R` | `application/results/2026_09_27_visn_estimates_long.csv` | `manuscript/tables/supp_va_table.tex` |
| Estimators (all) | `simulations/centralized/code/01_fedgee_v3.R`; package `R/fedgee.R` | | |

## Three levels of reproduction

1. **Tables and figures from tracked results** (`reproduce.sh`, seconds). Needs only the
   files in this repository.
2. **Per-scenario summaries from raw simulation output.** The table scripts recompute the
   summaries automatically when the raw output is present in `simulations/*/output/`.
   The raw output (113 MB and 160,000 rows) is not tracked because of its size. It can be
   regenerated with level 3.
3. **Raw simulation output from scratch.** Run the SLURM drivers described in
   `simulations/centralized/README.md` (117,000 tasks, about one hour on 499 array jobs)
   and `simulations/decentralized/README.md`.

## Software

Versions used for the centralized simulation and the tables (from
`simulations/centralized/results/RUN_PROVENANCE.txt` and the local build):

- R 4.4.1; dplyr 1.2.0, tidyr 1.3.2, purrr 1.2.1, tibble 3.3.1, Matrix 1.7.4,
  geepack 1.3.13, saws 0.9.7.0, clubSandwich 0.6.1, ggplot2 4.0.1, gridExtra 2.3, haven 2.5.5.
- The VHA analysis ran on R 4.5.2 inside the VA; its full `sessionInfo()` is stored in
  `application/results/2026_09_27_visn_results_central.RDS` (element `session`).
- LaTeX: pdflatex and BibTeX with `natbib`, `booktabs`, `longtable`, `lmodern`.

## Provenance checks

- The MD5 checksums of `01_fedgee_v3.R`, `00_simulation_config_v3.R`, and
  `02_run_cluster_v3.R` in `simulations/centralized/code/` equal those recorded in
  `RUN_PROVENANCE.txt` when the cluster results were collected. The tracked code is the
  code that produced the results.
- The VHA analysis used the same `01_fedgee_v3.R` (same checksum) and the script
  `application/visn_analysis.R` (byte-identical to the file run in the VA).
- The decentralized simulation (Web Appendix I) was run in August 2026 with an earlier
  version of the `FedGEE` package for its centralized comparator. That package version
  was not recorded. See `simulations/decentralized/README.md`.

## Data

The VHA data are protected health information and cannot be shared. Access requires VA
research approval. `application/results/` holds only coefficient-level output (odds
ratios, standard errors, degrees of freedom, confidence intervals, and model
diagnostics); every number in it appears in the manuscript or its supplement.
`application/README.md` describes how to rerun the analysis with VA data access.
