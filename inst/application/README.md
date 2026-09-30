# VHA application (manuscript Section 7)

Produces Table 3 and Figure 1 of the main text and Web Table S5 of the supplement.

## Files

| File | Purpose |
|---|---|
| `visn_analysis.R` | Full analysis. Run inside the VA research environment (VINCI) |
| `make_supp_va_table.R` | Web Table S5 from `results/` |
| `visn_contrast.R` | Clinician vs non-clinician contrast (refit with clinician as reference; run in VA; pending) |
| `visn_forest.R` | Figure 1 (manuscript version without title or notes) and a standalone version |
| `results/2026_09_27_visn_estimates_long.csv` | All estimates: every variance correction × reference distribution × coefficient |
| `results/2026_09_27_visn_headline.csv` | One row per coefficient, main comparisons |
| `results/2026_09_27_visn_meta.csv` | Meta-analytic GEE estimates |
| `results/2026_09_27_visn_results_central.RDS` | The above plus model diagnostics and `sessionInfo()` |

`results/` holds coefficient-level output only. It contains no patient-level data.

## Provenance

The analysis ran on 2026-09-29 in the VA environment (R 4.5.2). `visn_analysis.R` is
byte-identical to the script that was run. It sources `01_fedgee_v3.R`, whose checksum
equals that of `../simulations/centralized/code/01_fedgee_v3.R`.

## Rerun with VA data access

1. Place `visn_analysis.R`, `../simulations/centralized/code/01_fedgee_v3.R`, and
   `../simulations/decentralized/code/01_decentralized_fedgee.R` in one folder with the
   analytic file `01_cleanData.dta` (source cohort: Anderson et al., *Ann Intern Med*
   2026;179:794–803).
2. Edit the `setwd()` line at the top of `visn_analysis.R`.
3. Run `Rscript visn_analysis.R`. Output is written to `2026_09_27_output/`.
4. Copy the output files into `results/` and run `make_supp_va_table.R` and
   `visn_forest.R` from this directory.

The model: marginal logistic GEE, working independence, $p = 12$, $K = 18$ VISNs,
29,041 hospitalizations of 18,900 Veterans. Recommended inference: score-space KC
correction with Bell–McCaffrey degrees of freedom. `inst/check_main_tables.R` confirms that
Table 3 in `main.tex` equals `results/`.
