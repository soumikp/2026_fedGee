#!/usr/bin/env bash
# Rebuild every table, figure, and PDF of the Statistics in Medicine manuscript
# from the tracked results in this directory. Runs in about two minutes.
#
#   bash inst/reproduce.sh          (from the repository root)
#
# Steps (see inst/README.md for what each output is):
#   1. Centralized simulation tables  -> manuscript/tables/supp_sim_tables.tex
#                                        simulations/centralized/results/main_table_rows.txt
#   2. Decentralized simulation tables -> manuscript/tables/supp_dec_tables.tex
#   3. VHA sensitivity table           -> manuscript/tables/supp_va_table.tex
#   4. VHA forest plot                 -> manuscript/figures/visn_forest_*.pdf/.png
#   5. LaTeX                           -> manuscript/main.pdf, manuscript/supplement.pdf
# Tables 1-3 of the main text are typed into main.tex. Step 1 prints the rows of
# Tables 1-2; results/2026_09_27_visn_estimates_long.csv holds Table 3.
set -euo pipefail
INST="$(cd "$(dirname "$0")" && pwd)"

echo "[1/5] centralized simulation tables"
( cd "$INST/simulations/centralized" && Rscript code/make_statmed_tables_v2.R > results/main_table_rows.txt )

echo "[2/5] decentralized simulation tables"
( cd "$INST/simulations/decentralized" && Rscript code/make_dec_tables.R > results/dec_table_log.txt )

echo "[3/5] VHA sensitivity table"
( cd "$INST/application" && Rscript make_supp_va_table.R > /dev/null )

echo "[4/5] VHA forest plot"
( cd "$INST/application" && Rscript visn_forest.R > /dev/null )

echo "[5/5] LaTeX"
( cd "$INST/manuscript"
  for f in main supplement; do
    pdflatex -interaction=nonstopmode "$f.tex" > /dev/null
    bibtex "$f" > /dev/null
    pdflatex -interaction=nonstopmode "$f.tex" > /dev/null
    pdflatex -interaction=nonstopmode "$f.tex" > /dev/null
    if grep -qE "^!|undefined" "$f.log"; then echo "  $f.tex: LaTeX errors or undefined references, see $f.log"; exit 1; fi
  done )

echo "Build done."

echo "[check] main-text Tables 1-3 against script output"
( cd "$INST" && Rscript check_main_tables.R )
