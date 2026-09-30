# Check that the numbers typed into Tables 1-3 of manuscript/main.tex equal the
# script output. Run from inst/ after reproduce.sh:   Rscript check_main_tables.R
# Exits with an error if any row differs.

tex <- readLines("manuscript/main.tex")
table_rows <- function(label) {
  i <- grep(sprintf("\\\\label\\{%s\\}", label), tex)
  a <- i + grep("^\\\\midrule", tex[i:length(tex)])[1]
  b <- i + grep("^\\\\bottomrule", tex[i:length(tex)])[1] - 2
  x <- tex[a:b]
  x <- trimws(x)
  x[grepl("[0-9]\\.[0-9]", x) & !startsWith(x, "\\multicolumn")]  # data rows only
}
norm <- function(x) gsub("\\s+", " ", trimws(x))
ok <- TRUE
report <- function(name, want, got) {
  bad <- setdiff(norm(want), norm(got))
  if (length(bad) || length(want) != length(got)) {
    ok <<- FALSE
    cat(name, ": MISMATCH\n"); print(setdiff(norm(want), norm(got))); print(setdiff(norm(got), norm(want)))
  } else cat(name, ": OK (", length(want), "rows )\n")
}

# Tables 1-2: rows printed by simulations/centralized/code/make_statmed_tables_v2.R
out <- readLines("simulations/centralized/results/main_table_rows.txt")
t1 <- out[(grep("MAIN Table 1", out) + 1):(grep("MAIN Table 2", out) - 1)]
t2 <- out[(grep("MAIN Table 2", out) + 1):(grep("^%% ranges", out) - 2)]
t2 <- t2[nzchar(trimws(t2))]
report("Table 1", t1, table_rows("tab:simulation"))
report("Table 2", t2, table_rows("tab:comparator"))

# Table 3: application estimates
e <- read.csv("application/results/2026_09_27_visn_estimates_long.csv")
m <- read.csv("application/results/2026_09_27_visn_meta.csv")
g <- function(v, i) e[e$variant == v & e$inference == i, ]
P <- g("site", "z"); F <- g("kc", "t_bm")
lab <- c(female1 = "\\quad Female", age651 = "\\quad Age $\\ge$65 years",
  auditc_low1 = "\\quad AUDIT-C no or low use", oud1 = "\\quad Opioid use disorder",
  can_high1 = "\\quad CAN score $\\ge$60th percentile",
  frail1 = "\\quad Frailty, mild or worse", icu_any1 = "\\quad Any ICU stay",
  psych1 = "\\quad Psychiatry discharging service",
  consult31 = "\\qquad Clinician", consult32 = "\\qquad Non-clinician",
  ama1 = "\\quad Discharge against advice")
t3 <- vapply(names(lab), function(t) { p <- P[P$term == t, ]; f <- F[F$term == t, ]; mm <- m[m$term == t, ]
  sprintf("%s & %.2f & %.2f, %.2f & %.2f, %.2f & %.1f & %.2f (%.2f, %.2f)\\\\", lab[[t]], f$OR,
          p$CI_low, p$CI_high, f$CI_low, f$CI_high, f$df, mm$OR, mm$CI_low, mm$CI_high) }, character(1))
report("Table 3", unname(t3), table_rows("tab:application"))

if (!ok) stop("main.tex tables do not match the script output")
