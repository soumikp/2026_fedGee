###############################################################################
# Fed-GEE VHA application -- clinician vs non-clinician consultation contrast.
#
# RUN LOCATION: VA network / VINCI, same folder as visn_analysis.R
# (needs 01_cleanData.dta and 01_fedgee_v3.R). Writes coefficient-level
# output only.
#
# Question: do the two consultation types differ in their association with
# MAUD initiation? Target: beta_nonclinician - beta_clinician.
#
# Method: refit the SAME model with clinician consultation as the reference
# level of consult3. The non-clinician coefficient is then the contrast. A
# change of reference level is an invertible linear reparameterization of beta.
# The KC operator T_i = (I - B_i B^{-1})^{-1/2} and the Bell-McCaffrey vector
# a_i = T_i' B^{-1} c are equivariant under such a map, so the KC standard
# error and Bell-McCaffrey df of the new coefficient are exactly those of the
# contrast in the original fit. (FG is not equivariant: its factor is
# diagonal. FG rows are reported for information only.)
#
# Data preparation is copied from visn_analysis.R (function load_maud) so the
# model is identical apart from the reference level.
###############################################################################
setwd("P:/ORD_Anderson_202410002D/Experiments/FEDGEE paper")
rm(list = ls())

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

FEDGEE_SRC <- "01_fedgee_v3.R"
DATA_DTA   <- "01_cleanData.dta"
OUT_DIR    <- "2026_09_27_output"
dir.create(OUT_DIR, showWarnings = FALSE)

.fedgee_env <- new.env(parent = globalenv())
sys.source(FEDGEE_SRC, envir = .fedgee_env)
fedgee_v3 <- get("fedgee_v2", envir = .fedgee_env)
options(fedgee.df_coefs = NULL)

# ---- data: identical to load_maud() in visn_analysis.R ----------------------
df <- haven::read_dta(DATA_DTA) %>%
  mutate(facility = as.factor(sta6a_dis), patienticn = as.factor(patienticn),
         visn = as.factor(visnfy17)) %>%
  arrange(patienticn)
df$consult3 <- relevel(factor(df$hospital_addiction_consult), ref = "0")
df$psych    <- factor(as.integer(as.character(df$speciality_cat) == "1"), levels = c(0, 1))
df$oud      <- factor(df$oud_dx)
df$icu_any  <- factor(df$icu)
df$can_high <- factor(as.integer(as.character(df$CAN_quint) %in% c("4", "5")), levels = c(0, 1))
df$frail    <- factor(as.integer(as.character(df$frail_cat2) %in% c("3", "4", "5")), levels = c(0, 1))
df$ama      <- factor(as.integer(as.character(df$discharge_cat) == "2"), levels = c(0, 1))
df$age65    <- factor(as.integer(as.character(df$age_cat2) == "3"), levels = c(0, 1))
df$female   <- factor(as.integer(as.character(df$sex) == "1"), levels = c(0, 1))
df$auditc_low <- factor(as.integer(as.character(df$auditc_baseline_cat_d) %in% c("1", "2")),
                        levels = c(0, 1))

f_model <- new_during_aud_med ~ consult3 + psych + oud + icu_any + can_high + frail +
  ama + age65 + female + auditc_low

fit_one <- function(d) {
  d <- droplevels(d)
  fedgee_v3(data_list = split(d, d$visn), main_formula = f_model,
            family_obj = binomial("logit"), corstr = "independence",
            id_col = "patienticn", n_iter = 50, tol = 1e-8, verbose = FALSE)
}

# ---- original parameterization (reference: no consultation) -----------------
fit0 <- fit_one(df)
# ---- reparameterized (reference: clinician consultation) --------------------
df1 <- df
df1$consult3 <- relevel(factor(df1$hospital_addiction_consult), ref = "1")
fit1 <- fit_one(df1)

# Check: the reparameterized estimate equals the difference of the originals.
b0 <- fit0$coefficients
b1 <- fit1$coefficients
chk <- unname((b0["consult32"] - b0["consult31"]) - b1["consult32"])
cat(sprintf("Reparameterization check: |(b2 - b1) - b_contrast| = %.2e\n", abs(chk)))
stopifnot(abs(chk) < 1e-6)
# Check: KC standard error of the contrast from the original fit, sqrt(c' V c),
# equals the KC standard error of the reparameterized coefficient.
cvec <- setNames(rep(0, length(b0)), names(b0)); cvec[c("consult31", "consult32")] <- c(-1, 1)
se_direct <- sqrt(drop(t(cvec) %*% fit0$variants$kc$vcov %*% cvec))
se_chk <- abs(se_direct - fit1$variants$kc$se[["consult32"]])
cat(sprintf("KC SE check: |sqrt(c'Vc) - se_contrast| = %.2e\n", se_chk))
stopifnot(se_chk < 1e-8)

# ---- contrast rows: non-clinician vs clinician -------------------------------
term <- "consult32"   # in fit1 this is non-clinician vs clinician
j <- which(names(b1) == term)
rows <- list()
for (vn in c("site", "kc", "md", "fg")) {
  v <- fit1$variants[[vn]]
  refs <- list(z = Inf, t_k1 = v$df_k1, t_kp = v$df_kp, t_bm = unname(v$df_bm[j]))
  for (ref in names(refs)) {
    dfv <- refs[[ref]]
    se  <- unname(v$se[j])
    cm  <- if (is.finite(dfv)) qt(0.975, dfv) else qnorm(0.975)
    est <- b1[[term]]
    rows[[length(rows) + 1]] <- data.frame(
      contrast = "non-clinician vs clinician consultation",
      variant = vn, inference = ref, estimate = est, se = se, df = dfv,
      ratio_of_OR = exp(est), CI_low = exp(est - cm * se), CI_high = exp(est + cm * se),
      p_value = if (is.finite(dfv)) 2 * pt(-abs(est / se), dfv) else 2 * pnorm(-abs(est / se)),
      stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)
out$role <- ifelse(out$variant == "kc" & out$inference == "t_bm", "PRIMARY", "")
print(out, digits = 3, row.names = FALSE)

write.csv(out, file.path(OUT_DIR, "2026_09_27_visn_consult_contrast.csv"), row.names = FALSE)
saveRDS(list(contrast = out, reparam_check = chk, se_check = se_chk, run_time = Sys.time(),
             session = utils::sessionInfo()),
        file.path(OUT_DIR, "2026_09_27_visn_consult_contrast.RDS"))
message("saved contrast results to ", OUT_DIR)
