###############################################################################
# Fed-GEE — VHA application, VISN as the unit of federation (K = 18 nodes)
# 2026-09-27 rerun with the v3 estimator. Replaces 2026_09_02_visn_analysis.R.
#
# RUN LOCATION: VA network / VINCI (needs 01_cleanData.dta). PHI stays inside;
# only coefficient-level summaries are written out.
#
# WHAT CHANGED FROM 2026_09_02 (and why)
#   1. Estimator file: 01_fedgee_v3.R (copy it from
#      inst/simulations/sim_v3/code/ into the VA working folder).
#   2. "fg" is now TRUE Fay-Graubard at the VISN level. The old "fg" was an
#      observation-level HC2; it is kept as "hc2obs" for the record.
#   3. Reference distribution: PRIMARY = KC + t(Bell-McCaffrey df). BM df is
#      computed in score space from the VISN breads only (design-based,
#      federated). sim_v3 (117 scenarios x 1000 reps): KC + t_BM covers
#      0.944-0.962 in every K x imbalance cell.
#   4. The old df (nu_j, realized-score count) is kept ONLY as a legacy column
#      to reconcile with the ENAR draft. sim_v3: it over-covers (0.96-0.98);
#      it drove the non-clinician "flip". Do not report it as a result.
#   5. Meta-analysis: per-VISN pre-screen for separation before geeglm
#      (geepack can hang on a separated site).
#
# ESTIMATORS (one estimand; beta identical across all federated variants)
#   pooled_visn (c)  - uncorrected VISN sandwich, z        (naive benchmark)
#   centralized (a)  - KC / MD / FG corrections, t(K-1), t(K-p), t(BM)
#   CR2 comparator   - clubSandwich CR2 + Satterthwaite (pooled-data benchmark)
#   meta (d)         - optional; per-VISN GEE, inverse-variance combined
#   decentralized (b)- optional; gossip L-sweep (unchanged from 09_02)
###############################################################################
setwd("P:/ORD_Anderson_202410002D/Experiments/FEDGEE paper")
rm(list = ls())

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(purrr)
  library(geepack)
})

`%||%`     <- function(a, b) if (is.null(a)) b else a
FEDGEE_SRC <- "01_fedgee_v3.R"
DEC_SRC    <- "01_decentralized_fedgee.R"
DATA_DTA   <- "01_cleanData.dta"
OUT_DIR    <- "2026_09_27_output"
dir.create(OUT_DIR, showWarnings = FALSE)

# The estimator files share internal names with different signatures.
# Source each into its own environment.
.fedgee_env <- new.env(parent = globalenv())
sys.source(FEDGEE_SRC, envir = .fedgee_env)
fedgee_v3 <- get("fedgee_v2", envir = .fedgee_env)  # v3 file keeps the v2 wrapper name
HAVE_SAWS         <- requireNamespace("saws", quietly = TRUE)
HAVE_CLUBSANDWICH <- requireNamespace("clubSandwich", quietly = TRUE)
if (!HAVE_SAWS) message("NOTE: saws not installed -> FG's own df (t_fgdf) will be NA. ",
                        "FG + t_BM is unaffected.")

# FG's own df via saws: all coefficients (cheap at K = 18, p = 12).
options(fedgee.df_coefs = NULL)

set.seed(20260927)

###############################################################################
# 0. DATA + SHARED MODEL (unchanged from 2026_09_02)
###############################################################################
load_maud <- function(dta_path = DATA_DTA) {
  df <- haven::read_dta(dta_path)
  df <- df %>%
    mutate(
      facility   = as.factor(sta6a_dis),   # discharge facility (sta6a)
      patienticn = as.factor(patienticn),
      visn       = as.factor(visnfy17)      # FY2017 VISN vintage; the FED UNIT
    ) %>%
    arrange(patienticn)

  # predictors: significant main effects from AoIM Tables 2-3 (dagger rows)
  df$consult3 <- relevel(factor(df$hospital_addiction_consult), ref = "0")
  df$psych    <- factor(as.integer(as.character(df$speciality_cat) == "1"),
                        levels = c(0, 1))
  df$oud      <- factor(df$oud_dx)
  df$icu_any  <- factor(df$icu)
  df$can_high <- factor(as.integer(as.character(df$CAN_quint) %in% c("4", "5")),
                        levels = c(0, 1))
  df$frail    <- factor(as.integer(as.character(df$frail_cat2) %in% c("3","4","5")),
                        levels = c(0, 1))
  df$ama      <- factor(as.integer(as.character(df$discharge_cat) == "2"),
                        levels = c(0, 1))
  df$age65    <- factor(as.integer(as.character(df$age_cat2) == "3"),
                        levels = c(0, 1))
  df$female   <- factor(as.integer(as.character(df$sex) == "1"),
                        levels = c(0, 1))
  # AUDIT-C no/low use (1,2) vs unhealthy+ (3,4,5)
  df$auditc_low <- factor(as.integer(as.character(df$auditc_baseline_cat_d)
                                     %in% c("1", "2")),
                          levels = c(0, 1))
  df
}

Y <- "new_during_aud_med"

# p = 1 intercept + 2 (consult3) + 9 binaries = 12; K = 18 -> K-p = 6, K-1 = 17
shared_formula <- as.formula(paste(
  Y,
  "~ consult3 + psych + oud + icu_any + can_high + frail + ama + age65 +",
  "female + auditc_low"
))

spot_check <- function(df) {
  cat("\n===== CODING SPOT-CHECK (confirm before trusting results) =====\n")
  vs <- df %>% count(visn, name = "n_rows") %>% arrange(n_rows)
  print(vs)
  cat(sprintf("\n#VISNs = %d   #facilities = %d   #rows = %d   VISN size CV = %.2f\n",
              nlevels(droplevels(df$visn)), nlevels(droplevels(df$facility)),
              nrow(df), sd(vs$n_rows) / mean(vs$n_rows)))
  for (v in c("hospital_addiction_consult","speciality_cat","discharge_cat",
              "CAN_quint","frail_cat2","age_cat2","sex","oud_dx","icu",
              "auditc_baseline_cat_d")) {
    if (v %in% names(df)) { cat("\n", v, ":\n", sep=""); print(table(df[[v]], useNA="ifany")) }
  }
  cat("\nOutcome ", Y, ":\n", sep=""); print(table(df[[Y]], useNA="ifany"))
  cat("Model terms after expansion (p): ",
      ncol(model.matrix(shared_formula, df)), "\n")
  cat("=================================================================\n\n")
}

# Abort before any correction if the VISN fit is not well-posed:
# (1) full column rank; (2) p <= K-1 so the VISN meat is not singular.
rank_guard <- function(df, formula, cluster_col = "visn") {
  df <- droplevels(df)
  X  <- model.matrix(formula, df)
  p  <- ncol(X)
  K  <- nlevels(droplevels(as.factor(df[[cluster_col]])))
  rX <- qr(X)$rank
  message(sprintf("[rank_guard] p=%d  K=%d  rank(X)=%d  K-p=%d  meat_rank=min(p,K-1)=%d",
                  p, K, rX, K - p, min(p, K - 1)))
  if (rX < p) stop(sprintf(
    "RANK GUARD FAILED: design rank %d < p %d (collinear dummy?). Inspect spot_check().", rX, p))
  if (p > K - 1) stop(sprintf(
    "RANK GUARD FAILED: p=%d > K-1=%d; VISN meat is singular. Reduce the model.", p, K - 1))
  invisible(list(p = p, K = K, rank_X = rX))
}

###############################################################################
# 1. CENTRALIZED FED-GEE (a) + POOLED (c), VISN-clustered
#    One fedgee_v3 fit returns every variance variant. For each variant we
#    build one row per coefficient per reference distribution.
###############################################################################
# reference distributions per variant (legacy t_nu kept for reconciliation)
REFS <- list(
  site      = c("z", "t_k1", "t_kp", "t_bm", "t_nu"),
  kc        = c("z", "t_k1", "t_kp", "t_bm", "t_nu"),
  md        = c("z", "t_k1", "t_kp", "t_bm", "t_nu"),
  fg        = c("z", "t_k1", "t_kp", "t_bm", "t_fgdf", "t_nu"),
  hc2obs    = c("t_bm", "t_nu"),
  patient_z = c("z")
)

run_centralized <- function(df, formula, family_obj = binomial("logit")) {
  df <- droplevels(df)
  data_list <- split(df, df$visn)            # 18 VISN nodes
  K <- length(data_list)
  fit <- fedgee_v3(
    data_list = data_list, main_formula = formula, family_obj = family_obj,
    corstr = "independence", id_col = "patienticn",
    n_iter = 50, tol = 1e-8, verbose = FALSE)
  if (is.null(fit)) stop("fedgee_v3 returned NULL (singular aggregation).")
  if (!isTRUE(fit$converged)) warning("fedgee_v3 did not converge in 50 iterations.")

  beta  <- fit$coefficients
  terms <- names(beta)
  p     <- length(beta)

  rows <- list()
  for (vn in names(REFS)) {
    v <- fit$variants[[vn]]
    if (is.null(v)) next
    for (ref in REFS[[vn]]) {
      dfv <- switch(ref,
        z      = rep(Inf, p),
        t_k1   = rep(v$df_k1, p),
        t_kp   = rep(v$df_kp, p),
        t_bm   = v$df_bm,
        t_fgdf = v$df_fg,
        t_nu   = v$df_nu)
      cm <- ifelse(is.finite(dfv) & dfv > 0, qt(0.975, pmax(dfv, 1e-8)), qnorm(0.975))
      if (ref != "z") cm[!is.finite(dfv)] <- NA   # missing df -> no interval, not z
      rows[[length(rows) + 1]] <- data.frame(
        variant = vn, inference = ref, K = K, p = p,
        term = terms, estimate = as.vector(beta), se = v$se, df = dfv,
        OR = exp(as.vector(beta)),
        CI_low  = exp(beta - cm * v$se), CI_high = exp(beta + cm * v$se),
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  out$p_value <- with(out, ifelse(is.finite(df), 2 * pt(-abs(estimate / se), df),
                                  2 * pnorm(-abs(estimate / se))))
  out$p_value[out$inference != "z" & !is.finite(out$df)] <- NA
  out$significant <- with(out, CI_low > 1 | CI_high < 1)
  out$role <- with(out, ifelse(variant == "kc" & inference == "t_bm", "PRIMARY",
                        ifelse(variant == "fg" & inference == "t_bm", "robustness",
                        ifelse(variant == "site" & inference == "z", "naive (pooled GEE)",
                        ifelse(inference == "t_nu", "legacy (ENAR draft; do not report)", "")))))

  # ---- CR2 + Satterthwaite comparator (clubSandwich) on the pooled data ----
  # Independence GEE = logistic glm, so fit glm: clubSandwich's geeglm method
  # requires patient ids nested in VISNs, and a patient can cross VISNs.
  cr2 <- tryCatch({
    if (!HAVE_CLUBSANDWICH) stop("clubSandwich not installed")
    g  <- glm(formula, data = df, family = family_obj)
    vc <- clubSandwich::vcovCR(g, cluster = df$visn, type = "CR2")
    ct <- clubSandwich::coef_test(g, vcov = vc, test = "Satterthwaite")
    data.frame(
      variant = "cr2", inference = "t_satt", K = K, p = p,
      term = rownames(ct), estimate = ct$beta, se = ct$SE, df = ct$df_Satt,
      OR = exp(ct$beta),
      CI_low  = exp(ct$beta - qt(0.975, ct$df_Satt) * ct$SE),
      CI_high = exp(ct$beta + qt(0.975, ct$df_Satt) * ct$SE),
      p_value = ct$p_Satt,
      significant = NA, role = "comparator (pooled data)",
      stringsAsFactors = FALSE)
  }, error = function(e) { message("  CR2 comparator skipped: ", e$message); NULL })
  if (!is.null(cr2)) {
    cr2$significant <- with(cr2, CI_low > 1 | CI_high < 1)
    out <- rbind(out, cr2)
  }

  list(estimates = out, diagnostics = fit$diagnostics,
       iterations = fit$iterations, converged = fit$converged)
}

# Headline table: one row per coefficient, primary vs naive vs legacy.
headline_table <- function(est) {
  pick <- function(vn, ref, lab) {
    est %>% filter(variant == vn, inference == ref) %>%
      transmute(term,
                !!paste0(lab, "_CI")  := sprintf("%.2f (%.2f, %.2f)", OR, CI_low, CI_high),
                !!paste0(lab, "_df")  := round(df, 1),
                !!paste0(lab, "_sig") := significant)
  }
  pick("site", "z", "naive") %>%
    left_join(pick("kc", "t_bm", "KC_BM"), by = "term") %>%
    left_join(pick("fg", "t_bm", "FG_BM"), by = "term") %>%
    left_join(pick("kc", "t_k1", "KC_K1"), by = "term") %>%
    left_join(pick("kc", "t_nu", "legacy_KC_nu"), by = "term") %>%
    mutate(flip_naive_to_primary = naive_sig != KC_BM_sig,
           primary_vs_legacy_differ = KC_BM_sig != legacy_KC_nu_sig)
}

###############################################################################
# 2. META-ANALYSIS (d) — optional. Per-VISN GEE, inverse-variance combined.
###############################################################################
run_meta <- function(df, formula, family_obj = binomial("logit"), SEP_THRESH = 10) {
  df <- droplevels(df)
  by_visn <- split(df, df$visn)
  fits <- imap(by_visn, function(d, v) {
    d <- droplevels(d)
    # pre-screen before geeglm: geepack can loop forever on a separated VISN
    yv <- d[[Y]]
    if (min(sum(yv == 1), sum(yv == 0)) < 2) return(list(visn = v, ok = FALSE))
    g0 <- tryCatch(suppressWarnings(glm(formula, data = d, family = family_obj)),
                   error = function(e) NULL)
    if (is.null(g0) || !g0$converged) return(list(visn = v, ok = FALSE))
    c0 <- coef(g0)
    if (any(!is.finite(c0)) || max(abs(c0)) > SEP_THRESH) return(list(visn = v, ok = FALSE))
    fv <- fitted(g0)
    if (any(fv < 1e-8 | fv > 1 - 1e-8)) return(list(visn = v, ok = FALSE))

    fit <- tryCatch(
      geeglm(formula, data = d, id = patienticn, family = family_obj,
             corstr = "independence"),
      error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit)) return(list(visn = v, ok = FALSE))
    b  <- coef(fit)
    se <- tryCatch(setNames(sqrt(diag(fit$geese$vbeta)), names(b)),
                   error = function(e) setNames(rep(NA_real_, length(b)), names(b)))
    if (any(!is.finite(se)) || any(se <= 0)) return(list(visn = v, ok = FALSE, beta = b))
    list(visn = v, ok = TRUE, beta = b, se = se)
  })
  ok  <- Filter(function(x) isTRUE(x$ok), fits)
  bad <- setdiff(names(by_visn), vapply(ok, `[[`, character(1), "visn"))
  if (length(ok) < 2) return(list(combined = NULL, n_ok = length(ok), dropped = bad))

  terms <- Reduce(union, lapply(ok, function(x) names(x$beta)))
  comb <- lapply(terms, function(tm) {
    num <- 0; den <- 0; k <- 0
    for (x in ok) {
      b_tm  <- x$beta[tm]
      se_tm <- if (tm %in% names(x$se)) x$se[[tm]] else NA_real_
      if (!is.na(b_tm) && is.finite(se_tm) && se_tm > 0) {
        w <- 1 / se_tm^2
        num <- num + w * b_tm; den <- den + w; k <- k + 1
      }
    }
    if (den == 0) data.frame(term = tm, estimate = NA_real_, se = NA_real_, n_visn = 0L)
    else data.frame(term = tm, estimate = unname(num/den), se = sqrt(1/den), n_visn = k)
  })
  comb <- do.call(rbind, comb)
  cm <- qnorm(0.975)
  comb <- transform(comb, estimator = "meta_analysis (d)", K = length(ok),
                    OR = exp(estimate),
                    CI_low = exp(estimate - cm*se), CI_high = exp(estimate + cm*se))
  list(combined = comb, n_ok = length(ok), dropped = bad)
}

###############################################################################
# 3. DECENTRALIZED FED-GEE (b) — optional; unchanged from 2026_09_02.
###############################################################################
build_visn_backbone_W <- function(K, backbone = c("complete","star","ring"), mh_weights) {
  backbone <- match.arg(backbone)
  A <- matrix(0, K, K)
  if (backbone == "complete") {
    A[,] <- 1
  } else if (backbone == "star") {
    A[1, ] <- 1; A[, 1] <- 1
  } else {
    for (i in seq_len(K)) {
      j <- if (i == K) 1L else i + 1L
      A[i, j] <- 1; A[j, i] <- 1
    }
  }
  diag(A) <- 1
  list(W = mh_weights(A), backbone = backbone)
}
rho_of_W <- function(W) {
  ev <- sort(abs(eigen(W, symmetric = TRUE, only.values = TRUE)$values), decreasing = TRUE)
  ev[2]
}

run_decentralized <- function(df, formula,
                              L_grid = c(10, 25, 50, 116, 200, 400),
                              backbones = c("complete","star","ring"),
                              family_obj = binomial("logit")) {
  .dec_env <- new.env(parent = globalenv())
  sys.source(DEC_SRC, envir = .dec_env)
  decentralized_fedgee <- get("decentralized_fedgee", envir = .dec_env)
  mh_weights           <- get("mh_weights",           envir = .dec_env)

  df <- droplevels(df)
  data_list <- split(df, df$visn)
  K <- length(data_list)
  out <- list()
  for (bb in backbones) {
    wm  <- build_visn_backbone_W(K, bb, mh_weights)
    rho <- rho_of_W(wm$W)
    for (L in L_grid) {
      fit <- tryCatch(
        decentralized_fedgee(
          data_list = data_list, main_formula = formula, family_obj = family_obj,
          corstr = "independence", id_col = "patienticn", W = wm$W,
          L_beta = L, L_S = L, L_B = L, sandwich_level = "site",
          md_correction = FALSE, n_iter = 50, tol = 1e-8, verbose = FALSE),
        error = function(e) { message(sprintf("  [%s L=%d] %s", bb, L, e$message)); NULL })
      if (is.null(fit)) next
      dfr <- fit$df_residual %||% (K - length(fit$coefficients))
      cm  <- if (dfr > 0) qt(0.975, dfr) else qnorm(0.975)
      out[[length(out) + 1]] <- data.frame(
        estimator = "decentralized_KC (b)", backbone = bb, rho = rho, rho_L = rho^L, L = L,
        K = K, df = dfr, term = names(fit$coefficients),
        estimate = as.vector(fit$coefficients), se = fit$se,
        OR = exp(as.vector(fit$coefficients)),
        CI_low = exp(fit$coefficients - cm * fit$se),
        CI_high = exp(fit$coefficients + cm * fit$se),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

###############################################################################
# 4. DRIVER
###############################################################################
main <- function(run_meta_analysis = TRUE, run_decentralized_sweep = FALSE) {
  df <- load_maud()
  spot_check(df)
  rank_guard(df, shared_formula, cluster_col = "visn")

  message(">> centralized Fed-GEE v3 (K = 18 VISNs): all variants + CR2 ...")
  cen <- run_centralized(df, shared_formula)
  dg  <- cen$diagnostics
  message(sprintf(paste0("   converged=%s iter=%d | VISN size CV=%.2f | max VISN leverage=%.3f",
                         " | KC clipped=%d MD clipped=%d FG capped=%d"),
                  cen$converged, cen$iterations, dg$site_size_cv, dg$max_leverage,
                  dg$kc_clipped, dg$md_clipped, dg$fg_capped))

  head_tab <- headline_table(cen$estimates)
  cat("\n===== HEADLINE: naive vs KC + t_BM (primary) vs FG + t_BM vs KC + t(K-1) =====\n")
  print(as.data.frame(head_tab), row.names = FALSE)
  cat("\nBM df per coefficient (KC):\n")
  print(cen$estimates %>% filter(variant == "kc", inference == "t_bm") %>%
          transmute(term, df_bm = round(df, 2)), row.names = FALSE)

  meta <- NULL
  if (run_meta_analysis) {
    message(">> meta-analysis (18 per-VISN GEEs) ...")
    meta <- run_meta(df, shared_formula)
    message(sprintf("   meta: %d/%d VISNs usable; dropped: %s",
                    meta$n_ok, nlevels(droplevels(df$visn)),
                    if (length(meta$dropped)) paste(meta$dropped, collapse = ", ") else "none"))
  }

  res <- list(centralized = cen$estimates, headline = head_tab,
              diagnostics = dg, meta = meta,
              formula = deparse(shared_formula),
              estimator_file = FEDGEE_SRC, saws = HAVE_SAWS,
              run_time = Sys.time(), session = utils::sessionInfo())
  saveRDS(res, file.path(OUT_DIR, "2026_09_27_visn_results_central.RDS"))
  write.csv(cen$estimates, file.path(OUT_DIR, "2026_09_27_visn_estimates_long.csv"), row.names = FALSE)
  write.csv(head_tab,      file.path(OUT_DIR, "2026_09_27_visn_headline.csv"),      row.names = FALSE)
  if (!is.null(meta$combined))
    write.csv(meta$combined, file.path(OUT_DIR, "2026_09_27_visn_meta.csv"), row.names = FALSE)
  message("   saved results to ", OUT_DIR)

  dec <- NULL
  if (run_decentralized_sweep) {
    message(">> decentralized L-sweep (slow) ...")
    dec <- run_decentralized(df, shared_formula)
    saveRDS(list(decentralized = dec, formula = deparse(shared_formula), run_time = Sys.time()),
            file.path(OUT_DIR, "2026_09_27_visn_results_decentralized.RDS"))
    write.csv(dec, file.path(OUT_DIR, "2026_09_27_visn_decentralized.csv"), row.names = FALSE)
  }
  invisible(c(res, list(decentralized = dec)))
}

results <- main()
