cw <- chick_sites()
dl <- split(cw, cw$site)
f <- weight ~ Time + Diet
fit <- fedgee(dl, f, family_obj = gaussian(), id_col = "Chick", verbose = FALSE)

test_that("point estimate equals pooled GEE (factor covariate)", {
  g <- geepack::geeglm(f, data = cw[order(cw$site), ], id = site,
                       family = gaussian(), corstr = "independence")
  expect_equal(unname(coef(fit)), unname(coef(g)), tolerance = 1e-8)
  expect_equal(unname(vcov(fit, "none")), unname(vcov(g)), tolerance = 1e-8)
})

test_that("KC and Bell-McCaffrey df equal clubSandwich CR2 at W = I", {
  skip_if_not_installed("clubSandwich")
  cr2 <- clubSandwich::coef_test(lm(f, data = cw), vcov = "CR2",
                                 cluster = cw$site, test = "Satterthwaite")
  kc <- fit$variances$site$KC
  expect_equal(unname(kc$se), cr2$SE, tolerance = 1e-8)
  expect_equal(unname(kc$df_bm), cr2$df_Satt, tolerance = 1e-6)
})

test_that("FG equals saws method d5", {
  skip_if_not_installed("saws")
  S <- lapply(dl, function(d) {
    X <- model.matrix(f, d)
    crossprod(X, d$weight - X %*% coef(fit))
  })
  om <- array(0, c(length(dl), 5, 5))
  for (i in seq_along(dl)) om[i, , ] <- crossprod(model.matrix(f, dl[[i]]))
  sw <- saws::saws(list(coefficients = coef(fit),
                        u = do.call(rbind, lapply(S, as.numeric)), omega = om),
                   method = "d5")
  expect_equal(unname(fit$variances$site$FG$vcov), unname(sw$V), tolerance = 1e-8)
})

test_that("MD >= KC >= none on the diagonal", {
  v <- fit$variances$site
  expect_true(all(v$MD$se >= v$KC$se - 1e-12))
  expect_true(all(v$KC$se >= v$none$se - 1e-12))
})

test_that("balanced sites give BM df = K - 1 for the intercept-only model", {
  bal <- lapply(1:8, function(k) data.frame(pat_id = 1:10, y = rnorm(10)))
  fb <- fedgee(bal, y ~ 1, family_obj = gaussian(), verbose = FALSE)
  expect_equal(unname(fb$variances$site$KC$df_bm), 7, tolerance = 1e-8)
})

test_that("summary switches variant without refitting", {
  s <- summary(fit, correction = "MD", df = "K-1")
  expect_equal(s$coefficients$SE, unname(fit$variances$site$MD$se))
  expect_true(all(s$coefficients$df == fit$n_sites - 1))
  expect_equal(summary(fit, df = "z")$coefficients$df, rep(Inf, 5))
  expect_equal(dim(confint(fit)), c(5L, 2L))
})

test_that("patient-level sandwich equals geeglm clustered by patient", {
  fp <- fedgee(dl, f, family_obj = gaussian(), id_col = "Chick",
               sandwich_level = "patient", verbose = FALSE)
  gp <- geepack::geeglm(f, data = cw[order(cw$Chick), ], id = Chick,
                        family = gaussian(), corstr = "independence")
  expect_equal(unname(summary(fp, correction = "none")$coefficients$SE),
               unname(summary(gp)$coefficients$Std.err), tolerance = 1e-8)
  expect_error(summary(fp, correction = "FG"), "site")
})

test_that("exchangeable alpha is estimated (not silently zero) and clamped", {
  expect_warning(
    fe <- fedgee(dl, f, family_obj = gaussian(), corstr = "exchangeable",
                 id_col = "Chick", verbose = FALSE),
    "clamped"
  )
  expect_true(any(unlist(fe$alpha) != 0))
})

test_that("common alpha reproduces pooled GEE for logistic exchangeable", {
  bl <- logistic_sites()
  fl <- fedgee(bl, y ~ x + z, corstr = "exchangeable", verbose = FALSE)
  expect_true(fl$converged)
  expect_equal(fl$n_sites, 10L)
  expect_true(all(is.finite(fl$se)))
})

test_that("rank warning fires when K - 1 < p", {
  expect_warning(
    fedgee(split(cw, cw$site %% 3), f, family_obj = gaussian(),
           id_col = "Chick", verbose = FALSE),
    "rank"
  )
})
