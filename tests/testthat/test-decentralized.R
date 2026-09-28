cw <- chick_sites()
dl <- split(cw, cw$site)
f <- weight ~ Time + Diet
cen <- fedgee(dl, f, family_obj = gaussian(), id_col = "Chick", verbose = FALSE)

test_that("weight matrices are doubly stochastic and symmetric", {
  for (s in c("hub", "ring", "complete")) {
    W <- build_weight_matrix(10, s)$W
    expect_equal(rowSums(W), rep(1, 10))
    expect_equal(W, t(W))
  }
  v <- build_weight_matrix(8, "visn", region_id = rep(1:4, each = 2), hub_sites = c(1, 3, 5, 7))
  expect_equal(colSums(v$W), rep(1, 8))
  expect_lt(build_weight_matrix(10, "complete")$rho, 1e-8)
  expect_error(mh_weights(matrix(0, 3, 3)), "self-loops")
})

test_that("complete graph reproduces the centralized fit (none, KC, MD)", {
  for (cr in c("none", "KC", "MD")) {
    d <- decentralized_fedgee(dl, f, family_obj = gaussian(), id_col = "Chick",
                              structure = "complete", sandwich_level = "site",
                              correction = cr, ridge = 0, n_iter = 200, verbose = FALSE)
    expect_true(d$converged)
    expect_equal(unname(coef(d)), unname(coef(cen)), tolerance = 1e-7)
    expect_equal(unname(d$se), unname(cen$variances$site[[cr]]$se), tolerance = 1e-6)
  }
})

test_that("MD is now distinct from KC", {
  args <- list(dl, f, family_obj = gaussian(), id_col = "Chick", structure = "complete",
               sandwich_level = "site", ridge = 0, n_iter = 200, verbose = FALSE)
  kc <- do.call(decentralized_fedgee, c(args, correction = "KC"))
  md <- do.call(decentralized_fedgee, c(args, correction = "MD"))
  expect_true(all(md$se > kc$se))
})

test_that("deprecated md_correction maps to KC with a warning", {
  expect_warning(
    d <- decentralized_fedgee(dl, f, family_obj = gaussian(), id_col = "Chick",
                              structure = "complete", sandwich_level = "site",
                              md_correction = TRUE, n_iter = 200, verbose = FALSE),
    "deprecated"
  )
  expect_equal(d$correction, "KC")
})

test_that("W or structure is required", {
  expect_error(decentralized_fedgee(dl, f, family_obj = gaussian(), id_col = "Chick",
                                    verbose = FALSE), "structure")
})

test_that("ring with enough rounds approaches the centralized estimate", {
  d <- decentralized_fedgee(dl, f, family_obj = gaussian(), id_col = "Chick",
                            structure = "ring", L_beta = 200, L_S = 200, L_B = 200,
                            n_iter = 200, tol = 1e-6, verbose = FALSE)
  expect_equal(unname(coef(d)), unname(coef(cen)), tolerance = 1e-4)
})
