
ages <- data.frame(age = seq(2, 18, length.out = 40))
k <- 6


test_that("Xraw %*% R reproduces Z and X for every knot setting", {
  for (fk in list(NULL, "equidistant", c(2, 5, 9, 13, 18))) {
    res  <- build_spline(ages, k, fk)
    Xraw <- mgcv::PredictMat(res$sm, ages)
    expect_equal(Xraw %*% res$R_Z, res$Z, tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(Xraw %*% res$R_X, res$X, tolerance = 1e-6, ignore_attr = TRUE)
  }
})


test_that("equidistant knots are built from the observed age range with k knots", {
  res <- build_spline(ages, k, "equidistant")
  expect_equal(unname(res$sm$xp), seq(min(ages$age), max(ages$age), length.out = k))
  expect_equal(res$K, k - 2)
})

test_that("returns the expected components and rotation dimensions", {
  res <- build_spline(ages, k, fixed_knots = NULL)
  expect_named(res, c("sm", "X", "Z", "S_fixed", "K", "R_Z", "R_X"))
  expect_equal(res$S_fixed, ncol(res$X))
  expect_equal(res$K, ncol(res$Z))
  expect_equal(dim(res$R_Z), c(k, res$K))
  expect_equal(dim(res$R_X), c(k, res$S_fixed))
})

test_that("rotated unpenalized part still spans constant + linear at new ages", {
  res <- build_spline(ages, k, NULL)
  new <- data.frame(age = c(2.5, 7.3, 11, 17.9))
  Xn_fixed <- mgcv::PredictMat(res$sm, new) %*% res$R_X
  for (target in list(rep(1, nrow(new)), new$age)) {
    expect_lt(max(abs(lm.fit(Xn_fixed, target)$residuals)), 1e-6)
  }
})

test_that("custom fixed knots override k", {
  knots <- c(2, 5, 9, 13, 18)
  res <- build_spline(ages, k = 10, fixed_knots = knots)   # k = 10 is ignored
  expect_equal(unname(res$sm$xp), knots)
  expect_equal(res$K, length(knots) - 2)
})