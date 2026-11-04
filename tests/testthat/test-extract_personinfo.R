library(testthat)

J <- 3; D <- 2; I <- 4; ND <- 5
age <- c(4, 5, 6)

# substraits: ND x (I*J) matrix, column k has constant value k
make_fit <- function(J = 3, D = 2, I = 4, mu = 1, sd = 2, mu_cols = D * J) {
  substraits <- matrix(rep(seq_len(I * J), each = ND), ND, I * J)
  list(
    summary = function(variables, ...) {
      switch(variables,
             z_theta     = data.frame(variable = paste0("z_theta[", 1:J, "]"),
                                      mean = c(0.1, 0.2, 0.3), q5 = c(-1, -1, -1), q95 = c(1, 1, 1)),
             theta       = data.frame(variable = "theta", mean = c(0.5, 0.6, 0.7)),
             theta_mean  = data.frame(variable = "theta_mean", mean = c(10, 11, 12)),
             theta_sigma = data.frame(variable = "theta_sigma", mean = c(2, 3, 4)))
    },
    draws = function(variables, format = "draws_matrix") {
      switch(variables,
             substraits  = substraits,
             theta_mean  = matrix(mu, ND, mu_cols),
             theta_sigma = matrix(sd, ND, mu_cols))
    }
  )
}

make_data <- function(irt_mod, ii = rep(1:4, each = 3)) {   # item-major ordering
  list(irt_mod = irt_mod, I = I, J = J, D = D, itemD = c(1, 1, 2, 2),
       ii = ii, age_variab = age)
}

# ---- irt_mod <= 4 ------------------------------------------------------------
test_that("irt_mod <= 4 returns one row per person with expected columns", {
  for (m in 1:4) {
    res <- extract_personinfo(make_fit(), make_data(m))
    expect_equal(nrow(res), J)
    expect_named(res, c("age", "normscore_est", "normscore_est_upper",
                        "normscore_est_lower", "theta_est", "theta_mu",
                        "theta_sigma", "irt_mod"))
    expect_equal(res$age, age)
    expect_equal(res$normscore_est, c(0.1, 0.2, 0.3))
    expect_equal(res$normscore_est_upper, rep(1, J))
    expect_equal(res$normscore_est_lower, rep(-1, J))
    expect_equal(res$theta_est, c(0.5, 0.6, 0.7))
    expect_equal(res$theta_mu, c(10, 11, 12))
    expect_equal(res$theta_sigma, c(2, 3, 4))
    expect_true(all(res$irt_mod == m))
  }
})

# ---- irt_mod > 4 -------------------------------------------------------------
test_that("higher-order/testlet output has J*(2D+1) rows in the right layout", {
  res <- extract_personinfo(make_fit(), make_data(7))
  expect_equal(nrow(res), J * (2 * D + 1))
  expect_named(res, c("est_ability", "ability_est_upper", "ability_est_lower",
                      "trait", "theta_mu", "theta_sigma", "age", "irt_mod"))
  expect_equal(res$trait,
               c(rep(c("delta1", "delta2"), each = J),
                 rep(c("normeddelta1", "normeddelta2"), each = J),
                 rep("theta", J)))
  expect_equal(res$age, rep(age, 2 * D + 1))
  expect_true(all(res$irt_mod == 7))
})

test_that("delta rows pick the first item of each dimension (item-major data)", {
  # first item of dim 1 = item 1 -> cols 1:3; of dim 2 = item 3 -> cols 7:9
  res <- extract_personinfo(make_fit(), make_data(7))
  expect_equal(res$est_ability[1:6], c(1:3, 7:9))
  # normed = (delta - mu) / sd with mu = 1, sd = 2
  expect_equal(res$est_ability[7:12], (c(1:3, 7:9) - 1) / 2)
  expect_equal(res$est_ability[13:15], c(0.1, 0.2, 0.3))   # theta from z_theta
})

test_that("person-major ordering of ii still yields person order within dimension", {
  # Rows ordered person then item. This test documents the ordering assumption:
  # it FAILS with the current code if int_delta is not built dimension by dimension.
  ii_pm <- rep(1:4, times = 3)
  fit <- make_fit()
  res <- extract_personinfo(fit, make_data(7, ii = ii_pm))
  # expected: persons 1..3 within each dimension, so the delta values must
  # increase with person index inside each block of J rows
  d1 <- res$est_ability[res$trait == "delta1"]
  expect_equal(d1, sort(d1))
})

test_that("theta_mu / theta_sigma only filled for the delta rows", {
  res <- extract_personinfo(make_fit(), make_data(7))
  expect_false(anyNA(res$theta_mu[1:(D * J)]))
  expect_true(all(is.na(res$theta_mu[-(1:(D * J))])))
  expect_true(all(is.na(res$theta_sigma[-(1:(D * J))])))
  expect_equal(res$theta_mu[1:(D * J)], rep(1, D * J))
  expect_equal(res$theta_sigma[1:(D * J)], rep(2, D * J))
})

test_that("models 5/6 replicate per-person mu/sd across dimensions", {
  for (m in 5:6) {
    res_small <- extract_personinfo(make_fit(mu_cols = J), make_data(m))
    res_full  <- extract_personinfo(make_fit(mu_cols = D * J), make_data(7))
    expect_equal(res_small$est_ability, res_full$est_ability, info = m)
    expect_equal(nrow(res_small), J * (2 * D + 1))
  }
})

test_that("intervals bracket the estimate when draws vary", {
  fit <- make_fit()
  base <- fit$draws
  fit$draws <- function(variables, format = "draws_matrix") {
    x <- base(variables, format)
    if (variables == "substraits") x <- x + seq(-1, 1, length.out = ND)  # per-draw shift
    x
  }
  res <- extract_personinfo(fit, make_data(7))
  d <- res[res$trait %in% c("delta1", "delta2", "normeddelta1", "normeddelta2"), ]
  expect_true(all(d$ability_est_lower < d$est_ability))
  expect_true(all(d$est_ability < d$ability_est_upper))
})

test_that("output has no NA in estimates and no leaked row names", {
  res <- extract_personinfo(make_fit(), make_data(7))
  expect_false(anyNA(res[, c("est_ability", "ability_est_upper", "ability_est_lower", "trait", "age")]))
  expect_equal(rownames(res), as.character(seq_len(nrow(res))))
})