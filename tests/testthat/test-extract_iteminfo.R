
fake_fit <- function(iter = 3, chains = 2)
  list(metadata = function() list(iter_sampling = iter),
       num_chains = function() chains)

# get_value mock: returns recognisable values and records calls
mock_get_value <- function(fit, var, type) {
  if (type == "mean") seq_len(4) + switch(var, alpha = 0, beta = 10, gamma = 20, b = 30)
  else matrix(switch(var, alpha = 1, beta = 2, gamma = 3, b = 4), 6, 4)
}

run <- function(irt_mod, type, extra = list()) {
  data <- c(list(I = 4, irt_mod = irt_mod), extra)
  local_mocked_bindings(get_value = mock_get_value, .package = "IRTnorm")
  extract_iteminfo(fake_fit(), data, type)
}

test_that("type is validated", {
  expect_error(run(1, "bogus"))
})

test_that("mean: 2PL-type models estimate alpha and beta, gamma defaults to 0", {
  res <- run(2, "mean")
  expect_s3_class(res, "data.frame")
  expect_equal(res$alpha, 1:4)
  expect_equal(res$beta, 11:14)
  expect_equal(res$gamma, rep(0, 4))
})

test_that("mean: models without discrimination fix alpha = 1", {
  for (m in c(1, 5, 7, 9, 11)) expect_equal(run(m, "mean")$alpha, rep(1, 4), info = m)
})

test_that("mean: only irt_mod 3 estimates gamma", {
  expect_equal(run(3, "mean")$gamma, 21:24)
  for (m in c(1, 2, 5)) expect_equal(run(m, "mean")$gamma, rep(0, 4), info = m)
})

test_that("draws: all components are ndraws x nitem matrices", {
  for (m in c(1, 2, 3, 5, 9)) {
    res <- run(m, "draws")
    expect_type(res, "list")
    for (v in c("alpha", "beta", "gamma"))
      expect_equal(dim(res[[v]]), c(6, 4), info = paste(m, v))
  }
})

test_that("draws: defaults are constant matrices", {
  res <- run(1, "draws")
  expect_true(all(res$alpha == 1)); expect_true(all(res$gamma == 0))
})

test_that("fixed item parameters in data override estimation", {
  res <- run(2, "mean", list(alpha = c(1, 2, 3, 4), beta = c(5, 6, 7, 8)))
  expect_equal(res$alpha, c(1, 2, 3, 4))
  expect_equal(res$beta, c(5, 6, 7, 8))
})

test_that("fixed parameters in draws mode are replicated per draw, one row each", {
  res <- run(2, "draws", list(alpha = c(1, 2, 3, 4)))
  expect_equal(dim(res$alpha), c(6, 4))
  expect_true(all(apply(res$alpha, 1, identical, c(1, 2, 3, 4))))
})

test_that("GRM (irt_mod 4) returns alpha and threshold only", {
  res <- run(4, "draws")
  expect_named(res, c("alpha", "threshold"))
  expect_equal(res$threshold, matrix(4, 6, 4))
  expect_null(res$beta)
})

test_that("get_value is only called for parameters that are not fixed", {
  calls <- character()
  local_mocked_bindings(
    get_value = function(fit, var, type) { calls <<- c(calls, var); mock_get_value(fit, var, type) },
    .package = "IRTnorm")
  extract_iteminfo(fake_fit(), list(I = 4, irt_mod = 3, alpha = 1:4), "mean")
  expect_setequal(calls, c("beta", "gamma"))
})

test_that("ndraws = iter_sampling * num_chains", {
  res <- extract_iteminfo(fake_fit(iter = 5, chains = 3),
                          list(I = 4, irt_mod = 1, alpha = 1:4, beta = 1:4, gamma = rep(0, 4)),
                          "draws")
  expect_equal(nrow(res$alpha), 15)
})