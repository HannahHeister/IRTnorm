
skip_if_not_installed("posterior")

# ---- Helpers ---------------------------------------------------------------

all_vars <- c("z_theta[1]", "z_theta[2]", "p1_theta", "p2_theta", "p3_theta",
              "l1", "l2", "l3", "u1", "u2", "u3",
              "unique[1]", "unique[2]", "nu[1]", "nu[2]",
              "lambda[1]", "lambda[2]", "unique_sigma[1]", "unique_sigma[2]",
              "alpha[1]", "alpha[2]", "beta[1]", "beta[2]",
              "gamma[1]", "gamma[2]", "b[1]", "b[2]")

# Mock cmdstanr fit: only $draws() is used by the function.
# shift_var: variable whose chains are pulled apart (-> large rhat)
make_fit <- function(vars = all_vars, n_chain = 4, n_iter = 200, seed = 1,
                     shift_var = NULL) {
  set.seed(seed)
  arr <- array(rnorm(n_iter * n_chain * length(vars)),
               dim = c(n_iter, n_chain, length(vars)),
               dimnames = list(iteration = as.character(seq_len(n_iter)),
                               chain = as.character(seq_len(n_chain)),
                               variable = vars))
  if (!is.null(shift_var)) {
    arr[, , shift_var] <- arr[, , shift_var] +
      matrix(10 * (seq_len(n_chain) - 1), n_iter, n_chain, byrow = TRUE)
  }
  d <- posterior::as_draws_array(arr)
  list(draws = function() d)
}

params_of <- function(res) unique(res$variable)

# ---- Cmdstanr requirement --------------------------------------------------

test_that("errors with install instructions when cmdstanr is missing", {
  # Needs testthat >= 3.2 to mock a base function; skip otherwise
  skip_if(packageVersion("testthat") < "3.2.0")
  local_mocked_bindings(requireNamespace = function(...) FALSE, .package = "base")
  expect_error(fit_info_param(make_fit(), list(age_mod = 1, irt_mod = 1)),
               "cmdstanr")
})

# Everything below needs the real cmdstanr check to pass
skip_if_not_installed("cmdstanr")

# ---- Parameter selection ---------------------------------------------------

test_that("age_mod 1, irt_mod 1: theta, polynomial params and beta", {
  res <- fit_info_param(make_fit(), list(age_mod = 1, irt_mod = 1))
  expect_equal(params_of(res), c("z_theta", "p1_theta", "p2_theta", "beta"))
})

test_that("age_mod 2 uses the l/u parameters instead of p*_theta", {
  res <- fit_info_param(make_fit(), list(age_mod = 2, irt_mod = 1))
  expect_equal(params_of(res), c("z_theta", "l1", "l2", "u1", "u2", "beta"))
})

test_that("irt_mod > 8 adds the third age parameter(s)", {
  r1 <- fit_info_param(make_fit(), list(age_mod = 1, irt_mod = 9))
  r2 <- fit_info_param(make_fit(), list(age_mod = 2, irt_mod = 9))
  expect_true("p3_theta" %in% params_of(r1))
  expect_true(all(c("l3", "u3") %in% params_of(r2)))
  # not added for irt_mod <= 8
  expect_false("p3_theta" %in% params_of(
    fit_info_param(make_fit(), list(age_mod = 1, irt_mod = 8))))
})

test_that("each model family selects its own item and relation parameters", {
  get <- function(m) params_of(fit_info_param(make_fit(), list(age_mod = 1, irt_mod = m)))
  expect_true("gamma" %in% get(3))                       # 3PL
  expect_true(all(c("alpha", "b") %in% get(4)) && !"beta" %in% get(4))
  expect_true(all(c("unique", "nu") %in% get(5)))        # multidimensional
  expect_true("lambda" %in% get(9))
  expect_true("unique_sigma" %in% get(11))
  expect_false("alpha" %in% get(1))                      # Rasch-type: no alpha
})

test_that("parameters supplied as fixed data are not diagnosed", {
  data <- list(age_mod = 1, irt_mod = 3, alpha = 1, gamma = 0.2)
  res <- fit_info_param(make_fit(), data)
  expect_false(any(c("alpha", "gamma") %in% params_of(res)))
  expect_true("beta" %in% params_of(res))
  # data entries that are not model parameters are ignored
  res2 <- fit_info_param(make_fit(), list(age_mod = 1, irt_mod = 1, J = 100))
  expect_equal(params_of(res2), c("z_theta", "p1_theta", "p2_theta", "beta"))
})

# ---- Output structure ------------------------------------------------------

test_that("output has one row per diagnostic and parameter, in block order", {
  res <- fit_info_param(make_fit(), list(age_mod = 1, irt_mod = 1))
  p <- c("z_theta", "p1_theta", "p2_theta", "beta")
  expect_s3_class(res, "data.frame")
  expect_named(res, c("diagnostic", "variable",
                      "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."))
  expect_equal(nrow(res), 3 * length(p))
  expect_equal(res$diagnostic, rep(c("rhat", "ess_bulk", "ess_tail"), each = 4))
  expect_equal(res$variable, rep(p, times = 3))
  expect_equal(rownames(res), as.character(seq_len(nrow(res))))
})

# ---- Diagnostic values -----------------------------------------------------

test_that("well-mixed chains give rhat near 1 and large effective sample sizes", {
  res <- fit_info_param(make_fit(n_chain = 4, n_iter = 500),
                        list(age_mod = 1, irt_mod = 1))
  expect_true(all(res$Max.[res$diagnostic == "rhat"] < 1.05))
  expect_true(all(res$Min.[res$diagnostic == "ess_bulk"] > 400))   # of 2000 draws
})

test_that("diagnostics are summarised across the elements of a vector parameter", {
  res <- fit_info_param(make_fit(shift_var = "beta[1]"),
                        list(age_mod = 1, irt_mod = 1))
  beta_rhat <- res[res$diagnostic == "rhat" & res$variable == "beta", ]
  expect_gt(beta_rhat$Max., 1.5)        # beta[1]: chains far apart
  expect_lt(beta_rhat$Min., 1.05)       # beta[2]: fine
  expect_true(beta_rhat$Min. <= beta_rhat$Median && beta_rhat$Median <= beta_rhat$Max.)
})

test_that("a scalar parameter has identical summary values", {
  res <- fit_info_param(make_fit(), list(age_mod = 1, irt_mod = 1))
  row <- res[res$diagnostic == "rhat" & res$variable == "p1_theta", 3:8]
  expect_equal(as.numeric(row), rep(row$Min., 6))
})

test_that("a problem in one parameter does not leak into the others", {
  res <- fit_info_param(make_fit(shift_var = "p1_theta"),
                        list(age_mod = 1, irt_mod = 1))
  r <- res[res$diagnostic == "rhat", ]
  expect_gt(r$Max.[r$variable == "p1_theta"], 1.5)
  expect_lt(r$Max.[r$variable == "p2_theta"], 1.05)
})
