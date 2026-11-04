# ---- helpers ---------------------------------------------------------------

# Minimal stand-in for a loo object, with just the parts fit_one_model() reads.
# estimates[[1]] is the first element of the matrix, i.e. elpd_loo / Estimate.
# elpd_loo (-100) differs from p_loo (2) and looic (200), so reading the wrong
# row would be caught.
fake_loo <- function(n = 4, elpd = -100) {
  list(
    pointwise = cbind(elpd_loo = rep(elpd / n, n),
                      p_loo    = rep(0.5, n),
                      looic    = rep(-2 * elpd / n, n)),
    estimates = matrix(c(elpd, 2, -2 * elpd, 10, 0.5, 20), nrow = 3,
                       dimnames = list(c("elpd_loo", "p_loo", "looic"),
                                       c("Estimate", "SE"))),
    diagnostics = list(pareto_k = c(0.1, 0.2, 0.3, 0.4),
                       n_eff    = c(400, 410, 420, 430),
                       r_eff    = c(0.9, 0.8, 0.7, 0.6))
  )
}

# Replaces data_prep(), fit_IRTnorm() and loo_extr() so no Stan is involved.
# Returns an environment that records the arguments each stub received.
local_fake_pipeline <- function(.env = parent.frame()) {
  calls   <- new.env()
  newdata <- list(age_variable = c(1, 2, 3), J = 3L)
  fake_fit <- list(
    summary = function(variables, ...) {
      data.frame(mean = if (variables == "theta_mean") c(0.1, 0.2, 0.3)
                 else c(1.0, 1.1, 1.2))
    }
  )
  testthat::local_mocked_bindings(
    data_prep   = function(...) { calls$data_prep <- list(...); newdata },
    fit_IRTnorm = function(...) { calls$fit_IRTnorm <- list(...); fake_fit },
    loo_extr    = function(data, fit) {
      calls$loo_extr <- list(data = data, fit = fit)
      fake_loo()
    },
    .env = .env
  )
  calls
}

poly_spec <- list(irt_model = "2PLnorm", age_model = "poly", poly_mean = 1, poly_sd = 2)

run_fit <- function(spec = poly_spec, ...) {
  fit_one_model(spec, raw_data = raw, age_variable = "age",
                item_variables = c("item1", "item2"), ...)
}

# ---- return structure ------------------------------------------------------

test_that("return value has the expected elements and no trajectory by default", {
  local_fake_pipeline()
  out <- run_fit()
  expect_named(out, c("elpd_pointwise", "elpd_summary", "paretok", "neff", "reff"))
})

test_that("pointwise ELPD and diagnostics are passed through from loo", {
  local_fake_pipeline()
  out <- run_fit()
  expect_equal(out$elpd_pointwise, rep(-25, 4))
  expect_equal(out$paretok, c(0.1, 0.2, 0.3, 0.4))
  expect_equal(out$neff, c(400, 410, 420, 430))
  expect_equal(out$reff, c(0.9, 0.8, 0.7, 0.6))
})

test_that("ELPD summary uses the elpd_loo estimate, not another row", {
  local_fake_pipeline()
  out <- run_fit()
  expect_equal(out$elpd_summary$ELPD, -100)   # not p_loo (2) or looic (200)
  expect_equal(out$elpd_summary$irt_model, "2PLnorm")
  expect_equal(out$elpd_summary$age_model, "poly")
  expect_equal(nrow(out$elpd_summary), 1)
})

# ---- model label -----------------------------------------------------------

test_that("modelcombi is poly_mean_poly_sd for polynomial specs", {
  local_fake_pipeline()
  expect_equal(run_fit()$elpd_summary$modelcombi, "1_2")
})

test_that("modelcombi is irt_model_age_model when no poly settings are given", {
  local_fake_pipeline()
  spec <- list(irt_model = "2PLnorm", age_model = "spline")  # "spline" is a placeholder
  expect_equal(run_fit(spec)$elpd_summary$modelcombi, "2PLnorm_spline")
})

test_that("poly_mean and poly_sd are NA when absent", {
  local_fake_pipeline()
  spec <- list(irt_model = "2PLnorm", age_model = "spline")
  out <- run_fit(spec)
  expect_true(is.na(out$elpd_summary$poly_mean))
  expect_true(is.na(out$elpd_summary$poly_sd))
})

test_that("poly_mean and poly_sd are kept when given", {
  local_fake_pipeline()
  out <- run_fit()
  expect_equal(out$elpd_summary$poly_mean, 1)
  expect_equal(out$elpd_summary$poly_sd, 2)
})

# ---- trajectory ------------------------------------------------------------

test_that("traj_info = TRUE adds a correctly shaped trajectory", {
  local_fake_pipeline()
  out <- run_fit(traj_info = TRUE)
  expect_true("trajectory" %in% names(out))
  traj <- out$trajectory
  expect_equal(nrow(traj), 6)
  expect_equal(traj$est, c(0.1, 0.2, 0.3, 1.0, 1.1, 1.2))
  expect_equal(traj$age, rep(c(1, 2, 3), 2))
  expect_equal(traj$moment, rep(c("mean", "sd"), each = 3))
  expect_equal(unique(traj$modelcombi), "1_2")
})

# ---- argument handling -----------------------------------------------------

test_that("data_prep receives the data arguments and the model specification", {
  calls <- local_fake_pipeline()
  run_fit()
  args <- calls$data_prep
  expect_equal(args$raw_data, raw)
  expect_equal(args$age_variable, "age")
  expect_equal(args$item_variables, c("item1", "item2"))
  expect_equal(args$irt_model, "2PLnorm")
  expect_equal(args$age_model, "poly")
  expect_equal(args$poly_mean, 1)
})

test_that("defaults are filled in for unspecified options", {
  calls <- local_fake_pipeline()
  run_fit()
  expect_equal(calls$data_prep$k_mean, 10)
  expect_equal(calls$data_prep$k_sd, 8)
  expect_equal(calls$data_prep$k_sd2, 8)
  expect_true("K_item" %in% names(calls$data_prep))   # NULL default still passed
  expect_null(calls$data_prep$K_item)
})

test_that("user-specified options override the defaults", {
  calls <- local_fake_pipeline()
  run_fit(c(poly_spec, list(k_mean = 20)))
  expect_equal(calls$data_prep$k_mean, 20)
  expect_equal(calls$data_prep$k_sd, 8)               # other defaults untouched
})

test_that("fit_IRTnorm receives the prepared data and sampling settings", {
  calls <- local_fake_pipeline()
  run_fit(iter_warmup = 50, iter_sampling = 60, chains = 3,
          parallel_chains = 2, seed = 123)
  args <- calls$fit_IRTnorm
  expect_equal(args$data$J, 3L)
  expect_equal(args$iter_warmup, 50)
  expect_equal(args$iter_sampling, 60)
  expect_equal(args$chains, 3)
  expect_equal(args$parallel_chains, 2)
  expect_equal(args$seed, 123)
})

test_that("default sampling settings are 500/500 with 4 chains", {
  calls <- local_fake_pipeline()
  run_fit()
  expect_equal(calls$fit_IRTnorm$iter_warmup, 500)
  expect_equal(calls$fit_IRTnorm$iter_sampling, 500)
  expect_equal(calls$fit_IRTnorm$chains, 4)
})
