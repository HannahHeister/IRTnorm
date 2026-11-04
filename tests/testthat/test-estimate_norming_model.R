# ---- helpers ----------------------------------------------------------------
# Replaces the five helpers; records what each one received.
local_fake_pipeline <- function(.env = parent.frame()) {
  calls <- new.env()
  testthat::local_mocked_bindings(
    data_prep = function(...) {
      calls$data_prep <- list(...)
      "NEWDATA"
    },
    fit_IRTnorm = function(data, ...) {
      calls$fit_IRTnorm <- c(list(data = data), list(...))
      "FIT"
    },
    extract_info = function(fit, data) {
      calls$extract_info <- list(fit = fit, data = data)
      "MODEL_INFO"
    },
    info_new_norming = function(fit, data) {
      calls$info_new_norming <- list(fit = fit, data = data)
      "NORM_INFO"
    },
    fit_info_param = function(fit, data) {
      calls$fit_info_param <- list(fit = fit, data = data)
      "CONV_INFO"
    },
    .env = .env
  )
  calls
}

run_est <- function(...) {
  estimate_norming_model(raw_data = raw, age_variable = "age",
                         item_variables = c("item1", "item2"),
                         irt_model = "2PLnorm", age_model = "poly",
                         poly_mean = 1, poly_sd = 2, ...)
}

# ---- input check -------------------------------------------------------------
test_that("non-logical extract_stan_object errors before any work is done", {
  calls <- local_fake_pipeline()
  expect_error(run_est(extract_stan_object = "yes"), "has to be logical")
  expect_error(run_est(extract_stan_object = NULL), "has to be logical")
  expect_null(calls$data_prep)
})

# ---- return structure --------------------------------------------------------
test_that("result contains all components, including the Stan object by default", {
  local_fake_pipeline()
  out <- run_est()
  expect_named(out, c("model_info", "norm_info", "conv_info", "newdata",
                      "raw_data", "age_variable", "item_variables",
                      "model_stan_obj"))
})

test_that("Stan object is omitted when extract_stan_object = FALSE", {
  local_fake_pipeline()
  out <- run_est(extract_stan_object = FALSE)
  expect_false("model_stan_obj" %in% names(out))
  expect_named(out, c("model_info", "norm_info", "conv_info", "newdata",
                      "raw_data", "age_variable", "item_variables"))
})

test_that("each helper's result lands in the right slot", {
  local_fake_pipeline()
  out <- run_est()
  expect_equal(out$model_info, "MODEL_INFO")
  expect_equal(out$norm_info, "NORM_INFO")
  expect_equal(out$conv_info, "CONV_INFO")
  expect_equal(out$newdata, "NEWDATA")
  expect_equal(out$model_stan_obj, "FIT")
  expect_equal(out$raw_data, raw)
  expect_equal(out$age_variable, "age")
  expect_equal(out$item_variables, c("item1", "item2"))
})

# ---- data flow ---------------------------------------------------------------
test_that("fit_IRTnorm gets the prepared data; downstream helpers get fit and data", {
  calls <- local_fake_pipeline()
  run_est()
  expect_equal(calls$fit_IRTnorm$data, "NEWDATA")
  for (h in c("extract_info", "info_new_norming", "fit_info_param")) {
    expect_equal(calls[[h]]$fit, "FIT", info = h)
    expect_equal(calls[[h]]$data, "NEWDATA", info = h)
  }
})

# ---- argument forwarding -----------------------------------------------------
test_that("data_prep receives all model arguments, with defaults filled in", {
  calls <- local_fake_pipeline()
  run_est()
  a <- calls$data_prep
  expect_equal(a$raw_data, raw)
  expect_equal(a$irt_model, "2PLnorm")
  expect_equal(a$age_model, "poly")
  expect_equal(a$poly_mean, 1)
  expect_equal(a$poly_sd, 2)
  expect_equal(a$k_mean, 10)
  expect_equal(a$k_sd, 8)
  expect_equal(a$k_sd2, 8)
  expect_null(a$K_item)
  expect_null(a$prior_knowledge)
})

test_that("every optional model argument is forwarded to data_prep", {
  calls <- local_fake_pipeline()
  run_est(itemD = 1.7, K_item = 3, poly_sd_unique = 0.5, k_mean = 5, k_sd = 4,
          k_sd2 = 3, fixed_knots_mean = c(1, 2), fixed_knots_sd = c(3, 4),
          fixed_knots_sd2 = c(5, 6), prior_knowledge = "pk", parameter_fixed = "pf")
  a <- calls$data_prep
  expect_equal(a$itemD, 1.7);            expect_equal(a$K_item, 3)
  expect_equal(a$poly_sd_unique, 0.5);   expect_equal(a$k_mean, 5)
  expect_equal(a$k_sd, 4);               expect_equal(a$k_sd2, 3)
  expect_equal(a$fixed_knots_mean, c(1, 2))
  expect_equal(a$fixed_knots_sd, c(3, 4))
  expect_equal(a$fixed_knots_sd2, c(5, 6))
  expect_equal(a$prior_knowledge, "pk")
  expect_equal(a$parameter_fixed, "pf")
})

test_that("sampling settings are forwarded to fit_IRTnorm", {
  calls <- local_fake_pipeline()
  run_est(iter_warmup = 50, iter_sampling = 60, chains = 3,
          parallel_chains = 2, seed = 123)
  a <- calls$fit_IRTnorm
  expect_equal(a$iter_warmup, 50);   expect_equal(a$iter_sampling, 60)
  expect_equal(a$chains, 3);         expect_equal(a$parallel_chains, 2)
  expect_equal(a$seed, 123)
})


