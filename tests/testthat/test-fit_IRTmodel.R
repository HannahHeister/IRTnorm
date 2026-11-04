
# Mocks: record what get_compiled_model() and mod$sample() receive
setup_mocks <- function(env = parent.frame()) {
  rec <- new.env()
  rec$compile <- NULL; rec$sample <- NULL
  testthat::local_mocked_bindings(
    get_compiled_model = function(age_model, irt_model, prior_knowledge, parameter_fixed) {
      rec$compile <- list(age_model = age_model, irt_model = irt_model,
                          prior_knowledge = prior_knowledge,
                          parameter_fixed = parameter_fixed)
      list(sample = function(...) { rec$sample <- list(...); "FIT" })
    },
    .package = "IRTnorm", .env = env
  )
  testthat::local_mocked_bindings(
    check_cmdstan_toolchain = function(...) invisible(TRUE),
    .package = "cmdstanr", .env = env
  )
  rec
}

base_data <- function(irt_mod = 2, age_mod = 1, ...)
  c(list(age_mod = age_mod, irt_mod = irt_mod, sm1 = 1, sm2 = 2, sm3 = 3), list(...))

fit_it <- function(data, ...) fit_IRTnorm(data, chains = 1, parallel_chains = 1, ...)

# ---- model selection ---------------------------------------------------------
test_that("age_mod and irt_mod are mapped to the right model names", {
  rec <- setup_mocks()
  irt_names <- c("1PLnorm","2PLnorm","3PLnorm","GRMnorm","Testlet1PLnorm","Testlet2PLnorm",
                 "HO1PLnorm-C","HO2PLnorm-C","HO1PLnorm-P","HO2PLnorm-P",
                 "HO1PLnorm-F","HO2PLnorm-F")
  for (m in 1:12) {
    fit_it(base_data(irt_mod = m, age_mod = 1))
    expect_equal(rec$compile$irt_model, irt_names[m], info = m)
    expect_equal(rec$compile$age_model, "poly")
  }
  fit_it(base_data(age_mod = 2))
  expect_equal(rec$compile$age_model, "splines")
})

test_that("returns the fit object from mod$sample()", {
  setup_mocks()
  expect_identical(fit_it(base_data()), "FIT")
})

# ---- data handling -----------------------------------------------------------
test_that("sm1-sm3 are dropped only for age_mod 2", {
  rec <- setup_mocks()
  fit_it(base_data(age_mod = 2))
  expect_false(any(c("sm1", "sm2", "sm3") %in% names(rec$sample$data)))
  fit_it(base_data(age_mod = 1))
  expect_true(all(c("sm1", "sm2", "sm3") %in% names(rec$sample$data)))
})

test_that("fixed item parameters are detected from names(data)", {
  rec <- setup_mocks()
  fit_it(base_data())
  expect_null(rec$compile$parameter_fixed)
  
  fit_it(base_data(alpha = 1:3, gamma = rep(0, 3)))
  expect_setequal(rec$compile$parameter_fixed, c("alpha", "gamma"))
  
  fit_it(base_data(irt_mod = 4, b = matrix(0, 3, 2)))
  expect_equal(rec$compile$parameter_fixed, "b")
})

test_that("known priors are translated to parameter names", {
  rec <- setup_mocks()
  fit_it(base_data())
  expect_null(rec$compile$prior_knowledge)
  
  fit_it(base_data(alpha_prior_mu = 1, alpha_prior_sd = 1, beta_prior_mu = 0))
  expect_setequal(rec$compile$prior_knowledge, c("alpha", "beta"))
  
  fit_it(base_data(gamma_prior_mu = 0, b_prior_mu = 0, nu_prior_sd = 1))
  expect_setequal(rec$compile$prior_knowledge, c("gamma", "b", "nu"))
})

# ---- sampler arguments -------------------------------------------------------
test_that("sampling arguments are forwarded", {
  rec <- setup_mocks()
  fit_IRTnorm(base_data(), seed = 11, iter_warmup = 200, iter_sampling = 300,
              chains = 2, parallel_chains = 2)
  s <- rec$sample
  expect_equal(s$seed, 11)
  expect_equal(s$iter_warmup, 200)
  expect_equal(s$iter_sampling, 300)
  expect_equal(s$chains, 2)
  expect_equal(s$parallel_chains, 2)
  expect_equal(s$refresh, 50)           # (200 + 300) / 10
})

test_that("seed: random when NULL (and controlled by set.seed), integer-valued", {
  rec <- setup_mocks()
  set.seed(1); fit_it(base_data()); s1 <- rec$sample$seed
  set.seed(1); fit_it(base_data()); s2 <- rec$sample$seed
  set.seed(2); fit_it(base_data()); s3 <- rec$sample$seed
  expect_equal(s1, s2)
  expect_false(s1 == s3)
  expect_true(is.numeric(s1) && s1 >= 1 && s1 <= .Machine$integer.max)
})

test_that("parallel_chains defaults to the number of cores", {
  rec <- setup_mocks()
  fit_IRTnorm(base_data(), chains = 1)
  expect_equal(rec$sample$parallel_chains, parallel::detectCores())
})

test_that("init = 1 is set for the GRM only", {
  rec <- setup_mocks()
  fit_it(base_data(irt_mod = 4))
  expect_equal(rec$sample$init, 1)
  for (m in c(1, 2, 3, 5, 9, 12)) {
    fit_it(base_data(irt_mod = m))
    expect_null(rec$sample$init, info = m)
  }
})

test_that("extra arguments in ... are passed to sample()", {
  rec <- setup_mocks()
  fit_it(base_data(irt_mod = 2), adapt_delta = 0.95, init = 0.5)
  expect_equal(rec$sample$adapt_delta, 0.95)
  expect_equal(rec$sample$init, 0.5)    # user init works for non-GRM
})
