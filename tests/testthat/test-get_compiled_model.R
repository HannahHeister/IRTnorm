make_fake_model <- function(id, exe) {
  force(exe)
  structure(list(id = id, exe_file = function() exe), class = "fake_model")
}

setup_cache <- function(env = parent.frame()) {
  testthat::skip_if_not_installed("cmdstanr")   # only the R package is needed, not CmdStan itself
  
  rec <- new.env(); rec$compiled <- 0; rec$code <- NULL
  dir <- withr::local_tempdir(.local_envir = env)
  
  local_mocked_bindings(
    get_cache_dir = function() dir,
    create_stan_data = function(age_model, irt_model, prior_knowledge, parameter_fixed) {
      rec$code <- list(age = age_model, irt = irt_model,
                       prior = prior_knowledge, fixed = parameter_fixed)
      "parameters {}"
    },
    .package = PKG, .env = env)
  
  local_mocked_bindings(
    cmdstan_model = function(stan_file, ...) {
      rec$compiled <- rec$compiled + 1
      exe <- file.path(dir, paste0("fake_exe_", rec$compiled))
      file.create(exe)
      make_fake_model(rec$compiled, exe)
    },
    .package = "cmdstanr", .env = env)
  
  old <- as.list(.model_cache, all.names = TRUE)
  rm(list = ls(.model_cache, all.names = TRUE), envir = .model_cache)
  withr::defer({
    rm(list = ls(.model_cache, all.names = TRUE), envir = .model_cache)
    list2env(old, envir = .model_cache)
  }, envir = env)
  
  rec$dir <- dir
  rec
}
# test-get_compiled_model.R

test_that("a model is compiled once, then served from memory", {
  rec <- setup_cache()
  m1 <- get_compiled_model("poly", "2PLnorm")
  m2 <- get_compiled_model("poly", "2PLnorm")
  expect_equal(rec$compiled, 1)
  expect_identical(m1, m2)
})

test_that("a fresh session loads from the disk cache instead of recompiling", {
  rec <- setup_cache()
  get_compiled_model("poly", "2PLnorm")                       # compiles, writes .rds
  expect_true(file.exists(file.path(rec$dir, "poly_2PLnorm.rds")))
  
  rm(list = ls(.model_cache, all.names = TRUE), envir = .model_cache)   # simulate new session
  expect_message(m <- get_compiled_model("poly", "2PLnorm"), "Loading cached model")
  expect_equal(rec$compiled, 1)                               # not compiled again
  expect_true(exists("poly_2PLnorm", envir = .model_cache))   # memory cache repopulated
})

test_that("different model configurations never share a cache entry", {
  rec <- setup_cache()
  calls <- list(
    list("poly",    "2PLnorm"),
    list("splines", "2PLnorm"),
    list("poly",    "1PLnorm"),
    list("poly",    "2PLnorm", "alpha"),                      # prior
    list("poly",    "2PLnorm", NULL, "alpha"),                # fixed
    list("poly",    "2PLnorm", c("alpha", "beta")),
    list("poly",    "2PLnorm", NULL, c("alpha", "beta")))
  for (a in calls) do.call(get_compiled_model, a)
  expect_equal(rec$compiled, length(calls))                   # every one compiled separately
  expect_length(list.files(rec$dir, pattern = "\\.rds$"), length(calls))
})

test_that("arguments reach create_stan_data()", {
  rec <- setup_cache()
  get_compiled_model("splines", "3PLnorm", prior_knowledge = "alpha", parameter_fixed = "beta")
  expect_equal(rec$code, list(age = "splines", irt = "3PLnorm", prior = "alpha", fixed = "beta"))
})

test_that("the order of prior_knowledge / parameter_fixed does not change the cache key", {
  rec <- setup_cache()
  get_compiled_model("poly", "2PLnorm", prior_knowledge = c("alpha", "beta"))
  get_compiled_model("poly", "2PLnorm", prior_knowledge = c("beta", "alpha"))
  expect_equal(rec$compiled, 1)   # fails until the key uses sort()
})