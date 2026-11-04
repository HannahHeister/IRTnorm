# tests/testthat/test-compare_age_models.R
#
# Prerequisite: add this guard at the top of compare_age_models(), before the
# hierarchical-model check:
#
#   if (!is.character(irt_model) || length(irt_model) != 1 || is.na(irt_model)) {
#     stop("`irt_model` must be a single model name (a character string of length 1).")
#   }
#
# Tests marked "FAILS TODAY" are specifications for known bugs. They should
# fail on the current code and pass once the bug is fixed.

# ============================================================
# Helpers
# ============================================================

# Fakes that record the arguments of EVERY fit call (poly rows and splines).
# NOTE: the return shapes are inferred from how compare_age_models() uses them.
# Test 4 (integration) is what verifies them against the real helpers.
make_recording_fakes <- function(elpd_fun = function(mu, sd) -100 - mu - sd,
                                 splines_elpd = -150, n_obs = 20) {
  rec <- new.env()
  rec$calls <- list()   # normalised arguments of every fit
  rec$rows  <- list()   # (mu, sd) of every poly fit, in order
  
  common <- c("raw_data", "age_variable", "item_variables", "seed",
              "iter_warmup", "iter_sampling", "chains", "parallel_chains")
  
  fit_one_row <- function(row, raw_data, irt_model, age_variable, item_variables,
                          prior_knowledge, parameter_fixed, seed, iter_warmup,
                          iter_sampling, chains, parallel_chains) {
    a <- mget(c(common, "irt_model", "prior_knowledge", "parameter_fixed"))
    rec$calls[[length(rec$calls) + 1]] <- c(list(kind = "poly"), a)
    rec$rows[[length(rec$rows) + 1]]   <- data.frame(mu = row$mu, sd = row$sd)
    e  <- elpd_fun(row$mu, row$sd)
    id <- paste0("poly_", row$mu, "_", row$sd)
    list(
      elpd_pointwise = rep(e / n_obs, n_obs),
      elpd_summary   = data.frame(modelcombi = id, poly_mean = row$mu,
                                  poly_sd = row$sd, ELPD = e),
      trajectory     = data.frame(modelcombi = id, age = 1:3, value = 0:2),
      paretok        = rep(0.2, n_obs)
    )
  }
  
  fit_one_model <- function(model_specifications, raw_data, age_variable,
                            item_variables, traj_info, iter_warmup, iter_sampling,
                            chains, parallel_chains, seed) {
    a <- mget(common)
    a$irt_model       <- model_specifications$irt_model
    a$prior_knowledge <- model_specifications$prior_knowledge
    a$parameter_fixed <- model_specifications$parameter_fixed
    rec$calls[[length(rec$calls) + 1]] <- c(list(kind = "splines"), a)
    list(
      elpd_pointwise = rep(splines_elpd / n_obs, n_obs),
      elpd_summary   = data.frame(modelcombi = "splines", poly_mean = NA_real_,
                                  poly_sd = NA_real_, ELPD = splines_elpd),
      trajectory     = data.frame(modelcombi = "splines", age = 1:3, value = 0:2),
      paretok        = rep(0.1, n_obs)
    )
  }
  
  # Ranks purely by ELPD (a stand-in for the real loo_comp)
  loo_comp <- function(loo_info) {
    ev <- loo_info$elpd_value
    ev <- ev[order(-ev$ELPD), ]
    ev$rank <- seq_len(nrow(ev))
    ev$mu   <- ev$poly_mean
    ev$sd   <- ev$poly_sd
    ev
  }
  
  list(fit_one_row = fit_one_row, fit_one_model = fit_one_model,
       loo_comp = loo_comp, rec = rec)
}

# Install the fakes for the calling test
local_fakes <- function(fk, env = parent.frame()) {
  local_mocked_bindings(fit_one_row   = fk$fit_one_row,
                        fit_one_model = fk$fit_one_model,
                        loo_comp      = fk$loo_comp,
                        .env = env)
}

# Mocks that error with a distinctive message if fitting starts. Any test
# expecting a validation error must match the validation message, so it cannot
# pass because of this mock's own error.
no_fit <- function(...) stop("FITTING STARTED")

local_no_fitting <- function(env = parent.frame()) {
  local_mocked_bindings(fit_one_row = no_fit, fit_one_model = no_fit,
                        loo_comp = no_fit, .env = env)
}

call_cam <- function(...) {
  compare_age_models(raw_data = data.frame(), age_variable = "age",
                     item_variables = "i1", parallel_chains = 1, seed = 1, ...)
}

one_row_grid <- data.frame(mu = 1, sd = 0)

evaluated <- function(fk) {
  do.call(rbind, fk$rec$rows)
}

# ============================================================
# 1. Forward-selection control flow
# ============================================================

test_that("forward: stops immediately when the start model is a local optimum", {
  # Peak at (1, 0): both neighbours are worse
  fk <- make_recording_fakes(function(mu, sd) -100 - 10 * (mu - 1) - 10 * sd)
  local_fakes(fk)
  
  run <- call_cam(method = "forward", irt_model = "1PLnorm",
                  include_splines = FALSE)
  
  fitted <- evaluated(fk)
  expect_equal(nrow(fitted), 3)                       # start + 2 candidates
  expect_equal(fitted[1, ], data.frame(mu = 1, sd = 0), ignore_attr = TRUE)
  expect_setequal(paste(fitted$mu, fitted$sd), c("1 0", "2 0", "1 1"))
  expect_equal(run$best_configuration$poly_mean, 1)
  expect_equal(run$best_configuration$poly_sd, 0)
})

test_that("forward: follows the better neighbour and stops at a local optimum", {
  # Distinct values (no ties); peak at (2, 1)
  surface <- function(mu, sd) -10 * abs(mu - 2) - 3 * abs(sd - 1) - 0.01 * (mu + sd)
  fk <- make_recording_fakes(surface)
  local_fakes(fk)
  
  run <- call_cam(method = "forward", irt_model = "1PLnorm",
                  maxeval = 50, include_splines = FALSE)
  
  fitted <- evaluated(fk)
  
  # Exact search path: start, then (mu+1, sd) and (mu, sd+1) at each step
  expect_equal(
    paste(fitted$mu, fitted$sd),
    c("1 0",
      "2 0", "1 1",     # from (1,0): moves to (2,0)
      "3 0", "2 1",     # from (2,0): moves to (2,1)
      "3 1", "2 2")     # from (2,1): neither improves -> stop
  )
  expect_equal(anyDuplicated(fitted[, c("mu", "sd")]), 0)   # never refits
  
  expect_equal(run$best_configuration$poly_mean, 2)
  expect_equal(run$best_configuration$poly_sd, 1)
  expect_equal(run$best_configuration$age_model, "poly")
})

# ============================================================
# 2. Specification tests for known bugs
# ============================================================

test_that("forward: reports splines when splines beat every polynomial  [FAILS TODAY]", {
  # Polys are terrible (-500ish); fake splines ELPD is -150.
  # finalize() only switches to "splines" when best_configuration is NULL,
  # but the forward path always supplies it.
  fk <- make_recording_fakes(function(mu, sd) -500, splines_elpd = -150)
  local_fakes(fk)
  
  run <- call_cam(method = "forward", irt_model = "1PLnorm",
                  include_splines = TRUE)
  
  expect_equal(run$best_model$modelcombi[run$best_model$rank == 1], "splines")
  expect_equal(run$best_configuration$age_model, "splines")
})

test_that("grid: reports splines when splines beat every polynomial", {
  # Counterpart of the test above; this path already works.
  fk <- make_recording_fakes(function(mu, sd) -500, splines_elpd = -150)
  local_fakes(fk)
  
  run <- call_cam(method = "grid", grid = one_row_grid, irt_model = "1PLnorm",
                  include_splines = TRUE)
  
  expect_equal(run$best_configuration$age_model, "splines")
})

test_that("forward: never fits more than maxeval polynomial models  [FAILS TODAY]", {
  # ELPD always improves, so only maxeval can stop the search.
  # Each iteration adds 2 candidates and the loop checks n_eval < maxeval,
  # so even maxeval values overshoot by one (e.g. 10 -> 11 fits).
  for (me in c(4, 5, 10)) {
    fk <- make_recording_fakes(function(mu, sd) mu + sd)
    local_fakes(fk)
    
    call_cam(method = "forward", irt_model = "1PLnorm",
             maxeval = me, include_splines = FALSE)
    
    expect_lte(length(fk$rec$rows), me)
  }
})

# ============================================================
# 3. Comparability of fits
# ============================================================

test_that("all candidate fits, including splines, share identical data and settings", {
  dat <- data.frame(age = c(1, 2, 3), i1 = c(0, 1, 1))
  pk  <- list(a = 1)
  pf  <- list(b = 2)
  
  for (method in c("grid", "forward")) {
    fk <- make_recording_fakes()
    local_fakes(fk)
    
    compare_age_models(
      method = method,
      grid = if (method == "grid") data.frame(mu = c(1, 2, 3), sd = c(0, 1, 0)),
      raw_data = dat, irt_model = "1PLnorm",
      age_variable = "age", item_variables = "i1",
      prior_knowledge = pk, parameter_fixed = pf,
      seed = 123, iter_warmup = 77, iter_sampling = 88,
      chains = 3, parallel_chains = 2, include_splines = TRUE
    )
    
    calls <- fk$rec$calls
    expect_gt(length(calls), 2)
    expect_true("splines" %in% vapply(calls, `[[`, "", "kind"), info = method)
    
    # Strip the 'kind' label and compare every call against the first
    settings <- lapply(calls, function(x) x[setdiff(names(x), "kind")])
    for (i in seq_along(settings)[-1]) {
      expect_identical(settings[[i]], settings[[1]],
                       info = sprintf("%s: call %d differs from call 1", method, i))
    }
    
    # And they are the values we supplied, not silently replaced defaults
    expect_equal(settings[[1]]$seed, 123)
    expect_equal(settings[[1]]$iter_warmup, 77)
    expect_equal(settings[[1]]$iter_sampling, 88)
    expect_equal(settings[[1]]$chains, 3)
    expect_equal(settings[[1]]$parallel_chains, 2)
    expect_equal(settings[[1]]$prior_knowledge, pk)
    expect_equal(settings[[1]]$parameter_fixed, pf)
  }
})

test_that("a NULL seed is drawn once and reused for every fit", {
  fk <- make_recording_fakes()
  local_fakes(fk)
  
  compare_age_models(method = "grid", grid = data.frame(mu = 1:3, sd = 0),
                     raw_data = data.frame(), irt_model = "1PLnorm",
                     age_variable = "age", item_variables = "i1",
                     seed = NULL, parallel_chains = 1, include_splines = TRUE)
  
  seeds <- vapply(fk$rec$calls, function(x) x$seed, numeric(1))
  expect_length(seeds, 4)           # 3 poly + splines
  expect_length(unique(seeds), 1)
})

# ============================================================
# 5. Input validation
# ============================================================

test_that("each hierarchical/testlet model is rejected before any fitting", {
  local_no_fitting()
  hier <- c("Testlet1PLnorm", "Testlet2PLnorm", "HO1PLnorm-C", "HO2PLnorm-C",
            "HO1PLnorm-P", "HO2PLnorm-P", "HO1PLnorm-F", "HO2PLnorm-F")
  for (m in hier) {
    expect_error(call_cam(method = "grid", grid = one_row_grid, irt_model = m),
                 "unidimensional", info = m)
  }
})

test_that("a vector irt_model is rejected outright, before fitting", {
  local_no_fitting()
  vectors <- list(
    c("1PLnorm", "2PLnorm"),            # two valid models
    c("1PLnorm", "Testlet1PLnorm"),     # one valid, one hierarchical
    c("Testlet1PLnorm", "1PLnorm"),     # order must not matter
    c("2PLnorm", "HO2PLnorm-F")
  )
  for (v in vectors) {
    expect_error(call_cam(method = "grid", grid = one_row_grid, irt_model = v),
                 "length 1",
                 info = paste(v, collapse = " + "))
  }
})

test_that("a vector irt_model is rejected for the forward method too", {
  local_no_fitting()
  expect_error(call_cam(method = "forward", irt_model = c("1PLnorm", "2PLnorm")),
               "length 1")
})

test_that("empty, NA, unknown and non-character irt_model values are rejected", {
  local_no_fitting()
  
  # Wrong length
  expect_error(call_cam(method = "grid", grid = one_row_grid,
                        irt_model = character(0)),
               "length 1")
  
  # Not in the allowed set (NA and unknown names)
  for (b in list(NA_character_, "4PLnorm", "")) {
    expect_error(call_cam(method = "grid", grid = one_row_grid, irt_model = b),
                 "should be one of",
                 info = deparse(b))
  }
  
  # Non-character input: wording varies by R version, so only require that
  # it fails and that fitting never started.
  for (b in list(1, factor("1PLnorm"), TRUE)) {
    err <- tryCatch(
      call_cam(method = "grid", grid = one_row_grid, irt_model = b),
      error = function(e) conditionMessage(e)
    )
    expect_type(err, "character")
    expect_false(grepl("FITTING STARTED", err), info = deparse(b))
  }
})

test_that("unidimensional models pass validation and reach the fitting stage", {
  local_no_fitting()
  for (m in c("1PLnorm", "2PLnorm")) {
    expect_error(call_cam(method = "grid", grid = one_row_grid, irt_model = m,
                          include_splines = FALSE),
                 "FITTING STARTED", info = m)
  }
})

test_that("grid method without a grid errors before fitting", {
  local_no_fitting()
  expect_error(call_cam(method = "grid", irt_model = "1PLnorm"),
               "`grid` must be supplied")
})

test_that("the default method is 'grid', so omitting both method and grid errors", {
  local_no_fitting()
  expect_error(call_cam(irt_model = "1PLnorm"), "`grid` must be supplied")
})

test_that("forward method needs no grid and reaches the fitting stage", {
  local_no_fitting()
  expect_error(call_cam(method = "forward", irt_model = "1PLnorm",
                        include_splines = FALSE),
               "FITTING STARTED")
})

test_that("an unknown method is rejected before fitting", {
  # Kept only as an ordering check: method is validated before anything runs.
  local_no_fitting()
  expect_error(call_cam(method = "random", grid = one_row_grid,
                        irt_model = "1PLnorm"),
               "should be one of")
})


