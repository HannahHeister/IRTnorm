# ---- helpers (tests/testthat/helper-compare.R) ----------------------------
spec <- function(irt = "1PLnorm", age = "poly", ...) {
  c(list(irt_model = irt, age_model = age), list(...))
}

make_model_spec <- function(irt_model = "2PLnorm",
                            age_model = "poly",
                            poly_mean = 3,
                            poly_sd = 2) {
  list(
    irt_model = irt_model,
    age_model = age_model,
    poly_mean = poly_mean,
    poly_sd = poly_sd
  )
}

make_test_data <- function() {
  data.frame(
    age = c(10, 11, 12),
    item1 = c(0, 1, 1),
    item2 = c(1, 0, 1)
  )
}
n_person  <- 4
n_item    <- 2
n_obs     <- n_person * n_item
item_vars <- paste0("item", seq_len(n_item))
raw <- data.frame(age = c(1, 2, 3, 4), item1 = c(0, 1, 0, 1), item2 = c(1, 1, 0, 0))

# Per-person ELPD contributions: 1PL = -6, 2PL = -4 (best), 3PL = -5, GRM = -8
elpd_per_person <- list(
  "1PLnorm" = c(-2, -1, -2, -1),
  "2PLnorm" = c(-1, -1, -1, -1),
  "3PLnorm" = c(-1, -2, -1, -1),
  "GRMnorm" = c(-2, -2, -2, -2)
)
# Person-major layout: person 1 item 1, person 1 item 2, person 2 item 1, ...
pw <- lapply(elpd_per_person, function(x) rep(x / n_item, each = n_item))
stopifnot(all(vapply(pw, length, integer(1)) == n_obs))

fake_fit_result <- function(spec, pointwise) {
  list(
    elpd_pointwise = pointwise,
    elpd_summary = data.frame(
      ELPD       = sum(pointwise),
      irt_model  = spec$irt_model,
      age_model  = spec$age_model,
      poly_mean  = if (is.null(spec$poly_mean)) NA else spec$poly_mean,
      poly_sd    = if (is.null(spec$poly_sd)) NA else spec$poly_sd,
      modelcombi = paste0(spec$irt_model, "_", spec$age_model)
    ),
    paretok = rep(0.1, length(pointwise)),
    neff    = rep(500, length(pointwise)),
    reff    = rep(1, length(pointwise))
  )
}

local_fake_fits <- function(fail = character(), .env = parent.frame()) {
  calls <- new.env()
  calls$args <- list()
  testthat::local_mocked_bindings(
    fit_one_model = function(model_specifications, ...) {
      calls$args[[length(calls$args) + 1]] <-
        c(list(model_specifications = model_specifications), list(...))
      irt <- model_specifications$irt_model
      if (irt %in% fail) stop("convergence issues")
      fake_fit_result(model_specifications, pw[[irt]])
    },
    .env = .env
  )
  calls
}

run_compare <- function(specs, ...) {
  compare_IRTnorm_models(specs, raw_data = raw, age_variable = "age",
                         item_variables = item_vars, ...)
}

three_specs <- list(spec("1PLnorm", poly_mean = 1, poly_sd = 0),
                    spec("2PLnorm", poly_mean = 1, poly_sd = 0),
                    spec("3PLnorm", poly_mean = 1, poly_sd = 0))

# ---- input checks ----------------------------------------------------------

test_that("invalid specification lists are rejected before any fitting", {
  calls <- local_fake_fits()
  expect_error(run_compare(list(spec("1PLnorm"))), "at least two")
  expect_length(calls$args, 0)
})

test_that("missing columns in raw_data raise an informative error", {
  local_fake_fits()
  expect_error(
    compare_IRTnorm_models(three_specs, raw_data = raw, age_variable = "age",
                           item_variables = c("item1", "nope")),
    "Columns not found in raw_data: nope"
  )
  expect_error(
    compare_IRTnorm_models(three_specs, raw_data = raw, age_variable = "oops",
                           item_variables = c("item1", "item2")),
    "Columns not found in raw_data: oops"
  )
})
test_that("does not run if not all required model aspects are definded", {
  local_fake_fits()
  specs <- list(spec("1PLnorm", "poly", poly_mean = 2, poly_sd = 3),
                spec("2PLnorm", "poly"))
  expect_error(run_compare(specs),"missing poly_mean, poly_sd")
  
})
# ---- output ----------------------------------------------------------------

test_that("output contains best_model and diagnostics", {
  local_fake_fits()
  out <- run_compare(three_specs)
  expect_named(out, c("best_model", "diagnostics"))
})


test_that("fake fits have one pointwise value per person-item response", {
  # guards the stub itself, so a wrong stub fails here and not in a confusing place
  res <- fake_fit_result(spec("2PLnorm"), pw[["2PLnorm"]])
  expect_length(res$elpd_pointwise, n_person * n_item)
  expect_length(res$paretok, n_person * n_item)
})

test_that("diagnostics has one row per model, item and person", {
  local_fake_fits()
  diag <- run_compare(three_specs)$diagnostics
  expect_equal(nrow(diag), length(three_specs) * n_item * n_person)   # 3 * 2 * 4 = 24
  expect_equal(as.vector(table(diag$irt_model)), rep(n_obs, length(three_specs)))
  expect_setequal(diag$irt_model, c("1PLnorm", "2PLnorm", "3PLnorm"))
  expect_true(all(c("model", "irt_model", "age_model", "poly_mean", "poly_sd",
                    "pareto_k", "neff", "reff") %in% names(diag)))
})



test_that("models are ranked by ELPD and labelled by position and irt_model", {
  local_fake_fits()
  out <- run_compare(three_specs)$best_model
  expect_equal(rownames(out), c("2_2PLnorm", "3_3PLnorm", "1_1PLnorm"))
  #expect_equal(out$ELPD, c(-4, -5, -6))
  #expect_equal(out$ELPD_diff, c(0, -1, -2))
  expect_equal(out$rank, 1:3)
  expect_equal(out$SE[1], 0)
  
  d <- pw[["2PLnorm"]] - pw[["3PLnorm"]]
  expect_equal(out$SE[2], sqrt(length(d) * var(d)))   # = sqrt(3/7) with 8 values
})


test_that("diagnostics has one row per model and observation, labelled correctly", {
  local_fake_fits()
  diag <- run_compare(three_specs)$diagnostics
  expect_equal(nrow(diag), 3 * 4*2)
  expect_setequal(diag$irt_model, rep(c("1PLnorm", "2PLnorm", "3PLnorm"), each = 4*2))
  expect_true(all(c("model", "irt_model", "age_model", "poly_mean", "poly_sd",
                    "pareto_k", "neff", "reff") %in% names(diag)))
})

test_that("poly settings are reported for poly models and NA otherwise", {
  local_fake_fits()
  specs <- list(spec("1PLnorm", "poly", poly_mean = 2, poly_sd = 3),
                spec("2PLnorm", "splines", k_sd = 2, k_mean =3))             # "spline" is a placeholder
  diag <- run_compare(specs)$diagnostics
  d1 <- diag[diag$model == "1_1PLnorm", ]
  d2 <- diag[diag$model == "2_2PLnorm", ]
  expect_true(all(d1$poly_mean == 2 & d1$poly_sd == 3))
  expect_true(all(is.na(d2$poly_mean) & is.na(d2$poly_sd)))
})

# ---- failure handling ------------------------------------------------------

test_that("a failed model is dropped with a warning and the rest are compared", {
  local_fake_fits(fail = "2PLnorm")
  expect_warning(out <- run_compare(three_specs), "Model 2_2PLnorm failed: convergence issues")
  expect_equal(sort(rownames(out$best_model)), c("1_1PLnorm", "3_3PLnorm"))
  expect_false("2_2PLnorm" %in% out$diagnostics$model)
})

test_that("specs and labels stay aligned after a failure", {
  local_fake_fits(fail = "1PLnorm")      # the FIRST model fails
  expect_warning(out <- run_compare(three_specs), "1_1PLnorm failed")
  diag <- out$diagnostics
  expect_equal(unique(diag$irt_model[diag$model == "2_2PLnorm"]), "2PLnorm")
  expect_equal(unique(diag$irt_model[diag$model == "3_3PLnorm"]), "3PLnorm")
})

test_that("fewer than two successful fits is an error", {
  local_fake_fits(fail = c("2PLnorm", "3PLnorm"))
  suppressWarnings(
    expect_error(run_compare(three_specs), "Fewer than two models fitted successfully")
  )
})

# ---- argument handling -----------------------------------------------------

test_that("every model gets the same seed, generated when none is supplied", {
  calls <- local_fake_fits()
  run_compare(three_specs)
  seeds <- vapply(calls$args, function(a) a$seed, numeric(1))
  expect_length(unique(seeds), 1)
  expect_false(is.na(seeds[1]))
})

test_that("a supplied seed is passed through unchanged", {
  calls <- local_fake_fits()
  run_compare(three_specs, seed = 42)
  expect_true(all(vapply(calls$args, function(a) a$seed == 42, logical(1))))
})

test_that("parallel_chains defaults to at most the number of chains", {
  calls <- local_fake_fits()
  run_compare(three_specs, chains = 2)
  pc <- vapply(calls$args, function(a) a$parallel_chains, numeric(1))
  expect_true(all(pc >= 1 & pc <= 2))
})

test_that("an explicit parallel_chains is respected", {
  calls <- local_fake_fits()
  run_compare(three_specs, chains = 4, parallel_chains = 3)
  expect_true(all(vapply(calls$args, function(a) a$parallel_chains == 3, logical(1))))
})

test_that("data and sampling settings are forwarded to every fit", {
  calls <- local_fake_fits()
  run_compare(three_specs, iter_warmup = 10, iter_sampling = 20, chains = 3,
              parallel_chains = 1, seed = 1)
  expect_length(calls$args, 3)
  for (a in calls$args) {
    expect_equal(a$raw_data, raw)
    expect_equal(a$age_variable, "age")
    expect_equal(a$item_variables, c("item1", "item2"))
    expect_equal(a$iter_warmup, 10)
    expect_equal(a$iter_sampling, 20)
    expect_equal(a$chains, 3)
  }
})

test_that("each model receives its own specification", {
  calls <- local_fake_fits()
  run_compare(three_specs)
  irts <- vapply(calls$args, function(a) a$model_specifications$irt_model, character(1))
  expect_equal(irts, c("1PLnorm", "2PLnorm", "3PLnorm"))
})

test_that("missing columns in raw_data produce an error", {
  local_mocked_bindings(
    validate_model_specifications_list = function(x) NULL,
    check_dimension = function(x) NULL, .package = "IRTnorm"
  )
  
  expect_error(
    compare_IRTnorm_models(
      model_specifications = list(
        make_model_spec(),
        make_model_spec("3PLnorm")
      ),
      raw_data = make_test_data(),
      age_variable = "age",
      item_variables = c("item1", "missing_item")
    ),
    "Columns not found in raw_data: missing_item"
  )
})

test_that("missing age column produces an error", {
  local_mocked_bindings(
    validate_model_specifications_list = function(x) NULL,
    check_dimension = function(x) NULL
  )
  
  expect_error(
    compare_IRTnorm_models(
      model_specifications = list(
        make_model_spec(),
        make_model_spec("3PLnorm")
      ),
      raw_data = make_test_data(),
      age_variable = "missing_age",
      item_variables = c("item1", "item2")
    ),
    "Columns not found in raw_data: missing_age"
  )
})

