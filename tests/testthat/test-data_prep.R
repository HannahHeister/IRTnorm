library(testthat)

PKG <- "IRTnorm"

# ---- Fixture: 6 persons, 4 items, a few NAs ----------------------------------
raw <- data.frame(
  age = c(5, 6, 7, 8, 9, 10),
  i1  = c(1, 0, 1, NA, 1, 0),
  i2  = c(0, 1, NA, 1, 0, 1),
  i3  = c(1, 1, 0, 0, NA, 1),
  i4  = c(0, NA, 1, 1, 1, 0)
)
items <- paste0("i", 1:4)
I <- 4; J <- 6
itemD4 <- c(1, 1, 2, 2)

prep <- function(irt_model = "2PLnorm", age_model = "poly", ...) {
  args <- list(raw_data = raw, age_variable = "age", item_variables = items,
               irt_model = irt_model, age_model = age_model)
  if (age_model == "poly") args <- c(args, list(poly_mean = 2, poly_sd = 1))
  args <- modifyList(args, list(...))
  do.call(data_prep, args)
}

fake_spline <- function(age_df, k, fixed_knots) {
  list(S_fixed = 2, K = k - 2, X = matrix(1, nrow(age_df), 2),
       Z = matrix(2, nrow(age_df), k - 2), R_Z = matrix(3, k, k - 2),
       R_X = matrix(4, k, 2), sm = paste0("sm_k", k))
}

# ==============================================================================
# Core data layout
# ==============================================================================
test_that("long format: ii, jj, y are aligned, NA-free and item-major", {
  d <- prep()
  n_obs <- sum(!is.na(raw[, items]))
  expect_equal(d$N, n_obs)
  expect_length(d$ii, n_obs); expect_length(d$jj, n_obs); expect_length(d$y, n_obs)
  expect_false(anyNA(d$y))
  expect_equal(d$I, I); expect_equal(d$J, J)
  # the key invariant: every (jj, ii) pair points at the matching response
  expect_equal(d$y, as.numeric(as.matrix(raw[, items])[cbind(d$jj, d$ii)]))
  expect_false(is.unsorted(d$ii))   # item-major
})

test_that("age is centred/scaled with rounded mean and sd", {
  d <- prep()
  expect_equal(d$m_age, round(mean(raw$age), 2))
  expect_equal(d$s_age, round(sd(raw$age), 2))
  expect_equal(d$age_variable, raw$age)
})

test_that("unidimensional defaults: itemD, dd, no multidim fields", {
  d <- prep("1PLnorm")
  expect_equal(d$itemD, rep(1, I))
  expect_equal(d$dd, rep(1, d$N))
  expect_null(d$D); expect_null(d$JD_loadings)
})

test_that("irt_model strings map to the right irt_mod number", {
  uni <- c("1PLnorm" = 1, "2PLnorm" = 2, "3PLnorm" = 3)
  for (nm in names(uni)) expect_equal(prep(nm)$irt_mod, uni[[nm]], info = nm)
  
  multi <- c("Testlet1PLnorm" = 5, "Testlet2PLnorm" = 6, "HO1PLnorm-C" = 7,
             "HO2PLnorm-C" = 8, "HO1PLnorm-P" = 9, "HO2PLnorm-P" = 10,
             "HO1PLnorm-F" = 11, "HO2PLnorm-F" = 12)
  for (nm in names(multi))
    expect_equal(prep(nm, itemD = itemD4)$irt_mod, multi[[nm]], info = nm)
  
  expect_equal(suppressWarnings(prep("GRMnorm", K_item = rep(2, I)))$irt_mod, 4)
  expect_error(prep("NotAModel"))
})

# ==============================================================================
# Multidimensional models
# ==============================================================================
test_that("multidim fields: D, U, dd and JD_loadings", {
  d <- prep("HO1PLnorm-C", itemD = itemD4)
  expect_equal(d$D, 2); expect_equal(d$U, 2 * J)
  expect_equal(d$itemD, itemD4)
  expect_equal(d$dd, itemD4[d$ii])                 # domain per observation
  # loading index = person + (domain - 1) * J
  expect_equal(d$JD_loadings, d$jj + (d$dd - 1) * J)
  expect_true(all(d$JD_loadings >= 1 & d$JD_loadings <= d$U))
})

test_that("multidim models validate itemD", {
  expect_error(prep("HO1PLnorm-C"), "itemD missing")
  expect_error(prep("HO1PLnorm-C", itemD = c(1, 2)), "differ")
})

# ==============================================================================
# Polynomial age model
# ==============================================================================
test_that("poly: design matrices and dimensions", {
  d <- prep(poly_mean = 2, poly_sd = 1)
  expect_equal(d$age_mod, 1)
  expect_equal(d$S, 3); expect_equal(d$W, 2)
  expect_equal(dim(d$age), c(J, 3)); expect_equal(dim(d$age2), c(J, 2))
  expect_equal(unname(d$age[, 1]), rep(1, J))                  # intercept
  scaled <- (raw$age - d$m_age) / d$s_age
  expect_equal(unname(d$age[, 2]), scaled)
  expect_equal(unname(d$age[, 3]), scaled^2)
})

test_that("poly degree 0 gives an intercept-only matrix", {
  d <- prep(poly_mean = 0, poly_sd = 0)
  expect_equal(d$S, 1); expect_equal(d$W, 1)
  expect_equal(d$age, matrix(1, J, 1))
})

test_that("poly: unique-sd design only for P/F models and defaults to the sd design", {
  expect_null(prep("2PLnorm")$age3)
  expect_null(prep("HO1PLnorm-C", itemD = itemD4)$age3)
  
  d <- prep("HO1PLnorm-P", itemD = itemD4)                      # poly_sd_unique NULL
  expect_equal(d$age3, d$age2); expect_equal(d$M, d$W)
  
  d2 <- prep("HO2PLnorm-F", itemD = itemD4, poly_sd_unique = 2)
  expect_equal(d2$M, 3); expect_equal(ncol(d2$age3), 3)
})

test_that("poly with missing poly_mean/poly_sd gives an informative error", {
  # FAILS with current code: `if (poly_sd == 0)` errors with a cryptic
  # "argument is of length zero". Add an explicit check.
  expect_error(prep(poly_mean = NULL, poly_sd = NULL), "poly_mean|poly_sd")
})

# ==============================================================================
# Spline age model (build_spline mocked: tests the wiring only)
# ==============================================================================
test_that("splines: components are mapped to the right names", {
  local_mocked_bindings(build_spline = fake_spline, .package = PKG)
  d <- prep(age_model = "splines", k_mean = 7, k_sd = 5)
  expect_equal(d$age_mod, 2)
  expect_equal(d$K1, 5); expect_equal(d$K2, 3)                  # k - 2
  expect_equal(d$sm1, "sm_k7"); expect_equal(d$sm2, "sm_k5")
  expect_equal(dim(d$R1_Z), c(7, 5)); expect_equal(dim(d$R2_Z), c(5, 3))
  expect_true(all(d$X1 == 1) && all(d$Z1 == 2) && all(d$R1_X == 4))
  expect_null(d$sm3); expect_null(d$K3)
})

test_that("splines: third spline only for P/F models and uses k_sd2", {
  local_mocked_bindings(build_spline = fake_spline, .package = PKG)
  for (m in c("HO1PLnorm-P", "HO2PLnorm-P", "HO1PLnorm-F", "HO2PLnorm-F")) {
    d <- prep(m, "splines", itemD = itemD4, k_sd2 = 6)
    expect_equal(d$sm3, "sm_k6", info = m)
    expect_equal(d$K3, 4, info = m)
    expect_true(all(c("M_fixed", "X3", "Z3", "R3_Z", "R3_X") %in% names(d)), info = m)
  }
  expect_null(prep("HO1PLnorm-C", "splines", itemD = itemD4)$sm3)
})

test_that("splines: fixed knots and the standardised age reach build_spline", {
  seen <- list()
  local_mocked_bindings(
    build_spline = function(age_df, k, fixed_knots) {
      seen[[length(seen) + 1]] <<- list(age = age_df$age, k = k, fk = fixed_knots)
      fake_spline(age_df, k, fixed_knots)
    }, .package = PKG)
  d <- prep(age_model = "splines", k_mean = 7, k_sd = 5,
            fixed_knots_mean = "equidistant", fixed_knots_sd = c(-1, 0, 1))
  expect_equal(seen[[1]]$fk, "equidistant"); expect_equal(seen[[1]]$k, 7)
  expect_equal(seen[[2]]$fk, c(-1, 0, 1));   expect_equal(seen[[2]]$k, 5)
  expect_equal(seen[[1]]$age, (raw$age - d$m_age) / d$s_age)    # scaled, not raw
})

# ==============================================================================
# Prior knowledge and fixed parameters
# ==============================================================================
prior_alpha <- list(alpha = list(alpha_prior_mu = rep(1, I), alpha_prior_sd = rep(0.5, I)))

test_that("prior knowledge and fixed parameters are written into the data", {
  d <- prep(prior_knowledge = prior_alpha,
            parameter_fixed = list(beta = c(-1, 0, 1, 2)))
  expect_equal(d$alpha_prior_mu, rep(1, I))
  expect_equal(d$alpha_prior_sd, rep(0.5, I))
  expect_equal(d$beta, c(-1, 0, 1, 2))
})

test_that("a parameter cannot be both fixed and given a prior", {
  expect_error(prep(prior_knowledge = prior_alpha,
                    parameter_fixed = list(alpha = rep(1, I))),
               "only indicate one")
})

test_that("fixed alpha/beta/gamma must have one value per item", {
  expect_error(prep(parameter_fixed = list(alpha = c(1, 1))), "all items")
})

test_that("prior of wrong length for alpha/beta/gamma is rejected", {
  # FAILS with current code: the check combines `length(x) != 2` with `&`,
  # and each prior is a list of length 2, so the condition is never TRUE.
  bad <- list(alpha = list(alpha_prior_mu = c(1, 1), alpha_prior_sd = c(1, 1)))
  expect_error(prep(prior_knowledge = bad), "wrong format")
})

test_that("nu prior must have one value per domain for testlet/HO models", {
  bad <- list(nu = list(nu_prior_mu = c(0.5, 0.5, 0.5), nu_prior_sd = c(1, 1, 1)))
  expect_error(prep("HO1PLnorm-C", itemD = itemD4, prior_knowledge = bad), "D domains")
  ok <- list(nu = list(nu_prior_mu = c(0.5, 0.5), nu_prior_sd = c(1, 1)))
  expect_no_error(prep("HO1PLnorm-C", itemD = itemD4, prior_knowledge = ok))
})

# ==============================================================================
# GRM
# ==============================================================================
raw_grm <- data.frame(age = raw$age,
                      i1 = c(1, 2, 3, NA, 1, 2), i2 = c(2, 3, NA, 1, 2, 3),
                      i3 = c(1, 1, 2, 3, NA, 3), i4 = c(3, NA, 1, 2, 2, 1))

test_that("GRM sets K_max, K_item and irt_mod = 4; K_item is validated", {
  d <- data_prep(raw_grm, "age", items, "GRMnorm", "poly",
                 poly_mean = 1, poly_sd = 1, K_item = rep(3, I))
  expect_equal(d$irt_mod, 4); expect_equal(d$K_max, 3); expect_equal(d$K_item, rep(3, I))
  expect_error(data_prep(raw_grm, "age", items, "GRMnorm", "poly",
                         poly_mean = 1, poly_sd = 1), "K_item missing")
  expect_error(data_prep(raw_grm, "age", items, "GRMnorm", "poly",
                         poly_mean = 1, poly_sd = 1, K_item = c(3, 3)), "differ")
})

test_that("GRM responses coded from 0 are recoded to start at 1", {
  r0 <- raw_grm; r0[, items] <- r0[, items] - 1                 # now coded 0..2
  d <- suppressWarnings(data_prep(r0, "age", items, "GRMnorm", "poly",
                                  poly_mean = 1, poly_sd = 1, K_item = rep(3, I)))
  expect_equal(min(d$y), 1)
  expect_equal(max(d$y), 3)
})

test_that("GRM fixed thresholds must be an ordered list with one entry per item", {
  ok  <- list(b = lapply(1:I, function(i) c(-1, 1)))
  bad <- list(b = lapply(1:I, function(i) c(1, -1)))            # decreasing
  args <- list(raw_grm, "age", items, "GRMnorm", "poly",
               poly_mean = 1, poly_sd = 1, K_item = rep(3, I))
  expect_no_error(do.call(data_prep, c(args, list(parameter_fixed = ok))))
  expect_error(do.call(data_prep, c(args, list(parameter_fixed = bad))), "thresholds")
  expect_error(do.call(data_prep, c(args, list(parameter_fixed = list(b = c(-1, 1))))),
               "thresholds")
})