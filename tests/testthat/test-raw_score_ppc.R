# tests/testthat/test-raw_score_ppc.R

# ---- Helpers ---------------------------------------------------------------

# Mock IRTnorm fit. itemD gives the dimension of each item; K is the maximum
# score per item (K > 1 -> irt_mod 4).
make_ppc_fit <- function(J = 40, draws = 5, itemD = rep(1, 4), K = 1, seed = 1) {
  set.seed(seed)
  I <- length(itemD)
  multi <- any(itemD > 1)
  irt_mod <- if (multi) 5 else if (K > 1) 4 else 1
  item_variables <- paste0("item", seq_len(I))
  
  resp <- matrix(sample(0:K, J * I, replace = TRUE), J, I,
                 dimnames = list(NULL, item_variables))
  raw_data <- as.data.frame(resp)
  raw_data$age <- runif(J, 2.5, 11.5)            # strictly inside (2, 12)
  
  structure(list(
    newdata = list(itemD = itemD, D = max(itemD), irt_mod = irt_mod,
                   K_item = rep(K, I),
                   jj = rep(seq_len(J), each = I),   # person of each observation
                   dd = rep(itemD, times = J)),      # dimension of each observation
    model_info = list(y_rep = matrix(sample(0:K, draws * J * I, replace = TRUE),
                                     nrow = draws, ncol = J * I)),
    raw_data = raw_data,
    item_variables = item_variables,
    age_variable = "age"
  ), class = "IRTnorm")
}

set_age <- function(fit, person, value) {
  fit$raw_data$age[person] <- value
  fit
}

n_obs_total <- function(groups) sum(vapply(groups, function(g) sum(g$df$obs), numeric(1)))

# ---- Input validation ------------------------------------------------------

test_that("errors for non-IRTnorm input", {
  expect_error(raw_score_ppc(list()), "estimate_norming_model")
  expect_error(raw_score_ppc(NULL), "estimate_norming_model")
})

test_that("errors when age_range is not numeric", {
  expect_error(raw_score_ppc(make_ppc_fit(), age_range = c("2", "12")), "numeric")
})

# ---- Warning: individuals outside age_range --------------------------------

test_that("individuals outside age_range warn with count and bounds", {
  fit <- make_ppc_fit(J = 40)
  fit <- set_age(fit, 1, 1)      # below
  fit <- set_age(fit, 2, 13)     # above
  expect_warning(
    raw_score_ppc(fit, age_range = c(2, 6, 12)),
    "2 individual\\(s\\) fall outside 'age_range' \\(2,12\\]"
  )
})

test_that("lower bound is open and upper bound is closed", {
  fit <- set_age(make_ppc_fit(), 1, 2)               # exactly on lower bound
  expect_warning(raw_score_ppc(fit, c(2, 12)), "1 individual")
  
  fit <- set_age(make_ppc_fit(), 1, 12)              # exactly on upper bound
  expect_no_warning(raw_score_ppc(fit, c(2, 12)))
})

test_that("no warning when everyone is inside age_range or age_range is NULL", {
  expect_no_warning(raw_score_ppc(make_ppc_fit(), c(2, 6, 12)))
  expect_no_warning(raw_score_ppc(set_age(make_ppc_fit(), 1, 99)))   # NULL: no check
})

test_that("persons outside the range are excluded from the results", {
  fit <- set_age(set_age(make_ppc_fit(J = 40), 1, 1), 2, 13)
  res <- suppressWarnings(raw_score_ppc(fit, c(2, 6, 12)))
  expect_equal(n_obs_total(res[[1]]), 38)
})

test_that("the warning is raised once, also for multidimensional fits", {
  fit <- set_age(make_ppc_fit(itemD = c(1, 1, 2, 2)), 1, 1)
  w <- testthat::capture_warnings(raw_score_ppc(fit, c(2, 6, 12)))
  expect_length(w, 1)
})

# ---- Structure: unidimensional ---------------------------------------------

test_that("no age_range: one overall result with the expected fields", {
  fit <- make_ppc_fit(J = 40, draws = 5, itemD = rep(1, 4))
  res <- raw_score_ppc(fit)
  expect_length(res, 1)
  expect_equal(res[[1]]$age, "full age range")
  expect_equal(res[[1]]$dim, 0)
  expect_equal(res[[1]]$nitem, 4)
  expect_equal(res[[1]]$df$score, 0:4)
  expect_equal(dim(res[[1]]$replicated), c(5, 5))            # draws x (nitem + 1)
  expect_equal(unname(rowSums(res[[1]]$replicated)), rep(40, 5))
})

test_that("observed distribution equals the raw sum-score counts", {
  fit <- make_ppc_fit(J = 40, itemD = rep(1, 4))
  res <- raw_score_ppc(fit)[[1]]
  sums <- rowSums(fit$raw_data[fit$item_variables])
  expect_equal(res$df$obs, tabulate(sums + 1, nbins = 5))
})

test_that("missing item responses are ignored in the observed score", {
  fit <- make_ppc_fit(J = 40, itemD = rep(1, 4))
  fit$raw_data$item1[1:3] <- NA
  res <- raw_score_ppc(fit)[[1]]
  sums <- rowSums(fit$raw_data[fit$item_variables], na.rm = TRUE)
  expect_equal(res$df$obs, tabulate(sums + 1, nbins = 5))
  expect_equal(sum(res$df$obs), 40)
})

test_that("nitem is the number of items for dichotomous models", {
  expect_equal(raw_score_ppc(make_ppc_fit(itemD = rep(1, 6)))[[1]]$nitem, 6)
})

test_that("nitem is the sum of K_item for irt_mod 4", {
  res <- raw_score_ppc(make_ppc_fit(itemD = rep(1, 4), K = 3))   # max total score 12
  expect_equal(res[[1]]$nitem, 12)
  expect_equal(nrow(res[[1]]$df), 13)
  expect_equal(ncol(res[[1]]$replicated), 13)
})

# ---- Structure: age groups -------------------------------------------------

test_that("with age_range: one result per interval inside the dimension entry", {
  fit <- make_ppc_fit(J = 40)
  res <- raw_score_ppc(fit, c(2, 6, 12))
  expect_length(res, 1)                  # one dimension
  expect_length(res[[1]], 2)             # two age groups
  expect_equal(vapply(res[[1]], function(g) g$age, character(1)),
               c("age: (2,6]", "age: (6,12]"))
  expect_equal(n_obs_total(res[[1]]), 40)
})

test_that("each age group holds the right persons in observed and replicated data", {
  fit <- make_ppc_fit(J = 40)
  res <- raw_score_ppc(fit, c(2, 6, 12))[[1]]
  in_g1 <- fit$raw_data$age > 2 & fit$raw_data$age <= 6
  expect_equal(sum(res[[1]]$df$obs), sum(in_g1))
  expect_equal(unname(rowSums(res[[1]]$replicated)), rep(sum(in_g1), 5))
  
  sums <- rowSums(fit$raw_data[fit$item_variables])
  expect_equal(res[[1]]$df$obs, tabulate(sums[in_g1] + 1, nbins = 5))
})

# ---- Multidimensional ------------------------------------------------------

test_that("general result plus one result per dimension", {
  fit <- make_ppc_fit(itemD = c(1, 1, 2, 2, 2))
  res <- raw_score_ppc(fit)
  expect_length(res, 3)
  expect_equal(vapply(res, function(r) r$dim, numeric(1)), 0:2)
  expect_equal(vapply(res, function(r) r$nitem, numeric(1)), c(5, 2, 3))
  expect_equal(vapply(res, function(r) nrow(r$df), integer(1)), c(6, 3, 4))
})

test_that("observed scores use only the items of each dimension", {
  fit <- make_ppc_fit(itemD = c(1, 1, 2, 2, 2))
  res <- raw_score_ppc(fit)
  items <- fit$raw_data[fit$item_variables]
  expect_equal(res[[2]]$df$obs, tabulate(rowSums(items[, 1:2]) + 1, nbins = 3))
  expect_equal(res[[3]]$df$obs, tabulate(rowSums(items[, 3:5]) + 1, nbins = 4))
})

test_that("replicated scores use only the observations of each dimension", {
  fit <- make_ppc_fit(J = 40, itemD = c(1, 1, 2, 2, 2))
  res <- raw_score_ppc(fit)
  # one score per person and draw, with the maximum matching the dimension size
  expect_equal(unname(rowSums(res[[2]]$replicated)), rep(40, 5))
  expect_equal(unname(rowSums(res[[3]]$replicated)), rep(40, 5))
  expect_equal(ncol(res[[2]]$replicated), 3)
  expect_equal(ncol(res[[3]]$replicated), 4)
})

test_that("replicated dimension scores are the sums of that dimension's observations", {
  # deterministic y_rep: dimension 1 items always 1, dimension 2 items always 0
  fit <- make_ppc_fit(J = 10, draws = 3, itemD = c(1, 1, 2, 2))
  fit$model_info$y_rep[] <- rep(c(1, 1, 0, 0), times = 10)[col(fit$model_info$y_rep)]
  res <- raw_score_ppc(fit)
  expect_equal(unname(res[[2]]$replicated[, 3]), rep(10, 3))   # all persons score 2
  expect_equal(unname(res[[3]]$replicated[, 1]), rep(10, 3))   # all persons score 0
})

test_that("multidimensional with age_range: every dimension has all age groups", {
  res <- raw_score_ppc(make_ppc_fit(itemD = c(1, 1, 2, 2)), c(2, 6, 12))
  expect_length(res, 3)
  for (r in res) expect_length(r, 2)
  expect_equal(unique(vapply(res[[3]], function(g) g$dim, numeric(1))), 2)
})

# Fails with the current code: int_data[, data$itemD == d] drops to a vector
# for a one-item dimension and rowSums() errors. Use drop = FALSE.
test_that("a dimension with a single item works", {
  fit <- make_ppc_fit(itemD = c(1, 1, 1, 2))
  expect_no_error(res <- raw_score_ppc(fit))
  expect_equal(res[[3]]$nitem, 1)
  expect_equal(nrow(res[[3]]$df), 2)
})

# ---- Side effects ----------------------------------------------------------

test_that("the fit object is not modified", {
  fit <- make_ppc_fit(itemD = c(1, 1, 2, 2))
  before <- fit
  raw_score_ppc(fit, c(2, 6, 12))
  expect_identical(fit, before)
})