# tests/testthat/test-compute_nc_by_age.R

# 6 persons x 2 draws, nitem = 2
# ages: 3, 5, 5.5, 8, 10, 12
raw_data <- data.frame(id = 1:6, age = c(3, 5, 5.5, 8, 10, 12))
observed_scores <- c(0, 1, 1, 2, 2, 0)
score_rep <- cbind(c(0, 1, 1, 2, 2, 0),
                   c(1, 1, 2, 2, 0, 0))

run_by_age <- function(age_range, ...) {
  compute_nc_by_age(observed_scores, raw_data, "age", age_range,
                    score_rep, nitem = 2, draws = 2, ...)
}

test_that("returns one result per age group with the expected structure", {
  res <- run_by_age(c(2, 6, 12))
  expect_type(res, "list")
  expect_length(res, 2)
  for (r in res) {
    expect_true(all(c("df", "replicated", "nitem", "age", "dim") %in% names(r)))
  }
})

test_that("age labels follow the [lower,upper] format in order", {
  res <- run_by_age(c(2, 6, 12))
  expect_equal(vapply(res, function(r) r$age, character(1)),
               c("age: [2,6]", "age: [6,12]"))
})

test_that("intervals are closed on the left and right", {
  res <- run_by_age(c(2, 5.5, 12))
  expect_equal(sum(res[[1]]$df$obs), 3)   # ages 3, 5, 5.5
  expect_equal(sum(res[[2]]$df$obs), 4)   # ages 5.5, 8, 10, 12
})

test_that("observed scores are filtered to the persons of each group", {
  res <- run_by_age(c(2, 6, 12))
  # group 1: persons 1-3, scores 0, 1, 1 -> counts 1, 2, 0
  expect_equal(res[[1]]$df$obs, c(1, 2, 0))
  # group 2: persons 4-6, scores 2, 2, 0 -> counts 1, 0, 2
  expect_equal(res[[2]]$df$obs, c(1, 0, 2))
})

test_that("replicated data use the same persons as the observed data", {
  res <- run_by_age(c(2, 6, 12))
  # group 1, draw 1: persons 1-3 have scores 0, 1, 1 -> counts 1, 2, 0
  expect_equal(res[[1]]$replicated[1, ], c(1, 2, 0))
  # group 2, draw 2: persons 4-6 have scores 2, 0, 0 -> counts 2, 0, 1
  expect_equal(res[[2]]$replicated[2, ], c(2, 0, 1))
})

test_that("dim is passed to every group (default 0)", {
  expect_true(all(vapply(run_by_age(c(2, 6, 12)),
                         function(r) r$dim == 0, logical(1))))
  expect_true(all(vapply(run_by_age(c(2, 6, 12), dim = 2),
                         function(r) r$dim == 2, logical(1))))
})

test_that("a single age group (two-element age_range) works", {
  res <- run_by_age(c(2, 12))
  expect_length(res, 1)
  expect_equal(sum(res[[1]]$df$obs), 6)
})

test_that("persons outside age_range are excluded", {
  res <- run_by_age(c(4, 9))              # only ages 5, 5.5, 8
  expect_equal(sum(res[[1]]$df$obs), 3)
  expect_equal(nrow(res[[1]]$replicated), 2)
})

test_that("age_variable selects the correct column", {
  rd <- data.frame(age_months = c(3, 5, 5.5, 8, 10, 12), age = 99)
  res <- compute_nc_by_age(observed_scores, rd, "age_months", c(2, 6, 12),
                           score_rep, nitem = 2, draws = 2)
  expect_equal(sum(res[[1]]$df$obs), 3)
})

test_that("an empty age group gives zero observed counts", {
  res <- run_by_age(c(2, 6, 7, 12))       # no one in (6, 7]
  expect_equal(res[[2]]$df$obs, c(0, 0, 0))
  expect_equal(res[[2]]$replicated, matrix(0, nrow = 2, ncol = 3))
})