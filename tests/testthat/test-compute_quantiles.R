# tests/testthat/test-compute_quantiles.R

# 5 draws x 3 scores. Columns are chosen so the quantiles are easy to verify:
#   score 0: 1..5       -> Q5 = 1.2, Q50 = 3, Q95 = 4.8
#   score 1: 10,10,...  -> constant, all quantiles 10
#   score 2: 5..1       -> same values as score 0, in reversed order
replicated <- cbind(1:5,
                    rep(10, 5),
                    5:1)

test_that("returns one row per score and columns Q5, Q50, Q95", {
  res <- compute_quantiles(replicated)
  expect_s3_class(res, "data.frame")
  expect_equal(dim(res), c(3, 3))
  expect_named(res, c("Q5", "Q50", "Q95"))
})

test_that("quantiles are computed per column (score), not per row (draw)", {
  res <- compute_quantiles(replicated)
  expect_equal(res$Q5,  c(1.2, 10, 1.2))
  expect_equal(res$Q50, c(3,   10, 3))
  expect_equal(res$Q95, c(4.8, 10, 4.8))
})

test_that("row order follows column order of the input", {
  # an asymmetric matrix catches an accidental transposition or reordering
  m <- cbind(rep(1, 5), rep(2, 5), rep(3, 5))
  expect_equal(compute_quantiles(m)$Q50, c(1, 2, 3))
})

test_that("Q5 <= Q50 <= Q95 for every score", {
  set.seed(1)
  m <- matrix(rpois(200 * 6, lambda = 4), nrow = 200, ncol = 6)
  res <- compute_quantiles(m)
  expect_true(all(res$Q5 <= res$Q50))
  expect_true(all(res$Q50 <= res$Q95))
})

test_that("a single score (one column) still returns a one-row data frame", {
  res <- compute_quantiles(matrix(1:5, ncol = 1))
  expect_equal(dim(res), c(1, 3))
  expect_equal(res$Q50, 3)
})

test_that("a single draw (one row) gives identical quantiles", {
  res <- compute_quantiles(matrix(c(2, 4, 6), nrow = 1))
  expect_equal(res$Q5,  c(2, 4, 6))
  expect_equal(res$Q50, c(2, 4, 6))
  expect_equal(res$Q95, c(2, 4, 6))
})