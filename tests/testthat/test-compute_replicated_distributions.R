# tests/testthat/test-compute_replicated_distributions.R

# 4 persons x 2 draws, nitem = 2
#   draw 1 scores: 0, 1, 1, 2  -> counts for scores 0:2 = 1, 2, 1
#   draw 2 scores: 2, 2, 0, 1  -> counts for scores 0:2 = 1, 1, 2
score_rep <- matrix(c(0, 1, 1, 2,
                      2, 2, 0, 1), nrow = 4, ncol = 2)

test_that("returns draws x (nitem + 1) matrix of score counts", {
  res <- compute_replicated_distributions(score_rep, nitem = 2, draws = 2)
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(2, 3))
  expect_equal(res, rbind(c(1, 2, 1),
                          c(1, 1, 2)))
})

test_that("each row sums to the number of persons", {
  # also catches scores that fall outside 0:nitem and are silently dropped
  res <- compute_replicated_distributions(score_rep, nitem = 2, draws = 2)
  expect_equal(rowSums(res), rep(nrow(score_rep), 2))
})

test_that("age_filter restricts persons (logical and numeric index agree)", {
  expected <- rbind(c(1, 1, 0),
                    c(0, 0, 2))
  res_lgl <- compute_replicated_distributions(
    score_rep, age_filter = c(TRUE, TRUE, FALSE, FALSE), nitem = 2, draws = 2)
  res_idx <- compute_replicated_distributions(
    score_rep, age_filter = 1:2, nitem = 2, draws = 2)
  expect_equal(res_lgl, expected)
  expect_equal(res_idx, expected)
})

test_that("age_filter selecting no one gives all zeros", {
  res <- compute_replicated_distributions(
    score_rep, age_filter = rep(FALSE, 4), nitem = 2, draws = 2)
  expect_equal(res, matrix(0, nrow = 2, ncol = 3))
})

test_that("age_filter selecting a single person works", {
  res <- compute_replicated_distributions(
    score_rep, age_filter = 4, nitem = 2, draws = 2)
  expect_equal(res, rbind(c(0, 0, 1),    # draw 1: score 2
                          c(0, 1, 0)))   # draw 2: score 1
})

test_that("only the first `draws` columns are used", {
  res <- compute_replicated_distributions(score_rep, nitem = 2, draws = 1)
  expect_equal(res, matrix(c(1, 2, 1), nrow = 1))
})