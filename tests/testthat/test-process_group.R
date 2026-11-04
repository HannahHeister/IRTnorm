# tests/testthat/test-process_group.R

# 4 persons x 2 draws, nitem = 2
#   draw 1 counts for scores 0:2 = 1, 2, 1
#   draw 2 counts for scores 0:2 = 1, 1, 2
score_rep <- matrix(c(0, 1, 1, 2,
                      2, 2, 0, 1), nrow = 4, ncol = 2)
score_obs <- c(0, 1, 1, 2)

test_that("returns list with df, replicated and nitem", {
  res <- process_group(score_obs, score_rep, nitem = 2, draws = 2)
  expect_named(res, c("df", "replicated", "nitem"))
  expect_equal(res$nitem, 2)
})

test_that("df combines observed counts and replicated quantiles per score", {
  res <- process_group(score_obs, score_rep, nitem = 2, draws = 2)
  expect_named(res$df, c("score", "obs", "Q5", "Q50", "Q95"))
  expect_equal(res$df$score, 0:2)
  expect_equal(res$df$obs, c(1, 2, 1))
  # quantiles of the two draws per score (type 7 with n = 2)
  expect_equal(res$df$Q5,  c(1, 1.05, 1.05))
  expect_equal(res$df$Q50, c(1, 1.5,  1.5))
  expect_equal(res$df$Q95, c(1, 1.95, 1.95))
})

test_that("replicated matrix is returned unchanged from the helper", {
  res <- process_group(score_obs, score_rep, nitem = 2, draws = 2)
  expect_equal(res$replicated, rbind(c(1, 2, 1),
                                     c(1, 1, 2)))
})

test_that("age_filter is applied to the replicated data", {
  # Guards the positional call compute_replicated_distributions(score_rep,
  # age_filter, nitem, draws): swapped arguments would fail here.
  res <- process_group(score_obs = c(0, 1), score_rep = score_rep, nitem = 2,
                       draws = 2, age_filter = c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(res$replicated, rbind(c(1, 1, 0),
                                     c(0, 0, 2)))
  expect_equal(res$df$obs, c(1, 1, 0))
  expect_equal(nrow(res$df), 3)
})

test_that("empty observed scores give zero counts but keep replicated quantiles", {
  res <- process_group(numeric(0), score_rep, nitem = 2, draws = 2)
  expect_equal(res$df$obs, c(0, 0, 0))
  expect_equal(res$df$Q50, c(1, 1.5, 1.5))
})