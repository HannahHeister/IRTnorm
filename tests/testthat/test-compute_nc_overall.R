score_rep <- matrix(c(0, 1, 1, 2,
                      2, 2, 0, 1), nrow = 4, ncol = 2)
score_obs <- c(0, 1, 1, 2)

test_that("adds age label and dim to the process_group result", {
  res <- compute_nc_overall(score_obs, score_rep, nitem = 2, draws = 2)
  expect_equal(res$age, "full age range")
  expect_equal(res$dim, 0)                       # default
  expect_equal(compute_nc_overall(score_obs, score_rep, 2, 2, dim = 3)$dim, 3)
})
