test_that("empty input returns zero counts for all possible scores", {
  expect_equal(
    create_full_distribution(numeric(0), nitem = 3),
    data.frame(score = 0:3, obs = 0)
  )
})

test_that("missing scores are filled with zeros", {
  res <- create_full_distribution(c(1, 1, 3, 5, 5, 5), nitem = 5)
  expect_equal(res$score, 0:5)
  expect_equal(res$obs, c(0, 2, 0, 1, 0, 3))
})

test_that("missing scores at the edges (0 and nitem) are filled", {
  res <- create_full_distribution(c(2, 3, 3), nitem = 5)
  expect_equal(res$obs, c(0, 0, 1, 2, 0, 0))
})

test_that("complete distribution is returned unchanged in structure", {
  res <- create_full_distribution(c(0, 1, 1, 2, 2, 2), nitem = 2)
  expect_named(res, c("score", "obs"))
  expect_equal(res$score, 0:2)
  expect_equal(res$obs, c(1, 2, 3))
  expect_type(res$score, "integer")  
})