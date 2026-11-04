test_that("matches a hand-computed case", {
  p <- grm_probs(0, 1, c(-1, 1))
  expect_equal(as.vector(p), c(0.2689414, 0.4621172, 0.2689414), tolerance = 1e-6)
})

test_that("rows are probability distributions for ordered thresholds", {
  p <- grm_probs(seq(-4, 4, length.out = 25), a = 2, b = c(-1.5, -0.3, 0.8, 2))
  expect_equal(rowSums(p), rep(1, 25))
  expect_true(all(p >= 0))
})

test_that("a single threshold reduces to the 2PL", {
  theta <- c(-1, 0, 2)
  expect_equal(grm_probs(theta, 1.4, 0.5)[, 2], plogis(1.4 * (theta - 0.5)))
})