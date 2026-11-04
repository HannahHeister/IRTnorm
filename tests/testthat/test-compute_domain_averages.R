# tests/testthat/test-compute_domain_averages.R

# 4 items, domains 1,1,2,2. Only the upper triangle is read.
#   within pairs:  (1,2) = 0.1, (3,4) = 0.3
#   between pairs: (1,3) = 0.2, (1,4) = 0.4, (2,3) = 0.6, (2,4) = 0.8
make_q3 <- function() {
  q <- matrix(0, 4, 4)
  q[1, 2] <- 0.1; q[3, 4] <- 0.3
  q[1, 3] <- 0.2; q[1, 4] <- 0.4; q[2, 3] <- 0.6; q[2, 4] <- 0.8
  q
}
itemD <- c(1, 1, 2, 2)

test_that("returns within and between averages with pair counts", {
  res <- compute_domain_averages(make_q3(), itemD)
  expect_equal(res$domains, c("within", "between"))
  expect_equal(res$mean_correlation, c(mean(c(0.1, 0.3)), mean(c(0.2, 0.4, 0.6, 0.8))))
  expect_equal(res$n_pairs, c(2, 4))
})

test_that("only the upper triangle is used", {
  q <- make_q3()
  q[lower.tri(q)] <- 99          # garbage below the diagonal must be ignored
  diag(q) <- 99                  # and on the diagonal
  expect_equal(compute_domain_averages(q, itemD)$mean_correlation,
               c(0.2, 0.5))
})

test_that("NA pairs are skipped in both means and counts", {
  q <- make_q3()
  q[1, 2] <- NA                  # within
  q[2, 4] <- NA                  # between
  res <- compute_domain_averages(q, itemD)
  expect_equal(res$mean_correlation, c(0.3, mean(c(0.2, 0.4, 0.6))))
  expect_equal(res$n_pairs, c(1, 3))
})

test_that("domain labels need not be consecutive or numeric order", {
  res <- compute_domain_averages(make_q3(), c(5, 2, 2, 5))
  # within: (2,3) = 0.6, (1,4) = 0.4; between: the other four pairs
  expect_equal(res$mean_correlation[1], 0.5)
  expect_equal(res$n_pairs, c(2, 4))
})


test_that("pair counts add up to J * (J - 1) / 2", {
  set.seed(1)
  q <- matrix(rnorm(25), 5, 5)
  res <- compute_domain_averages(q, c(1, 1, 2, 2, 3))
  expect_equal(sum(res$n_pairs), 10)
})