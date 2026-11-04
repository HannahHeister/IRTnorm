test_that("matches a hand-built dummy matrix", {
  info <- matrix(c(1, 2, 3,     # item 1 (K = 3), persons 1-3
                   2, 1, 2),    # item 2 (K = 2)
                 nrow = 3)
  out <- make_dummy_matrix(info, K_item = c(3, 2))
  
  expected <- matrix(c(1, 0, 0,   # Item1_Cat1
                       0, 1, 0,   # Item1_Cat2
                       0, 0, 1,   # Item1_Cat3
                       0, 1, 0,   # Item2_Cat1
                       1, 0, 1),  # Item2_Cat2
                     nrow = 3,
                     dimnames = list(NULL, c("Item1_Cat1", "Item1_Cat2", "Item1_Cat3",
                                             "Item2_Cat1", "Item2_Cat2")))
  expect_equal(out, expected)
})

test_that("each person has exactly one 1 per item block", {
  set.seed(1)
  K    <- c(4, 2, 3)
  info <- sapply(K, function(k) sample(k, 20, replace = TRUE))
  out  <- make_dummy_matrix(info, K)
  blocks <- rep(seq_along(K), times = K)
  for (i in seq_along(K)) {
    expect_equal(unname(rowSums(out[, blocks == i, drop = FALSE])), rep(1, 20))
  }
})

test_that("a data.frame gives the same result as a matrix", {
  info <- matrix(c(1, 2, 2, 1), nrow = 2)
  expect_equal(make_dummy_matrix(as.data.frame(info), c(2, 2)),
               make_dummy_matrix(info, c(2, 2)))
})

test_that("a single item and a single person work", {
  out <- make_dummy_matrix(matrix(2, nrow = 1, ncol = 1), K_item = 3)
  expect_equal(dim(out), c(1, 3))
  expect_equal(as.vector(out), c(0, 1, 0))
})