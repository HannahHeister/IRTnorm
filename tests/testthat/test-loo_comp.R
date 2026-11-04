# Helper: build the list(elpd_summary, pointwise_matrix) that loo_comp() expects.
# ELPD is derived from the pointwise matrix so the two are always consistent.
# NOTE: the summary needs >= 2 columns, otherwise `loo_info[[1]][best_list, ]`
# drops to a vector and loo_comp() errors (true of real input, which has 6 columns).
make_loo_info <- function(pw) {
  df <- data.frame(
    ELPD      = colSums(pw),
    irt_model = colnames(pw),
    row.names = colnames(pw)
  )
  list(elpd_value = df, loo_ELPD = pw)
}

# 4 observations, 3 models, deliberately NOT in ranked order.
# ELPD: m1 = -6, m2 = -4 (best), m3 = -5
pw_basic <- cbind(
  m1 = c(-2, -1, -2, -1),
  m2 = c(-1, -1, -1, -1),
  m3 = c(-1, -2, -1, -1)
)

test_that("models are ordered by decreasing ELPD", {
  out <- loo_comp(make_loo_info(pw_basic))
  expect_equal(rownames(out), c("m2", "m3", "m1"))
  expect_false(is.unsorted(rev(out$ELPD)))
})

test_that("ELPD_diff is relative to the best model", {
  out <- loo_comp(make_loo_info(pw_basic))
  expect_equal(out$ELPD_diff, c(0, -1, -2))
})

test_that("best model has SE of 0 and rank 1", {
  out <- loo_comp(make_loo_info(pw_basic))
  expect_equal(out$SE[1], 0)
  expect_equal(out$rank, 1:3)
})

test_that("SE of the difference matches sqrt(n * var(d)) by hand", {
  out <- loo_comp(make_loo_info(pw_basic))
  # m2 - m3 = c(0,1,0,0): var = 0.25, n = 4 -> SE = 1
  expect_equal(out$SE[out$irt_model == "m3"], 1)
  # m2 - m1 = c(1,0,1,0): var = 1/3, n = 4 -> SE = sqrt(4/3)
  expect_equal(out$SE[out$irt_model == "m1"], sqrt(4 / 3))
})

test_that("result does not depend on input column order", {
  a <- loo_comp(make_loo_info(pw_basic))
  b <- loo_comp(make_loo_info(pw_basic[, c("m3", "m1", "m2")]))
  expect_equal(a, b)
})

test_that("works when the best model is already first", {
  out <- loo_comp(make_loo_info(pw_basic[, c("m2", "m3", "m1")]))
  expect_equal(rownames(out), c("m2", "m3", "m1"))
})

test_that("two models are handled", {
  out <- loo_comp(make_loo_info(pw_basic[, c("m1", "m2")]))
  expect_equal(rownames(out), c("m2", "m1"))
  expect_equal(out$rank, 1:2)
})

test_that("original columns of the summary are preserved", {
  out <- loo_comp(make_loo_info(pw_basic))
  expect_true(all(c("ELPD", "irt_model", "ELPD_diff", "SE", "rank") %in% names(out)))
})

test_that("ties in best ELPD are handled", {
  pw_tie <- cbind(
    m1 = c(-1, -1, -1, -1),
    m2 = c(-2, 0, -2, 0),
    m3 = c(-2, -2, -1, -1)
  )
  expect_equal(colSums(pw_tie)[c("m1", "m2")], c(m1 = -4, m2 = -4))
  out <- NULL
  expect_no_error(out <- loo_comp(make_loo_info(pw_tie)))
  expect_equal(out$rank, c(1,1,2))
  expect_equal(out$ELPD_diff[1:2], c(0, 0))
  expect_equal(out$SE[1], 0)
})