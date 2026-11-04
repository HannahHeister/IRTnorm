make_compare_info <- function() {
  best_model <- data.frame(
    rank      = 1:4,
    poly_mean = c(1, 2, NA, 3),
    poly_sd   = c(1, 1, NA, 2),
    ELPD_diff = c(0, -2, -5, -8),
    SE        = c(0, 1.5, 2, 2.5)   # SE == 0 for the best model on purpose
  )
  
  trajectory <- expand.grid(
    age         = 1:5,
    moment      = c("mean", "sd"),
    modelcombi  = c("a", "b", "c"),
    stringsAsFactors = FALSE
  )
  trajectory$rank <- match(trajectory$modelcombi, c("a", "b", "c"))
  trajectory$est  <- trajectory$age * trajectory$rank + (trajectory$moment == "sd")
  structure(
    list(best_model = best_model, trajectory = trajectory),
    class = "compare_info"
  )
}

test_that("both perspectives build without error or warning", {
  ci <- make_compare_info()
  
  for (persp in c("ELPD", "trajectory")) {
    p <- visual_age_models(ci, perspective = persp)
    expect_s3_class(p, "ggplot")
    expect_no_error(expect_no_warning(ggplot2::ggplot_build(p)))
  }
})

test_that("ELPD: SE == 0 is replaced and NA degrees become 'splines'", {
  ci <- make_compare_info()
  p  <- visual_age_models(ci, perspective = "ELPD")
  d  <- p$data
  
  # SE == 0 handling: no Inf/NaN in the plotted ratio
  expect_true(all(is.finite(d$ELPD_diff / d$SE)))
  expect_equal(d$SE[d$rank == 1], 1)
  
  # NA handling
  expect_false(anyNA(d$poly_mean))
  expect_false(anyNA(d$poly_sd))
  expect_true("splines" %in% d$poly_mean)
  expect_s3_class(d$poly_sd, "factor")
  expect_true("splines" %in% levels(d$poly_sd))
})

test_that("highlight_ranks controls what is highlighted, including edge cases", {
  ci <- make_compare_info()
  n_models <- nrow(ci$best_model)
  
  n_labels <- function(k) {
    p <- visual_age_models(ci, "ELPD", highlight_ranks = k)
    nrow(p$layers[[2]]$data)                 # the geom_text layer
  }
  n_highlighted <- function(k) {
    p <- visual_age_models(ci, "trajectory", highlight_ranks = k)
    length(unique(p$layers[[2]]$data$modelcombi))   # the highlight layer
  }
  
  # Typical values: both perspectives agree
  for (k in 1:3) {
    expect_equal(n_labels(k), k)
    expect_equal(n_highlighted(k), k)
  }
  
  # highlight_ranks = 0: nothing highlighted, but still builds
  expect_equal(n_labels(0), 0)
  expect_equal(n_highlighted(0), 0)
  for (persp in c("ELPD", "trajectory")) {
    expect_no_error(
      ggplot2::ggplot_build(visual_age_models(ci, persp, highlight_ranks = 0))
    )
  }
  
  # More highlights than models: everything highlighted, no error
  expect_equal(n_labels(10), n_models)
  expect_equal(n_highlighted(10), length(unique(ci$trajectory$modelcombi)))
  for (persp in c("ELPD", "trajectory")) {
    expect_no_error(
      ggplot2::ggplot_build(visual_age_models(ci, persp, highlight_ranks = 10))
    )
  }
})

