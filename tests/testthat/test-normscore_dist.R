
# ---- Helpers ---------------------------------------------------------------

# Unidimensional mock (irt_mod <= 4)
make_fit_uni <- function(n = 300, seed = 1, irt_mod = 1,
                         age_min = 2, age_max = 12) {
  set.seed(seed)
  info_sample <- data.frame(
    age = runif(n, age_min, age_max),
    normscore_est = rnorm(n),
    irt_mod = irt_mod
  )
  structure(list(model_info = list(info_sample = info_sample)),
            class = "IRTnorm")
}

# Multidimensional mock (irt_mod > 4)
make_fit_multi <- function(n = 200, seed = 2, irt_mod = 5) {
  set.seed(seed)
  traits <- c("theta", "normeddelta1", "normeddelta2")
  info_sample <- do.call(rbind, lapply(traits, function(tr) {
    data.frame(
      age = runif(n, 2, 12),
      est_ability = rnorm(n),
      trait = tr,
      irt_mod = irt_mod
    )
  }))
  structure(list(model_info = list(info_sample = info_sample)),
            class = "IRTnorm")
}

# Wrapper to build the plot silently
get_plot <- function(...) {
  suppressWarnings(normscore_dist(...))
}

# ---- Input validation ------------------------------------------------------

test_that("errors when fit_info is not of class IRTnorm", {
  expect_error(
    normscore_dist(list(model_info = list()), age_range = c(2, 7, 12)),
    "estimate_norming_model"
  )
  expect_error(normscore_dist(NULL, age_range = c(2, 7, 12)))
})

test_that("errors when all norm scores are identical", {
  fit <- make_fit_uni()
  fit$model_info$info_sample$normscore_est <- 0.5
  expect_error(
    normscore_dist(fit, age_range = c(2, 7, 12)),
    "identical"
  )
})

test_that("multidimensional model requires dim_plot", {
  fit <- make_fit_multi()
  expect_error(
    normscore_dist(fit, age_range = c(2, 7, 12)),
    "dim_plot"
  )
})

test_that("identical est_ability in multidimensional model triggers error", {
  fit <- make_fit_multi()
  fit$model_info$info_sample$est_ability <- 1
  expect_error(
    normscore_dist(fit, age_range = c(2, 7, 12), dim_plot = "general"),
    "identical"
  )
})

# ---- Return value ----------------------------------------------------------

test_that("returns a ggplot object for unidimensional models", {
  fit <- make_fit_uni()
  p <- normscore_dist(fit, age_range = c(2, 7, 12))
  expect_s3_class(p, "ggplot")
})

test_that("works for every unidimensional irt_mod value (1 to 4)", {
  for (m in 1:4) {
    fit <- make_fit_uni(irt_mod = m)
    expect_s3_class(normscore_dist(fit, age_range = c(2, 7, 12)), "ggplot")
  }
})

test_that("returns a ggplot for multidimensional model with 'general'", {
  fit <- make_fit_multi()
  p <- normscore_dist(fit, age_range = c(2, 7, 12), dim_plot = "general")
  expect_s3_class(p, "ggplot")
})

test_that("returns a ggplot for multidimensional model with specific dimension", {
  fit <- make_fit_multi()
  # substr(dim_plot, 7, 7) extracts the dimension number, e.g. "domain1" -> "1"
  p1 <- normscore_dist(fit, age_range = c(2, 7, 12), dim_plot = "domain1")
  p2 <- normscore_dist(fit, age_range = c(2, 7, 12), dim_plot = "domain2")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

# ---- Plot structure --------------------------------------------------------

test_that("plot has histogram and line layers", {
  p <- normscore_dist(make_fit_uni(), age_range = c(2, 7, 12))
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomBar" %in% geoms)   # geom_histogram uses GeomBar
  expect_true("GeomLine" %in% geoms)
  expect_length(p$layers, 2)
})

test_that("plot has correct axis labels", {
  p <- normscore_dist(make_fit_uni(), age_range = c(2, 7, 12))
  expect_equal(p$labels$x, "norm score")
  expect_equal(p$labels$y, "frequency of individuals")
})

test_that("plot uses facets with one panel per age group", {
  fit <- make_fit_uni()
  p <- normscore_dist(fit, age_range = c(2, 5, 8, 12))
  expect_s3_class(p$facet, "FacetGrid")
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$layout$layout$PANEL)), 3)
})

test_that("facet order follows age_range", {
  fit <- make_fit_uni()
  p <- normscore_dist(fit, age_range = c(2, 5, 8, 12))
  built <- ggplot2::ggplot_build(p)
  labs <- as.character(built$layout$layout$`factor(...)`)
  # Fall back: just check the facet columns exist in the expected order
  expect_equal(
    labs,
    c("age: (2,5]", "age: (5,8]", "age: (8,12]")
  )
})

# ---- Binwidth and curve scaling -------------------------------------------

test_that("histogram uses the requested binwidth", {
  fit <- make_fit_uni()
  p <- normscore_dist(fit, age_range = c(2, 12), binwidth = 0.25)
  built <- ggplot2::ggplot_build(p)
  hist_data <- built$data[[1]]
  widths <- unique(round(hist_data$xmax - hist_data$xmin, 6))
  expect_equal(widths, 0.25)
})

test_that("normal curve data has 300 points per age group", {
  fit <- make_fit_uni()
  p <- normscore_dist(fit, age_range = c(2, 7, 12))
  curve_data <- p$layers[[2]]$data
  expect_equal(nrow(curve_data), 300 * 2)
  expect_equal(as.vector(table(curve_data$age_group)), c(300, 300))
})

test_that("curve x range is symmetric around zero and spans max |score|", {
  fit <- make_fit_uni()
  p <- normscore_dist(fit, age_range = c(2, 12))
  curve_data <- p$layers[[2]]$data
  extrem <- max(abs(fit$model_info$info_sample$normscore_est[
    fit$model_info$info_sample$age > 2 & fit$model_info$info_sample$age <= 12
  ]))
  expect_equal(min(curve_data$x), -extrem)
  expect_equal(max(curve_data$x), extrem)
})

test_that("curve height scales with binwidth and group size", {
  fit <- make_fit_uni()
  age_range <- c(2, 7, 12)
  bw <- 0.4
  p <- normscore_dist(fit, age_range = age_range, binwidth = bw)
  curve_data <- p$layers[[2]]$data
  
  s <- fit$model_info$info_sample
  n1 <- sum(s$age > 2 & s$age <= 7)
  group1 <- curve_data[curve_data$age_group == "age: (2,7]", ]
  expected <- dnorm(group1$x, 0, 1) * bw * n1
  expect_equal(group1$y, expected)
})

test_that("curve peak is at x = 0 region (standard normal shape)", {
  fit <- make_fit_uni()
  p <- normscore_dist(fit, age_range = c(2, 12))
  cd <- p$layers[[2]]$data
  expect_lt(abs(cd$x[which.max(cd$y)]), 0.05)
})

# ---- Age grouping ----------------------------------------------------------

test_that("age group intervals are left-open and right-closed", {
  fit <- make_fit_uni(n = 6)
  fit$model_info$info_sample$age <- c(2, 5, 5.0001, 7, 7.5, 12)
  fit$model_info$info_sample$normscore_est <- c(-1, -0.5, 0, 0.5, 1, 1.5)
  # age == 2 is excluded (lower bound open), so a warning is expected
  expect_warning(
    p <- normscore_dist(fit, age_range = c(2, 5, 12), binwidth = 0.5),
    "not within the plotted age groups"
  )
  hist_data <- p$layers[[1]]$data
  expect_equal(nrow(hist_data), 5)
  expect_equal(sum(hist_data$age_group == "age: (2,5]"), 1)   # age 5
  expect_equal(sum(hist_data$age_group == "age: (5,12]"), 4)
})

test_that("warns and removes individuals outside age_range", {
  fit <- make_fit_uni()
  expect_warning(
    p <- normscore_dist(fit, age_range = c(4, 6, 8)),
    "not within the plotted age groups"
  )
  hist_data <- p$layers[[1]]$data
  expect_false(any(is.na(hist_data$age_group)))
  expect_true(all(hist_data$age > 4 & hist_data$age <= 8))
})

test_that("no warning when all individuals fall within age_range", {
  fit <- make_fit_uni(age_min = 3, age_max = 11)
  expect_no_warning(normscore_dist(fit, age_range = c(2, 7, 12)))
})

test_that("age group labels match expected format", {
  fit <- make_fit_uni()
  p <- normscore_dist(fit, age_range = c(2, 7, 12))
  expect_setequal(
    unique(p$layers[[1]]$data$age_group),
    c("age: (2,7]", "age: (7,12]")
  )
})

test_that("a single age group (two-element age_range) works", {
  fit <- make_fit_uni()
  p <- normscore_dist(fit, age_range = c(2, 12))
  expect_s3_class(p, "ggplot")
  expect_equal(length(unique(p$layers[[2]]$data$age_group)), 1)
})

test_that("group counts sum to the number of retained individuals", {
  fit <- make_fit_uni()
  p <- normscore_dist(fit, age_range = c(2, 5, 8, 12))
  expect_equal(nrow(p$layers[[1]]$data), 300)
})

# ---- Multidimensional data selection ---------------------------------------

test_that("'general' uses only rows with trait == 'theta'", {
  fit <- make_fit_multi(n = 200)
  p <- normscore_dist(fit, age_range = c(2, 12), dim_plot = "general")
  d <- p$layers[[1]]$data
  expect_equal(nrow(d), 200)
  expect_true(all(d$trait == "theta"))
})

test_that("specific dimension selects the matching normeddelta trait", {
  fit <- make_fit_multi(n = 200)
  p <- normscore_dist(fit, age_range = c(2, 12), dim_plot = "domain2")
  d <- p$layers[[1]]$data
  expect_equal(nrow(d), 200)
  expect_true(all(d$trait == "normeddelta2"))
})

test_that("multidimensional plot uses est_ability as normscore_est", {
  fit <- make_fit_multi(n = 200)
  p <- normscore_dist(fit, age_range = c(2, 12), dim_plot = "general")
  d <- p$layers[[1]]$data
  expect_equal(d$normscore_est, d$est_ability)
})

test_that("dim_plot is ignored for unidimensional models", {
  fit <- make_fit_uni()
  expect_s3_class(
    normscore_dist(fit, age_range = c(2, 12), dim_plot = "general"),
    "ggplot"
  )
})

# ---- Edge cases ------------------------------------------------------------

test_that("function does not modify the input object", {
  fit <- make_fit_uni()
  before <- fit$model_info$info_sample
  invisible(normscore_dist(fit, age_range = c(2, 7, 12)))
  expect_identical(fit$model_info$info_sample, before)
})

test_that("plot can be built without error", {
  p <- normscore_dist(make_fit_uni(), age_range = c(2, 7, 12))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("works with a very small sample in a group", {
  fit <- make_fit_uni(n = 20)
  expect_s3_class(get_plot(fit, age_range = c(2, 3, 12)), "ggplot")
})

