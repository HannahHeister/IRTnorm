
# ---- Helpers ---------------------------------------------------------------

# One result as produced by compute_nc_overall() / compute_nc_by_age().
# Observed counts of dimension d are shifted by 1000 * d, so tests can tell
# which dimension ended up in the plot.
make_result <- function(nitem = 3, draws = 6, age = "full age range",
                        dim = 0, seed = 1) {
  set.seed(seed)
  rep_mat <- matrix(rpois(draws * (nitem + 1), 5), nrow = draws, ncol = nitem + 1)
  q <- function(p) as.numeric(apply(rep_mat, 2, quantile, probs = p))
  list(
    df = data.frame(score = 0:nitem,
                    obs = rpois(nitem + 1, 5) + 1000 * dim,
                    Q5 = q(0.05), Q50 = q(0.5), Q95 = q(0.95)),
    replicated = rep_mat,
    nitem = nitem,
    age = age,
    dim = dim
  )
}

# raw_score_ppc() without age_range: list(general, dim1, dim2, ...)
# nitem has one value per entry (general first)
make_info_overall <- function(nitem = 3) {
  structure(lapply(seq_along(nitem) - 1L, function(d)
    make_result(nitem[d + 1], dim = d, seed = d + 1)), class = "IRTnorm_ppc")
}

# raw_score_ppc() with age_range: list(general = list(groups), dim1 = list(groups), ...)
age_labels <- c("age: [2,6]", "age: [6,12]")
make_info_age <- function(nitem = 3, ages = age_labels) {
  structure(lapply(seq_along(nitem) - 1L, function(d)
    lapply(seq_along(ages), function(i)
      make_result(nitem[d + 1], age = ages[i], dim = d, seed = 10 * d + i))), 
    class = "IRTnorm_ppc")
}

# The function prints the plot, so run it on a null device
run_vis <- function(...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  suppressWarnings(visual_raw_score_ppc(...))
}

layer_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))
panels <- function(p) ggplot2::ggplot_build(p)$layout$layout

# ---- Return value and fixed plot structure ---------------------------------

test_that("returns a ggplot invisibly", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(visual_raw_score_ppc(make_info_overall()))
  expect_s3_class(visual_raw_score_ppc(make_info_overall()), "ggplot")
})

test_that("plot has jittered replications, three quantile lines and observed points", {
  p <- run_vis(make_info_overall())
  expect_equal(unname(layer_geoms(p)),
               c("GeomPoint", "GeomLine", "GeomLine", "GeomLine", "GeomPoint"))
  linetypes <- vapply(p$layers[2:4], function(l) l$aes_params$linetype, character(1))
  expect_equal(unname(linetypes), c("dashed", "dotted", "dotted"))
})

test_that("axis labels are set", {
  p <- run_vis(make_info_overall())
  expect_equal(p$labels$x, "raw score")
  expect_equal(p$labels$y, "frequency of individuals")
})

test_that("plot can be built without error for every input shape", {
  for (info in list(make_info_overall(),
                    make_info_overall(c(5, 2, 3)),
                    make_info_age(),
                    make_info_age(c(5, 2, 3)))) {
    p <- run_vis(info, dim_plot = if (length(info) > 1 && is.list(info[[1]][[1]]) &&
                                      !is.null(info[[1]][[1]]$age)) "general")
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

# ---- Shape 1: single result (unidimensional, no age_range) ------------------

test_that("single result: data are reshaped to long format correctly", {
  res <- structure(list(make_result(nitem = 3, draws = 6)), class = "IRTnorm_ppc")
  p <- run_vis(res)
  d <- p$layers[[1]]$data
  
  expect_equal(nrow(d), 6 * 4)                       # draws x (nitem + 1)
  expect_equal(sort(unique(d$score)), 0:3)
  expect_type(d$score, "double")                     # converted from column names
  expect_equal(sort(unique(d$id)), 1:6)
  
  # orientation: row = draw, column = score + 1
  for (i in c(1, 4, 6)){
    for (s in c(0, 2, 3)) {
      expect_equal(d$frequency[d$id == i & d$score == s], res[[1]]$replicated[i, s + 1])
    }
  } 
})

test_that("single result: quantile lines and observed points use the right columns", {
  res <- structure(list(make_result()), class = "IRTnorm_ppc")
  built <- ggplot2::ggplot_build(run_vis(res))
  expect_equal(built$data[[2]]$y, res[[1]]$df$Q50)
  expect_equal(built$data[[3]]$y, res[[1]]$df$Q5)
  expect_equal(built$data[[4]]$y, res[[1]]$df$Q95)
  expect_equal(built$data[[5]]$y, res[[1]]$df$obs)
  expect_equal(built$data[[5]]$x, res[[1]]$df$score)
})

test_that("single result: no facets", {
  p <- run_vis(structure(list(make_result()), class = "IRTnorm_ppc"))
  expect_s3_class(p$facet, "FacetNull")
})

# ---- Shape 2: age groups, unidimensional ------------------------------------

test_that("age groups: one facet panel per group with matching labels", {
  p <- run_vis(make_info_age())
  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(nrow(panels(p)), 2)
  expect_setequal(as.character(panels(p)$age_group), age_labels)
})

test_that("age groups: data of all groups are stacked", {
  p <- run_vis(make_info_age(nitem = 3))                # 2 groups, 6 draws
  expect_equal(nrow(p$layers[[1]]$data), 2 * 6 * 4)
  expect_equal(nrow(p$layers[[5]]$data), 2 * 4)
  expect_setequal(unique(p$layers[[5]]$data$age_group), age_labels)
})

test_that("each age group keeps its own observed counts", {
  info <- make_info_age()[[1]]
  built <- ggplot2::ggplot_build(run_vis(structure(list(info),class = "IRTnorm_ppc")))
  expect_equal(built$data[[5]]$y, c(info[[1]]$df$obs, info[[2]]$df$obs))
})

test_that("wrapped and unwrapped age-group input give the same plot data", {
  as_ppc <- function(x) structure(x, class = "IRTnorm_ppc")
  groups <- make_info_age()[[1]]            # list of age-group results (unwrapped)
  
  p_wrapped   <- run_vis(as_ppc(list(groups)))
  p_unwrapped <- run_vis(as_ppc(groups))
  
  expect_equal(p_wrapped$layers[[1]]$data, p_unwrapped$layers[[1]]$data)
  expect_equal(p_wrapped$layers[[5]]$data, p_unwrapped$layers[[5]]$data)
})

test_that("a single age group is plotted without facets", {
  p <- run_vis(make_info_age(ages = "age: [2,12]"))
  expect_s3_class(p$facet, "FacetNull")
  expect_equal(nrow(p$layers[[5]]$data), 4)
})

# ---- Shape 3: multidimensional, no age_range --------------------------------

test_that("multidimensional: facets by dimension with readable labels", {
  p <- run_vis(make_info_overall(c(5, 2, 3)))
  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(nrow(panels(p)), 3)
  expect_equal(levels(p$layers[[5]]$data[["domain"]]),
               c("general", "domain1", "domain2"))
  expect_equal(levels(p$layers[[1]]$data[["domain"]]),
               c("general", "domain1", "domain2"))
})

test_that("multidimensional: dimensions with different numbers of items are stacked", {
  p <- run_vis(make_info_overall(c(5, 2, 3)))
  expect_equal(nrow(p$layers[[1]]$data), 6 * (6 + 3 + 4))       # draws x (nitem + 1)
  expect_equal(nrow(p$layers[[5]]$data), 6 + 3 + 4)
  d <- p$layers[[5]]$data
  expect_equal(max(d$score[d[["domain"]] == "domain1"]), 2)
})

test_that("multidimensional: each domain holds its own observed counts", {
  info <- make_info_overall(c(3, 3, 3))
  d <- run_vis(info)$layers[[5]]$data
  for (k in 0:2) {
    expect_equal(d$obs[d[["domain"]] == c("general", "domain1", "domain2")[k + 1]],
                 info[[k + 1]]$df$obs)
  }
})

# ---- Shape 4: multidimensional with age groups ------------------------------

test_that("multidimensional with age groups requires dim_plot", {
  expect_error(run_vis(make_info_age(c(3, 3, 3))), "dim_plot")
})

test_that("dim_plot = 'general' plots the general results by age group", {
  p <- run_vis(make_info_age(c(3, 3, 3)), dim_plot = "general")
  expect_equal(nrow(panels(p)), 2)
  expect_true(all(p$layers[[5]]$data$obs < 1000))               # dimension 0 only
  expect_equal(p$labels$title, "general")
})

# These two describe the intended behaviour and FAIL with the current indexing
# (list position 1 is the general result, so domain d is at position d + 1).
test_that("dim_plot = 'domain1' plots dimension 1", {
  p <- run_vis(make_info_age(c(3, 3, 3)), dim_plot = "domain1")
  obs <- p$layers[[5]]$data$obs
  expect_true(all(obs >= 1000 & obs < 2000))
  expect_equal(nrow(panels(p)), 2)
  expect_equal(p$labels$title, "domain1")
})

test_that("dim_plot = 'domain2' plots dimension 2", {
  p <- run_vis(make_info_age(c(3, 3, 3)), dim_plot = "domain2")
  obs <- p$layers[[5]]$data$obs
  expect_true(all(obs >= 2000 & obs < 3000))
})

test_that("malformed dim_plot gives a clear error", {
  info <- make_info_age(c(3, 3, 3))
  expect_error(suppressWarnings(run_vis(info, dim_plot = "domainX")), "dim_plot has to be one of the following values")
  expect_error(run_vis(info, dim_plot = "domain"), "dim_plot has to be one of the following values")
  expect_error(run_vis(info, dim_plot = "domain9"), "dim_plot has to be one of the following values")
  })



# ---- dim_plot and title ----------------------------------------------------

test_that("no title when dim_plot is NULL", {
  expect_null(run_vis(make_info_overall())$labels$title)
  expect_null(run_vis(make_info_age())$labels$title)
})

test_that("dim_plot only sets the title when dimension selection is needed", {
  p <- run_vis(make_info_overall(c(3, 3, 3)), dim_plot = "general")
  expect_equal(p$labels$title,NULL)
  
  p <- run_vis(make_info_age(c(3, 3, 3)), dim_plot = "general")
  expect_equal(p$labels$title,"general")
  })

# ---- Side effects ----------------------------------------------------------

test_that("the input is not modified", {
  for (info in list(make_info_overall(c(3, 3, 3)), make_info_age(c(3, 3, 3)))) {
    before <- info
    run_vis(info, dim_plot = "general")
    expect_identical(info, before)
  }
})