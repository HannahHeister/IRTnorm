# tests/testthat/test-effect_age.R

# ---- Helpers ---------------------------------------------------------------

make_effect_uni <- function(n = 50, irt_mod = 1, seed = 1) {
  set.seed(seed)
  age <- runif(n, 2, 12)
  info_sample <- data.frame(
    age = age,
    theta_mu = 0.3 * (age - 7),
    theta_sigma = 1 + 0.02 * age,
    irt_mod = irt_mod
  )
  structure(list(model_info = list(info_sample = info_sample)),
            class = "IRTnorm")
}

# Multidimensional mock (irt_mod > 4). Rows are blocked by trait; persons share
# age / theta_mu / theta_sigma across traits. Traits: "theta", "normeddelta<d>",
# "delta<d>". nu: one value per dimension. lambda, unique_sigma: one per person.
make_effect_multi <- function(n = 30, irt_mod = 5, n_dim = 3, seed = 2,
                              nu = c(0.5, 0.8, 0.3)[seq_len(n_dim)],
                              lambda = NULL, unique_sigma = NULL) {
  set.seed(seed)
  age <- runif(n, 2, 12)
  mu <- 0.3 * (age - 7)
  sigma <- runif(n, 0.8, 1.2)
  if (is.null(lambda)) lambda <- runif(n, 0.3, 0.9)
  if (is.null(unique_sigma)) unique_sigma <- rep(runif(n, 0.4, 1),n_dim)
  rows <- function(trait) {
    data.frame(age = age, theta_mu = mu, theta_sigma = sigma,
               trait = trait, irt_mod = irt_mod)
  }
  info_sample <- do.call(rbind, c(
    list(rows("theta")),
    lapply(seq_len(n_dim), function(d) rows(paste0("normeddelta", d))),
    lapply(seq_len(n_dim), function(d) rows(paste0("delta", d)))
  ))
  structure(
    list(model_info = list(
      info_sample = info_sample,
      info_relation = list(nu = nu, lambda = lambda, unique_sigma = unique_sigma)
    )),
    class = "IRTnorm"
  )
}

# The function prints the plot, so run it on a null device
run_effect <- function(...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  effect_age(...)
}

layout_of <- function(p) ggplot2::ggplot_build(p)$layout$layout
# facet column lookup that does not depend on the exact column name
facet_col <- function(lay, name) lay[[grep(name, names(lay), value = TRUE)[1]]]

delta_rows <- function(fit) {
  s <- fit$model_info$info_sample
  s[substr(s$trait, 1, 1) == "d", ]
}

# ---- Input validation ------------------------------------------------------

test_that("errors when fit_info is not of class IRTnorm", {
  expect_error(run_effect(list()), "estimate_norming_model")
  expect_error(run_effect(NULL), "estimate_norming_model")
})

test_that("multidimensional fit without 'd' traits gives an informative error", {
  fit <- make_effect_multi()
  s <- fit$model_info$info_sample
  fit$model_info$info_sample <- s[s$trait %in% c("theta", "normeddelta1"), ]
  expect_error(run_effect(fit), "No traits starting with 'd'")
})

# ---- Return value ----------------------------------------------------------

test_that("returns plot and data invisibly", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_invisible(effect_age(make_effect_uni()))
  res <- effect_age(make_effect_uni())
  expect_named(res, c("plot", "data"))
  expect_s3_class(res$plot, "ggplot")
  expect_s3_class(res$data, "data.frame")
})

# ---- Unidimensional --------------------------------------------------------

test_that("uni: data are in long format with mean and sd stacked", {
  fit <- make_effect_uni(n = 50)
  d <- run_effect(fit)$data
  s <- fit$model_info$info_sample
  
  expect_named(d, c("age", "value", "moment"))
  expect_equal(nrow(d), 100)
  expect_equal(levels(d$moment), c("mean", "standard deviation"))
  expect_equal(d$value[d$moment == "mean"], s$theta_mu)
  expect_equal(d$value[d$moment == "standard deviation"], s$theta_sigma)
  expect_equal(d$age[d$moment == "mean"], s$age)
  expect_equal(d$age[d$moment == "standard deviation"], s$age)
})

test_that("uni: works for every irt_mod from 1 to 4", {
  for (m in 1:4) expect_no_error(run_effect(make_effect_uni(irt_mod = m)))
})

test_that("uni: one line layer, two free-scale panels, correct label", {
  p <- run_effect(make_effect_uni())$plot
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomLine")
  expect_equal(p$labels$y, "latent trait")
  expect_s3_class(p$facet, "FacetWrap")
  lay <- layout_of(p)
  expect_equal(nrow(lay), 2)
  expect_equal(as.character(facet_col(lay, "moment")), c("mean", "standard deviation"))
  expect_true(p$facet$params$free$x && p$facet$params$free$y)
  expect_no_error(ggplot2::ggplot_build(p))
})

# ---- Multidimensional: structure -------------------------------------------

test_that("multi: only delta traits are used (not theta or normeddelta)", {
  d <- run_effect(make_effect_multi())$data
  expect_setequal(unique(d$trait), c("delta1", "delta2", "delta3"))
})

test_that("multi: three moments for every person and trait", {
  n <- 30; D <- 2
  d <- run_effect(make_effect_multi(n = n, n_dim = D))$data
  expect_named(d, c("age", "value", "moment", "trait"))
  expect_equal(nrow(d), 3 * D * n)
  expect_equal(levels(d$moment), c("mean", "standard deviation", "correlation"))
  expect_equal(as.vector(table(d$moment, d$trait)), rep(n, 3 * D))
})

test_that("multi: scales with the number of dimensions", {
  d <- run_effect(make_effect_multi(n = 20, n_dim = 3))$data
  expect_equal(length(unique(d$trait)), 3)
  expect_equal(nrow(d), 3 * 3 * 20)
})

test_that("multi: mean and sd come from theta_mu and theta_sigma of the delta rows", {
  fit <- make_effect_multi()          # irt_mod 5
  d <- run_effect(fit)$data
  s <- delta_rows(fit)
  expect_equal(d$value[d$moment == "mean"], s$theta_mu)
  expect_equal(d$value[d$moment == "standard deviation"], s$theta_sigma)
  expect_equal(d$age[d$moment == "correlation"], s$age)
})

# ---- Multidimensional: correlation per model family ------------------------

test_that("irt_mod 5-8: correlation equals nu of each dimension", {
  for (m in 5:8) {
    fit <- make_effect_multi(n = 30, irt_mod = m, nu = c(0.5, 0.8), n_dim = 2)
    d <- run_effect(fit)$data
    expect_equal(d$value[d$moment == "correlation"], rep(c(0.5, 0.8), each = 30))
  }
})

test_that("irt_mod 9-10: correlation is nu (per dimension) times lambda (per person)", {
  n <- 30; D <- 2
  for (m in 9:10) {
    fit <- make_effect_multi(n = n, irt_mod = m, n_dim = D, nu = c(0.5, 0.8))
    d <- run_effect(fit)$data
    lambda <- fit$model_info$info_relation$lambda
    expect_equal(d$value[d$moment == "correlation"],
                 rep(c(0.5, 0.8), each = n) * rep(lambda, times = D))
  }
})

test_that("irt_mod 9-10: the same person has the same lambda in every dimension", {
  n <- 30
  fit <- make_effect_multi(n = n, irt_mod = 9, nu = c(0.5, 0.5), n_dim = 2)
  d <- run_effect(fit)$data
  cor <- d$value[d$moment == "correlation"]
  expect_equal(cor[1:n], cor[(n + 1):(2 * n)])
})

test_that("irt_mod 11-12: sd combines theta_sigma and unique_sigma (per person)", {
  n <- 30; D <- 3
  for (m in 11:12) {
    fit <- make_effect_multi(n = n, irt_mod = m, n_dim = D)
    d <- run_effect(fit)$data
    s <- delta_rows(fit)
    us <- fit$model_info$info_relation$unique_sigma
    expected_sd <- sqrt(s$theta_sigma^2 + us^2)
    expect_equal(d$value[d$moment == "standard deviation"], expected_sd)
    expect_equal(d$value[d$moment == "correlation"], s$theta_sigma / expected_sd)
  }
})

test_that("irt_mod 11-12: correlation lies in (0, 1) and sd exceeds theta_sigma", {
  fit <- make_effect_multi(irt_mod = 11)
  d <- run_effect(fit)$data
  cor <- d$value[d$moment == "correlation"]
  expect_true(all(cor > 0 & cor < 1))
  expect_true(all(d$value[d$moment == "standard deviation"] >
                    delta_rows(fit)$theta_sigma))
})

test_that("info_relation is read from model_info", {
  fit <- make_effect_multi(irt_mod = 5)
  fit$model_info$info_relation <- NULL
  expect_error(run_effect(fit))
})

# ---- Multidimensional: plot ------------------------------------------------

test_that("multi: one row per moment and column per trait with free y scales", {
  p <- run_effect(make_effect_multi(n = 30, n_dim = 2))$plot
  lay <- layout_of(p)
  expect_equal(nrow(lay), 6)
  expect_setequal(as.character(facet_col(lay, "moment")),
                  c("mean", "standard deviation", "correlation"))
  expect_setequal(as.character(facet_col(lay, "trait")), c("delta1", "delta2"))
  expect_true(p$facet$params$free$y)
  expect_false(p$facet$params$free$x)
  expect_equal(max(lay$ROW), 3)
})

test_that("multi: panel count scales with dimensions", {
  expect_equal(nrow(layout_of(run_effect(make_effect_multi(n_dim = 3))$plot)), 9)
})

test_that("multi: one line layer, builds without error for all model families", {
  for (m in c(5, 9, 11)) {
    p <- run_effect(make_effect_multi(irt_mod = m))$plot
    expect_equal(class(p$layers[[1]]$geom)[1], "GeomLine")
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

# ---- Side effects ----------------------------------------------------------

test_that("the input is not modified", {
  for (fit in list(make_effect_uni(), make_effect_multi(irt_mod = 9))) {
    before <- fit
    run_effect(fit)
    expect_identical(fit, before)
  }
})




