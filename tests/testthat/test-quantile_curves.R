

# ---- Helpers ---------------------------------------------------------------

# Unidimensional mock (irt_mod <= 4)
make_qfit_uni <- function(n = 60, seed = 1, irt_mod = 1) {
  set.seed(seed)
  age <- runif(n, 2, 12)
  info_sample <- data.frame(
    age = age,
    theta_mu = 0.3 * (age - 7),            # trait grows with age
    theta_sigma = runif(n, 0.8, 1.2),
    est_ability = rnorm(n, 0.3 * (age - 7)),
    irt_mod = irt_mod
  )
  structure(list(model_info = list(info_sample = info_sample)),
            class = "IRTnorm")
}

# Multidimensional mock (irt_mod > 4): general trait "theta" + domain traits
# Multidimensional mock (irt_mod > 4)
# traits: "theta", "normeddelta<d>", "delta<d>" for d = 1..n_dim
# with delta = normeddelta * theta_sigma + theta_mu
make_qfit_multi <- function(n = 60, seed = 2, irt_mod = 5, n_dim = 2) {
  set.seed(seed)
  age         <- runif(n, 2, 12)
  theta_mu    <- 0.3 * (age - 7)
  theta_sigma <- runif(n, 0.8, 1.2)
  
  make_rows <- function(trait, est_ability) {
    data.frame(
      age = age,
      theta_mu = theta_mu,
      theta_sigma = theta_sigma,
      est_ability = est_ability,
      trait = trait,
      irt_mod = irt_mod
    )
  }
  
  normed <- lapply(seq_len(n_dim), function(d) rnorm(n))
  
  rows <- c(
    list(make_rows("theta", rnorm(n, theta_mu, theta_sigma))),
    lapply(seq_len(n_dim), function(d)
      make_rows(paste0("normeddelta", d), normed[[d]])),
    lapply(seq_len(n_dim), function(d)
      make_rows(paste0("delta", d), normed[[d]] * theta_sigma + theta_mu))
  )
  
  structure(
    list(model_info = list(info_sample = do.call(rbind, rows)),
         newdata = list(J = n)),
    class = "IRTnorm"
  )
}

# Run the function on a null graphics device (it calls print(p))
run_qc <- function(...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  quantile_curves(...)
}

# ---- Input validation ------------------------------------------------------

test_that("errors when fit_info is not of class IRTnorm", {
  expect_error(run_qc(list(model_info = list())), "estimate_norming_model")
  expect_error(run_qc(NULL))
})

test_that("errors when probs are outside (0, 1)", {
  fit <- make_qfit_uni()
  expect_error(run_qc(fit, probs = c(0.1, 1)),   "probs")
  expect_error(run_qc(fit, probs = c(0, 0.5)),   "probs")
  expect_error(run_qc(fit, probs = c(-0.1, 0.5)), "probs")
  expect_error(run_qc(fit, probs = 1.5),         "probs")
})

test_that("probs validation happens before any data access", {
  # a bad IRTnorm object should still give the probs error first
  bad <- structure(list(), class = "IRTnorm")
  expect_error(run_qc(bad, probs = 2), "probs")
})

# ---- Unidimensional: return value ------------------------------------------

test_that("returns the quantile data invisibly", {
  fit <- make_qfit_uni()
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_invisible(quantile_curves(fit))
})

test_that("returned data has expected columns and class", {
  res <- run_qc(make_qfit_uni())
  expect_s3_class(res, "data.frame")
  expect_named(res, c("age", "quan", "qname"))
  expect_type(res$quan, "double")
  expect_type(res$qname, "integer")
})

test_that("returned data has one row per person and probability", {
  n <- 60
  fit <- make_qfit_uni(n = n)
  expect_equal(nrow(run_qc(fit)), n * 7)                       # default probs
  expect_equal(nrow(run_qc(fit, probs = c(0.25, 0.5, 0.75))), n * 3)
  expect_equal(nrow(run_qc(fit, probs = 0.5)), n)
})

test_that("qname labels follow the probs and are repeated for each person", {
  fit <- make_qfit_uni(n = 10)
  res <- run_qc(fit, probs = c(0.1, 0.5, 0.9))
  expect_setequal(unique(res$qname), c("p0.1", "p0.5", "p0.9"))
  expect_equal(as.character(res$qname), rep(c("p0.1", "p0.5", "p0.9"), 10))
})

test_that("default probs are the documented seven quantiles", {
  res <- run_qc(make_qfit_uni(n = 5))
  expect_setequal(as.character(unique(res$qname)),
                  paste0("p", c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)))
})

test_that("age is repeated once per probability, in person order", {
  fit <- make_qfit_uni(n = 10)
  probs <- c(0.25, 0.5, 0.75)
  res <- run_qc(fit, probs = probs)
  expect_equal(res$age, rep(fit$model_info$info_sample$age, each = length(probs)))
})

# ---- Unidimensional: quantile correctness ----------------------------------

test_that("quantiles equal qnorm(probs, theta_mu, theta_sigma) per person", {
  fit <- make_qfit_uni(n = 10)
  probs <- c(0.05, 0.5, 0.95)
  res <- run_qc(fit, probs = probs)
  s <- fit$model_info$info_sample
  
  expected <- as.vector(sapply(seq_len(nrow(s)), function(i) {
    qnorm(probs, s$theta_mu[i], s$theta_sigma[i])
  }))
  expect_equal(res$quan, expected)
})

test_that("the median quantile equals theta_mu", {
  fit <- make_qfit_uni(n = 20)
  res <- run_qc(fit, probs = 0.5)
  expect_equal(res$quan, fit$model_info$info_sample$theta_mu)
})

test_that("quantiles increase with probability within each person", {
  fit <- make_qfit_uni(n = 20)
  res <- run_qc(fit)
  by_person <- split(res$quan, rep(seq_len(20), each = 7))
  expect_true(all(vapply(by_person, function(q) !is.unsorted(q, strictly = TRUE),
                         logical(1))))
})

test_that("symmetric probabilities are symmetric around theta_mu", {
  fit <- make_qfit_uni(n = 15)
  lo <- run_qc(fit, probs = 0.1)$quan
  hi <- run_qc(fit, probs = 0.9)$quan
  mu <- fit$model_info$info_sample$theta_mu
  expect_equal(mu - lo, hi - mu)
})

test_that("larger theta_sigma widens the quantile spread", {
  fit <- make_qfit_uni(n = 2)
  fit$model_info$info_sample$theta_mu <- c(0, 0)
  fit$model_info$info_sample$theta_sigma <- c(1, 2)
  res <- run_qc(fit, probs = c(0.05, 0.95))
  spread <- c(diff(res$quan[1:2]), diff(res$quan[3:4]))
  expect_equal(spread[2], 2 * spread[1])
})

test_that("works for every unidimensional irt_mod value (1 to 4)", {
  for (m in 1:4) {
    expect_no_error(run_qc(make_qfit_uni(irt_mod = m)))
  }
})

# ---- Unidimensional: plot --------------------------------------------------

test_that("unidimensional plot has points and smooth layers, no facets", {
  run_qc(make_qfit_uni())
  p <- ggplot2::last_plot()
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_equal(as.character(geoms), c("GeomPoint", "GeomSmooth"))
  expect_s3_class(p$facet, "FacetNull")
  expect_equal(p$labels$y, "latent trait")
})

test_that("plot data come from the right sources", {
  fit <- make_qfit_uni(n = 30)
  res <- run_qc(fit)
  p <- ggplot2::last_plot()
  expect_equal(nrow(p$layers[[1]]$data), 30)
  expect_equal(p$layers[[2]]$data, res)
})

test_that("plot can be built without error", {
  run_qc(make_qfit_uni())
  expect_no_error(suppressWarnings(ggplot2::ggplot_build(ggplot2::last_plot())))
})

# ---- Multidimensional ------------------------------------------------------

test_that("multidimensional: only delta traits are used (not theta or normeddelta)", {
  res <- run_qc(make_qfit_multi())
  expect_setequal(unique(res$domain), c("delta1", "delta2"))
  expect_false(any(grepl("normeddelta|theta", res$domain)))
})

test_that("multidimensional: scales with the number of domains", {
  n <- 25
  res <- run_qc(make_qfit_multi(n = n, n_dim = 3))
  expect_equal(length(unique(res$domain)), 3)
  expect_equal(nrow(res), n * 7 * 3)
})

test_that("multidimensional: quantiles are computed per domain", {
  fit <- make_qfit_multi(n = 10)
  probs <- c(0.1, 0.5, 0.9)
  res <- run_qc(fit, probs = probs)
  s <- fit$model_info$info_sample
  d1 <- s[s$trait == "delta1", ]
  
  expected <- as.vector(sapply(seq_len(nrow(d1)), function(i) {
    qnorm(probs, d1$theta_mu[i], d1$theta_sigma[i])
  }))
  expect_equal(res$quan[res$domain == "delta1"], expected)
})

test_that("multidimensional: the median equals theta_mu for each domain", {
  fit <- make_qfit_multi(n = 20)
  res <- run_qc(fit, probs = 0.5)
  s <- fit$model_info$info_sample
  for (d in c("delta1", "delta2")) {
    expect_equal(res$quan[res$domain == d], s$theta_mu[s$trait == d])
  }
})

test_that("multidimensional plot is faceted by domain and uses est_ability", {
  fit <- make_qfit_multi(n = 30)
  run_qc(fit)
  p <- ggplot2::last_plot()
  expect_s3_class(p$facet, "FacetGrid")
  
  pts <- p$layers[[1]]$data
  expect_setequal(unique(pts$trait), c("delta1", "delta2"))
  expect_equal(nrow(pts), 30 * 2)

  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$layout$layout$PANEL)), 2)
})

test_that("multidimensional: returns data with a domain column", {
  res <- run_qc(make_qfit_multi())
  expect_s3_class(res, "data.frame")
  expect_named(res, c("age", "quan", "qname", "domain"))
})

test_that("multidimensional: one row per person x prob x domain", {
  n <- 40
  fit <- make_qfit_multi(n = n)
  expect_equal(nrow(run_qc(fit)), n * 7 * 2)
  expect_equal(nrow(run_qc(fit, probs = c(0.25, 0.75))), n * 2 * 2)
})

# ---- Edge cases ------------------------------------------------------------

test_that("function does not modify the input object", {
  fit_u <- make_qfit_uni(); before_u <- fit_u$model_info$info_sample
  run_qc(fit_u)
  expect_identical(fit_u$model_info$info_sample, before_u)
  
  fit_m <- make_qfit_multi(); before_m <- fit_m$model_info$info_sample
  run_qc(fit_m)
  expect_identical(fit_m$model_info$info_sample, before_m)
})

test_that("a single probability works", {
  res <- run_qc(make_qfit_uni(n = 30), probs = 0.5)
  expect_equal(as.character(unique(res$qname)), "p0.5")
})

test_that("unsorted probs keep the order given", {
  res <- run_qc(make_qfit_uni(n = 3), probs = c(0.9, 0.1))
  expect_equal(as.character(res$qname), rep(c("p0.9", "p0.1"), 3))
})
