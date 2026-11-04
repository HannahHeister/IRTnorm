
# ---- Fixtures ----------------------------------------------------------------

make_fit_info <- function(irt_mod = 5, M = 3, I = 5, D = 2,
                          multidim = FALSE, K_max = 3) {
  data <- list(
    irt_mod = irt_mod, m_age = 5, s_age = 2,
    I = I, K_max = K_max, K_item = rep(if (irt_mod == 4) K_max else 2, I),
    itemD = if (multidim) rep(seq_len(D), length.out = I) else rep(1, I)
  )
  if (irt_mod > 4) data$D <- D   # irt_mod 4 deliberately has no D
  
  ni <- list(alpha = matrix(seq(0.5, 2, length.out = M * I), M, I))
  if (irt_mod == 4) {
    ni$beta <- matrix(seq(-1, 1, length.out = M * I * (K_max - 1)), M, I * (K_max - 1))
  } else {
    ni$beta  <- matrix(seq(-1, 1, length.out = M * I), M, I)
    ni$gamma <- matrix(seq(0.1, 0.3, length.out = M * I), M, I)
    ni$nu    <- if (multidim) matrix(0.4, M, D) else matrix(0.3, M, 1)
  }
  list(newdata = data, norm_info = ni)
}

make_age_effects <- function(M, D, multidim) {
  if (!multidim) {
    return(list(mu_post = matrix(0.2, M, 1), sigma_post = matrix(1.1, M, 1),
                lambda_post = NULL, unique_post = NULL))
  }
  list(
    mu_post     = array(rep(seq(0.1, by = -0.3, length.out = D), each = M), c(M, 1, D)),
    sigma_post  = array(rep(seq(1.1, 0.9, length.out = D), each = M), c(M, 1, D)),
    unique_post = array(0.4, c(M, 1, D)),
    lambda_post = matrix(0.8, M, D)
  )
}

# Install mocks for the duration of the calling test
use_mocks <- function(lp, multidim = FALSE, M = 3, D = 2, env = parent.frame()) {
  testthat::local_mocked_bindings(
    compute_age_effects = function(...) make_age_effects(M, D, multidim),
    log_post_draw_theta = lp,
    .package = "IRTnorm", .env = env
  )
}

# ---- Mock log-posteriors (signatures match the positional calls) -----------------

lp_uni_target <- function(target = 1, sd = 0.5)
  function(data, theta, a_vec, b_vec, c_vec, obs_K_item, y_vec, mu, sigma, ...)
    dnorm(theta, target, sd, log = TRUE)

lp_uni_flat <- function(data, theta, ...) 0

lp_uni_only_start <- function(data, theta, a_vec, b_vec, c_vec, obs_K_item,
                              y_vec, mu, sigma, ...)
  if (isTRUE(all.equal(theta, mu))) 0 else -Inf

lp_md <- function(m_theta = 0, s_theta = 1)
  function(data, theta, unique = NULL, ...)
    dnorm(theta, m_theta, s_theta, log = TRUE) + sum(dnorm(unique, 0, 1, log = TRUE))

make_recorder <- function() {
  rec <- new.env(); rec$first <- NULL
  rec$fn <- function(data, theta, a_vec, b_vec, c_vec, obs_K_item, y_vec, mu, sigma, ...) {
    if (is.null(rec$first))
      rec$first <- list(a_vec = a_vec, b_vec = b_vec, c_vec = c_vec,
                        obs_K_item = obs_K_item, y_vec = y_vec, mu = mu, sigma = sigma)
    0
  }
  rec
}

fast_mh <- list(n_iter = 200, burnin = 100, thin = 1, proposal_sd = 0.6)
long_mh <- list(n_iter = 5000, burnin = 500, thin = 1, proposal_sd = 1)

# ==============================================================================
# UNIDIMENSIONAL
# ==============================================================================

test_that("unidimensional: output structure and lengths", {
  fi <- make_fit_info(M = 3); use_mocks(lp_uni_target(), M = 3)
  res <- norm_new_individual(c(1, 0, 1, 1, 0), 7, fi, fast_mh)
  
  expect_named(res, c("theta_samples", "theta_zsample", "by_draw_samples",
                      "acceptance_rates", "summary", "mh_args", "obs_items"),
               ignore.order = TRUE)
  n_keep <- length(seq(101, 200, by = 1))
  expect_length(res$theta_samples, 3 * n_keep)
  expect_length(res$theta_zsample, 3 * n_keep)
  expect_length(res$by_draw_samples, 3)
  expect_true(all(lengths(res$by_draw_samples) == n_keep))
  expect_length(res$acceptance_rates, 3)
  expect_true(all(res$acceptance_rates >= 0 & res$acceptance_rates <= 1))
  expect_named(res$summary, c("mean", "sd", "median", "ci"))
  expect_length(res$summary$ci, 2)
  expect_equal(res$mh_args, fast_mh)
  expect_false(anyNA(res$theta_samples))
})

test_that("unidimensional: summary matches the z-sample", {
  fi <- make_fit_info(); use_mocks(lp_uni_target())
  set.seed(1)
  res <- norm_new_individual(c(1, 0, 1, 1, 0), 7, fi, fast_mh)
  expect_equal(res$summary$mean, mean(res$theta_zsample))
  expect_equal(res$summary$sd, sd(res$theta_zsample))
  expect_equal(res$summary$median, median(res$theta_zsample))
  expect_equal(unname(res$summary$ci),
               unname(quantile(res$theta_zsample, c(0.025, 0.975))))
  expect_lt(res$summary$ci[1], res$summary$ci[2])
})

test_that("thin and burnin control which samples are kept", {
  fi <- make_fit_info(M = 2); use_mocks(lp_uni_target(), M = 2)
  mh <- list(n_iter = 200, burnin = 100, thin = 5, proposal_sd = 0.6)
  res <- norm_new_individual(c(1, 0, 1, 1, 0), 7, fi, mh)
  expect_true(all(lengths(res$by_draw_samples) == length(seq(101, 200, by = 5))))
  
  mh2 <- list(n_iter = 50, burnin = 0, thin = 1, proposal_sd = 0.6)
  res2 <- norm_new_individual(c(1, 0, 1, 1, 0), 7, fi, mh2)
  expect_true(all(lengths(res2$by_draw_samples) == 50))
})

test_that("results are reproducible with the same seed and differ across seeds", {
  fi <- make_fit_info(); use_mocks(lp_uni_target())
  y <- c(1, 0, 1, 1, 0)
  set.seed(42); r1 <- norm_new_individual(y, 7, fi, fast_mh)
  set.seed(42); r2 <- norm_new_individual(y, 7, fi, fast_mh)
  set.seed(43); r3 <- norm_new_individual(y, 7, fi, fast_mh)
  expect_identical(r1$theta_samples, r2$theta_samples)
  expect_false(identical(r1$theta_samples, r3$theta_samples))
})

test_that("NA responses are dropped and parameters subset consistently", {
  fi <- make_fit_info(M = 2, I = 5); rec <- make_recorder()
  use_mocks(rec$fn, M = 2)
  res <- norm_new_individual(c(1, NA, 0, NA, 1), 7, fi, fast_mh)
  obs <- c(1, 3, 5)
  
  expect_equal(res$obs_items, obs)
  expect_equal(rec$first$y_vec, c(1, 0, 1))
  expect_equal(rec$first$a_vec, fi$norm_info$alpha[1, obs])
  expect_equal(rec$first$b_vec, fi$norm_info$beta[1, obs])
  expect_equal(rec$first$c_vec, fi$norm_info$gamma[1, obs])
  expect_equal(rec$first$obs_K_item, fi$newdata$K_item[obs])
})

test_that("starting values are the population mean/sd for the first draw", {
  fi <- make_fit_info(); rec <- make_recorder(); use_mocks(rec$fn)
  norm_new_individual(c(1, 0, 1, 1, 0), 7, fi, fast_mh)
  expect_equal(rec$first$mu, 0.2)
  expect_equal(rec$first$sigma, 1.1)
})

test_that("each posterior draw uses its own parameter row", {
  # Target of the mock = first observed item's alpha, which differs by draw.
  fi <- make_fit_info(M = 3)
  fi$norm_info$alpha[, 1] <- c(1, 2, 3)
  use_mocks(function(data, theta, a_vec, ...) dnorm(theta, a_vec[1], 0.5, log = TRUE), M = 3)
  set.seed(2)
  res <- norm_new_individual(c(1, 0, 1, 1, 0), 7, fi, long_mh)
  means <- vapply(res$by_draw_samples, mean, numeric(1))
  expect_equal(means, c(1, 2, 3), tolerance = 0.15)
})

test_that("sampler recovers a known target and z-transform is correct", {
  fi <- make_fit_info(M = 2); use_mocks(lp_uni_target(1, 0.5), M = 2)
  set.seed(3)
  res <- norm_new_individual(c(1, 0, 1, 1, 0), 7, fi, long_mh)
  expect_equal(mean(res$theta_samples), 1, tolerance = 0.1)
  expect_equal(sd(res$theta_samples), 0.5, tolerance = 0.1)
  # mu = 0.2, sigma = 1.1 from the mocked age effects
  expect_equal(res$theta_zsample, (res$theta_samples - 0.2) / 1.1)
  expect_equal(res$summary$mean, (1 - 0.2) / 1.1, tolerance = 0.1)
})

test_that("acceptance rate behaves sensibly with proposal size", {
  fi <- make_fit_info(M = 2)
  y <- c(1, 0, 1, 1, 0)
  
  use_mocks(lp_uni_flat, M = 2)
  set.seed(4)
  expect_true(all(norm_new_individual(y, 7, fi, fast_mh)$acceptance_rates == 1))
  
  use_mocks(lp_uni_target(1, 0.5), M = 2)
  set.seed(4)
  small <- norm_new_individual(y, 7, fi, modifyList(fast_mh, list(proposal_sd = 0.001)))
  huge  <- norm_new_individual(y, 7, fi, modifyList(fast_mh, list(proposal_sd = 100)))
  expect_true(all(small$acceptance_rates > 0.9))
  expect_true(all(huge$acceptance_rates < 0.2))
})

test_that("non-finite log posterior proposals are never accepted", {
  fi <- make_fit_info(M = 2); use_mocks(lp_uni_only_start, M = 2)
  res <- norm_new_individual(c(1, 0, 1, 1, 0), 7, fi, fast_mh)
  expect_true(all(res$acceptance_rates == 0))
  expect_true(all(res$theta_samples == 0.2))   # stuck at the mean
  expect_true(all(res$theta_zsample == 0))
})

test_that("burnin >= n_iter errors rather than silently returning garbage", {
  fi <- make_fit_info(M = 1); use_mocks(lp_uni_target(), M = 1)
  expect_error(norm_new_individual(c(1, 0, 1, 1, 0), 7, fi,
                                   list(n_iter = 100, burnin = 100, thin = 1, proposal_sd = 0.6)))
})

test_that("all-NA responses run and report no observed items", {
  # With the real log_post_draw_theta this should sample from the prior.
  fi <- make_fit_info(M = 2); use_mocks(lp_uni_target(), M = 2)
  res <- norm_new_individual(rep(NA_real_, 5), 7, fi, fast_mh)
  expect_length(res$obs_items, 0)
  expect_false(anyNA(res$theta_samples))
})

test_that("GRM (irt_mod 4) indexes threshold columns per observed item", {
  fi <- make_fit_info(irt_mod = 4, M = 2, I = 4, K_max = 3)  # no D, no gamma/nu
  fi$norm_info$threshold <- fi$norm_info$beta
  rec <- make_recorder(); use_mocks(rec$fn, M = 2)
  res <- norm_new_individual(c(1, NA, 2, NA), 7, fi, fast_mh)
  expect_equal(res$obs_items, c(1, 3))
  # beta has I*(K_max-1) columns; items 1 and 3 -> cols 1,5 and 3,7 (in loop order)
  expect_equal(rec$first$b_vec, fi$norm_info$beta[1, c(1, 5, 3, 7)])
  expect_equal(rec$first$a_vec, fi$norm_info$alpha[1, c(1, 3)])
  expect_null(rec$first$c_vec)
})

# ==============================================================================
# MULTIDIMENSIONAL
# ==============================================================================

test_that("multidimensional: output structure and dimensions", {
  M <- 2; D <- 2
  fi <- make_fit_info(irt_mod = 9, M = M, I = 6, D = D, multidim = TRUE)
  use_mocks(lp_md(), multidim = TRUE, M = M, D = D)
  res <- norm_new_individual(c(1, 0, 1, 1, 0, 1), 7, fi, fast_mh)
  
  n_keep <- 100
  expect_named(res, c("theta_samples", "zdelta_sample", "delta_sample",
                      "by_draw_samples", "summary", "mh_args", "obs_items"),
               ignore.order = TRUE)
  expect_length(res$theta_samples, M * n_keep)
  expect_equal(dim(res$zdelta_sample), c(M * n_keep, D))
  expect_equal(dim(res$delta_sample),  c(M * n_keep, D))
  expect_length(res$by_draw_samples, M)
  expect_equal(dim(res$by_draw_samples[[1]]), c(n_keep, 2 * D + 1))
  
  expect_s3_class(res$summary, "data.frame")
  expect_equal(nrow(res$summary), D + 1)
  expect_equal(res$summary$domain, c("general", "domain1", "domain2"))
  expect_named(res$summary, c("mean", "sd", "median", "lower2_5", "upper97_5", "domain"))
  expect_false(anyNA(res$summary[, 1:5]))
  expect_true(all(res$summary$lower2_5 <= res$summary$median))
  expect_true(all(res$summary$median <= res$summary$upper97_5))
})

test_that("multidimensional: summary matches the z-delta samples", {
  fi <- make_fit_info(irt_mod = 9, M = 2, I = 6, multidim = TRUE)
  use_mocks(lp_md(), multidim = TRUE, M = 2)
  set.seed(5)
  res <- norm_new_individual(c(1, 0, 1, 1, 0, 1), 7, fi, fast_mh)
  expect_equal(res$summary$mean[-1], unname(colMeans(res$zdelta_sample)))
  expect_equal(res$summary$mean[1], mean(res$theta_samples))
})

test_that("multidimensional: z-delta equals (delta - mu) / sigma (irt_mod < 9, 9, 10)", {
  mu <- c(0.1, -0.2); sigma <- c(1.1, 0.9)
  for (m in c(5, 9, 10)) {
    fi <- make_fit_info(irt_mod = m, M = 2, I = 6, multidim = TRUE)
    use_mocks(lp_md(), multidim = TRUE, M = 2)
    set.seed(6)
    res <- norm_new_individual(c(1, 0, 1, 1, 0, 1), 7, fi, fast_mh)
    s <- res$by_draw_samples[[1]]
    expected_z <- t((t(s[, 4:5]) - mu) / sigma)
    expect_equal(unname(s[, 2:3]), unname(expected_z), info = paste("irt_mod", m))
  }
})

test_that("multidimensional: loadings give unit-variance domains and correct correlation", {
  # theta ~ N(0,1), unique ~ N(0,1)  =>  z_delta ~ N(0,1), cor(theta, z) = nu * lambda
  M <- 2; D <- 2
  cases <- list(`9` = 0.4 * 0.8, `10` = 0.4 * 0.8, `5` = 0.4)  # model 5: nu only
  for (m in names(cases)) {
    fi <- make_fit_info(irt_mod = as.integer(m), M = M, I = 6, D = D, multidim = TRUE)
    use_mocks(lp_md(0, 1), multidim = TRUE, M = M, D = D)
    set.seed(7)
    res <- norm_new_individual(c(1, 0, 1, 1, 0, 1), 7, fi, long_mh)
    expect_equal(unname(colMeans(res$zdelta_sample)), c(0, 0), tolerance = 0.12, info = m)
    expect_equal(unname(apply(res$zdelta_sample, 2, sd)), c(1, 1), tolerance = 0.12, info = m)
    expect_equal(cor(res$theta_samples, res$zdelta_sample[, 1]), cases[[m]],
                 tolerance = 0.12, info = m)
  }
})

test_that("multidimensional: recovers the general-trait target", {
  fi <- make_fit_info(irt_mod = 9, M = 2, I = 6, multidim = TRUE)
  use_mocks(lp_md(0.5, 0.5), multidim = TRUE, M = 2)
  set.seed(8)
  res <- norm_new_individual(c(1, 0, 1, 1, 0, 1), 7, fi, long_mh)
  expect_equal(mean(res$theta_samples), 0.5, tolerance = 0.12)
  expect_equal(sd(res$theta_samples), 0.5, tolerance = 0.12)
})

test_that("multidimensional: irt_mod 11/12 z-delta uses the population unique sd", {
  # NOTE: the source divides by sqrt(sigma^2 + unique_curr^2), i.e. the *sampled*
  # unique effect. This test encodes the expected (unique_post) version and may
  # fail -- if so, check whether that line is a typo.
  mu <- c(0.1, -0.2); sigma <- c(1.1, 0.9); uniq <- 0.4
  for (m in c(11, 12)) {
    fi <- make_fit_info(irt_mod = m, M = 2, I = 6, multidim = TRUE)
    use_mocks(lp_md(), multidim = TRUE, M = 2)
    set.seed(9)
    res <- norm_new_individual(c(1, 0, 1, 1, 0, 1), 7, fi, fast_mh)
    s <- res$by_draw_samples[[1]]
    expected_z <- t((t(s[, 4:5]) - mu) / sqrt(sigma^2 + uniq^2))
    expect_equal(unname(s[, 2:3]), unname(expected_z), info = paste("irt_mod", m))
  }
})

test_that("multidimensional NA handling subsets items and itemD", {
  fi <- make_fit_info(irt_mod = 9, M = 2, I = 6, multidim = TRUE)
  use_mocks(lp_md(), multidim = TRUE, M = 2)
  res <- norm_new_individual(c(1, NA, 1, NA, 0, 1), 7, fi, fast_mh)
  expect_equal(res$obs_items, c(1, 3, 5, 6))
})