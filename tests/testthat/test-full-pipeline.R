# Simulate responses from a known model so the pipeline has something real to fit.
# 40 persons with continuous (unique) ages, so splines work and fits are stable.
sim_pipeline_data <- function(m, J = 40, I = 6, D = 2, K = 3, seed = 1, na_frac = 0) {
  set.seed(seed)
  multi <- m > 4
  itemD <- if (multi) rep(seq_len(D), each = I / D) else NULL
  age   <- runif(J, 1, 10)
  
  # general trait with a mild age trend; domain abilities load 0.7 on it
  g    <- 0.3 * as.numeric(scale(age)) + rnorm(J)
  abil <- if (multi) sapply(seq_len(D), function(d) 0.7 * g + sqrt(1 - 0.49) * rnorm(J))
  else matrix(g, J, 1)
  dom  <- if (multi) itemD else rep(1, I)
  
  has_alpha <- m %in% c(2, 3, 4, 6, 8, 10, 12)
  a <- if (has_alpha) runif(I, 0.8, 1.8) else rep(1, I)
  b <- seq(-1, 1, length.out = I)
  
  resp <- sapply(seq_len(I), function(i) {
    theta_i <- abil[, dom[i]]
    if (m == 4) {                                   # GRM, categories coded 1..K
      thr <- seq(-0.7, 0.7, length.out = K - 1) + b[i] / 2
      cum <- plogis(a[i] * outer(theta_i, thr, "-"))
      1 + rowSums(matrix(runif(length(cum)), nrow(cum)) < cum)
    } else {
      p <- plogis(a[i] * (theta_i - b[i]))
      if (m == 3) p <- 0.1 + 0.9 * p                # 3PL guessing
      rbinom(J, 1, p)
    }
  })
  colnames(resp) <- paste0("item", seq_len(I))
  if (na_frac > 0) resp[sample(length(resp), round(na_frac * length(resp)))] <- NA
  
  list(raw = data.frame(resp, age = age), itemD = itemD, theta = g, I = I, K = K,
       multi = multi, D = D)
}

# Arguments for estimate_norming_model(); extra args override
fit_args <- function(m, a, sim, ...) {
  args <- list(raw_data = sim$raw, age_variable = "age",
               item_variables = paste0("item", seq_len(sim$I)),
               irt_model = ALL_MODELS[m], age_model = a,
               iter_warmup = 100, iter_sampling = 100, chains = 2, seed = 1)
  if (a == "poly")    args <- c(args, list(poly_mean = 1, poly_sd = 1))
  if (a == "splines") args <- c(args, list(k_mean = 5, k_sd = 4, k_sd2 = 4))
  if (sim$multi)      args$itemD  <- sim$itemD
  if (m == 4)         args$K_item <- rep(sim$K, sim$I)
  modifyList(args, list(...))
}

skip_if_no_pipeline <- function() {
  skip_on_cran()
  skip_if_not(cmdstan_ok(), "cmdstan not available")
}

# One shared battery of checks for a fitted model
pipeline_check <- function(out, sim, n_draws = 200) {
  I <- sim$I; m <- out$newdata$irt_mod
  
  # structure
  expect_true(all(c("model_info", "norm_info", "conv_info", "newdata",
                    "model_stan_obj") %in% names(out)))
  
  # draw matrices: what norm_new_individual() indexes as x[m, items]
  ni <- out$norm_info
  expect_equal(nrow(ni$alpha), n_draws)
  expect_equal(ncol(ni$alpha), I)
  if (m == 4) {
    expect_equal(length(unlist(ni$threshold[1,])), I * (sim$K - 1))
    expect_equal(nrow(ni$threshold), n_draws)
  } else {
    expect_equal(dim(ni$gamma), c(n_draws, I))
  }
  expect_false(anyNA(ni$alpha))
  
  # downstream functions
  dom1 <- if (sim$multi) "domain1" else NULL
  dom2 <- if (sim$multi) "domain2" else NULL
  vis  <- function(x, ...) if (sim$multi){visual_raw_score_ppc(x, dim_plot = dom1)}else{ visual_raw_score_ppc(x)}
  
  ppc <- raw_score_ppc(fit_info = out)
  expect_no_error(vis(ppc))
  expect_no_error(vis(raw_score_ppc(fit_info = out, age_range = c(0, 5, 10))))
  expect_no_error(item_person_fit(fit_info = out, parameter = "item"))
  expect_no_error(item_person_fit(fit_info = out, parameter = "person"))
  expect_no_error(item_residual_cor(fit_info = out))
  expect_no_error(quantile_curves(fit_info = out))
  expect_no_error(effect_age(fit_info = out))
  if (sim$multi) {
    expect_no_error(normscore_dist(fit_info = out, age_range = seq(0, 10, 5),
                                   binwidth = 0.5, dim_plot = dom2))
  } else {
    expect_no_error(normscore_dist(fit_info = out, age_range = seq(0, 10, 5),
                                   binwidth = 0.5))
  }
  
  # norming a new person: finite output, also with a missing response
  y_new <- if (m == 4) sample(1:sim$K, I, replace = TRUE) else sample(0:1, I, replace = TRUE)
  res <- norm_new_individual(new_y = y_new, new_age = 5, fit_info = out)
  expect_true(all(is.finite(res$theta_samples)))
  y_na <- y_new; y_na[2] <- NA
  res_na <- norm_new_individual(new_y = y_na, new_age = 5, fit_info = out)
  expect_true(all(is.finite(res_na$theta_samples)))
  expect_equal(res_na$obs_items, setdiff(seq_len(I), 2))
  
  # sanity recovery: estimated ability is positively related to the true one
  # (weak on purpose, the chains are short)
  info <- out$model_info$info_sample
  est  <- if (sim$multi) info$est_ability[info$trait == "theta"] else info$theta_est
  expect_length(est, length(sim$theta))
  expect_gt(cor(est, sim$theta), 0)
}

# Specifications: how to pass the option, and how to check the result.
# `fixed`/`prior` use the nested formats data_prep() expects.
option_cases <- function(I = 6, D = 2, K = 3) {
  beta_fix  <- seq(-1, 1, length.out = I)
  b_fix     <- lapply(seq_len(I), function(i) c(-0.5, 0.5) + i / 20)
  tight     <- function(par, mu, n = I)
    setNames(list(rep(mu, n), rep(0.01, n)),
             paste0(par, c("_prior_mu", "_prior_sd")))
  
  list(
    list(name = "2PL fixed beta", m = 2,
         args = list(parameter_fixed = list(beta = beta_fix)),
         check = function(out) expect_equal(unname(colMeans(out$norm_info$beta)), beta_fix)),
    list(name = "2PL fixed alpha+beta", m = 2,
         args = list(parameter_fixed = list(alpha = rep(1.2, I), beta = beta_fix)),
         check = function(out) {
           expect_equal(unname(colMeans(out$norm_info$alpha)), rep(1.2, I))
           expect_equal(unname(colMeans(out$norm_info$beta)), beta_fix)
         }),
    list(name = "2PL prior alpha", m = 2,
         args = list(prior_knowledge = list(alpha = tight("alpha", 1.5))),
         check = function(out)
           expect_equal(unname(colMeans(out$norm_info$alpha)), rep(1.5, I), tolerance = 0.05)),
    list(name = "3PL fixed gamma", m = 3,
         args = list(parameter_fixed = list(gamma = rep(0.1, I))),
         check = function(out)
           expect_equal(unname(colMeans(out$norm_info$gamma)), rep(0.1, I))),
    list(name = "GRM fixed b", m = 4,
         args = list(parameter_fixed = list(b = b_fix)),
         check = function(out) {
           expect_equal(length(unlist(out$norm_info$threshold))/nrow(out$norm_info$threshold), I * (K - 1))
           expect_equal(length(unique(out$norm_info$threshold[, 1])), 1)  # constant over draws
         }),
    list(name = "GRM prior b", m = 4,
         args = list(prior_knowledge = list(b = list(
           b_prior_mu = lapply(seq_len(I), function(i) c(-0.5, 0.5)),
           b_prior_sd = lapply(seq_len(I), function(i) c(0.01, 0.01))))),
         check = function(out)
           expect_equal(mean(out$norm_info$threshold[, 1]), -0.5, tolerance = 0.05)),
    list(name = "HO1PL-C prior nu", m = 7,
         args = list(prior_knowledge = list(nu = tight("nu", 0.5, D))),
         check = function(out)
           expect_equal(unname(colMeans(out$norm_info$nu)), rep(0.5, D), tolerance = 0.05)),
    list(name = "2PL fixed beta + prior alpha", m = 2,
         args = list(parameter_fixed = list(beta = beta_fix),
                     prior_knowledge = list(alpha = tight("alpha", 1.5))),
         check = function(out) {
           expect_equal(unname(colMeans(out$norm_info$beta)), beta_fix)
           expect_equal(unname(colMeans(out$norm_info$alpha)), rep(1.5, I), tolerance = 0.05)
         })
  )
}


# Representative subset by default; IRTNORM_SLOW=all runs every model x age model
# Fixed subset of 13 fits (model position in ALL_MODELS, age model)
pipeline_grid <- function() {
  data.frame(
    m = c(1, 2, 3, 4, 6, 7, 9, 11,     # poly
          2, 3, 4, 7, 9),              # splines
    a = c(rep("poly", 8),
          rep("splines", 5)),
    stringsAsFactors = FALSE)
}

grid <- pipeline_grid()
for (r in seq_len(nrow(grid))) {
  local({
    m <- grid$m[r]; a <- grid$a[r]
    test_that(paste("real pipeline:", ALL_MODELS[m], a), {
      skip_if_no_pipeline()
      sim <- sim_pipeline_data(m)
      out <- suppressWarnings(do.call(estimate_norming_model, fit_args(m, a, sim)))
      pipeline_check(out, sim)
    })
  })
}

test_that("real pipeline: missing responses are handled in fitting", {
  skip_if_no_pipeline()
  sim <- sim_pipeline_data(2, na_frac = 0.1)
  out <- suppressWarnings(do.call(estimate_norming_model, fit_args(2, "poly", sim)))
  pipeline_check(out, sim)
})

test_that("real pi peline: fixed parameters and priors reach the fit", {
  for (case in option_cases()) {
    local({
      case <- case
      test_that(paste("real pipeline, option:", case$name), {
        skip_if_no_pipeline()
        sim <- sim_pipeline_data(case$m)
        out <- suppressWarnings(do.call(estimate_norming_model,
                       do.call(fit_args, c(list(case$m, "poly", sim), case$args))))
        case$check(out)
        pipeline_check(out, sim)     # draw shapes + downstream functions still work
      })
    })
  }
})

test_that("real pipeline: same seed gives identical draws", {
  skip_if_no_pipeline()
  sim <- sim_pipeline_data(1)
  o1 <- suppressWarnings(do.call(estimate_norming_model, fit_args(1, "poly", sim)))
  o2 <- suppressWarnings(do.call(estimate_norming_model, fit_args(1, "poly", sim)))
  expect_equal(o1$norm_info$beta, o2$norm_info$beta)
})

