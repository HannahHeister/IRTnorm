make_mock_ipf <- function(irt_mod = 3, J = 60, I = 6, D = 2,
                          K_item = NULL, seed = 1) {
  set.seed(seed)
  multi <- irt_mod > 4
  if (!multi) D <- 1
  itemD <- if (multi) rep(seq_len(D), length.out = I) else rep(1L, I)
  
  a <- runif(I, 0.8, 2); b <- rnorm(I, sd = 0.7); g <- runif(I, 0, 0.2)
  
  if (multi) {
    ability <- matrix(rnorm(J * D), ncol = D)
    info_sample <- data.frame(trait = rep(paste0("delta", seq_len(D)), each = J),
                              est_ability = as.vector(ability))
  } else {
    ability <- matrix(rnorm(J), ncol = 1)
    info_sample <- list(theta_est = ability[, 1])
  }
  theta_i <- ability[, itemD, drop = FALSE]        # J x I, ability of each item's domain
  
  if (irt_mod == 4) {
    if (is.null(K_item)) K_item <- rep(3L, I)
    K_max <- max(K_item)
    thr <- lapply(seq_len(I), function(i) sort(rnorm(K_item[i] - 1)))
    
    # Stan order: position i + I * (k - 1); unused slots stay NA on purpose
    threshold <- rep(NA_real_, I * (K_max - 1))
    for (i in seq_len(I)) for (k in seq_len(K_item[i] - 1))
      threshold[i + I * (k - 1)] <- thr[[i]][k]
    
    P <- lapply(seq_len(I), function(i) {
      cum <- plogis(a[i] * outer(theta_i[, i], thr[[i]], "-"))
      cbind(1, cum) - cbind(cum, 0)                 # J x K_i category probabilities
    })
    resp <- sapply(seq_len(I), function(i)
      apply(P[[i]], 1, function(p) sample.int(K_item[i], 1, prob = p)))  # coded 1..K
    info_item <- list(alpha = a, threshold = threshold)
    truth <- list(P = P)
  } else {
    K_item <- rep(2L, I); K_max <- 2L
    p <- sapply(seq_len(I), function(i)
      g[i] + (1 - g[i]) * plogis(a[i] * (theta_i[, i] - b[i])))
    resp <- matrix(rbinom(J * I, 1, p), nrow = J, ncol = I)
    info_item <- list(alpha = a, beta = b, gamma = g)
    truth <- list(p = p)
  }
  colnames(resp) <- paste0("item", seq_len(I))
  truth$resp <- resp
  
  structure(
    list(
      newdata = list(irt_mod = irt_mod, J = J, I = I, D = D, itemD = itemD,
                     K_item = K_item, K_max = K_max),
      raw_data = data.frame(resp, age = runif(J)),     # extra column on purpose
      item_variables = colnames(resp),
      model_info = list(info_sample = info_sample, info_item = info_item)
    ),
    class = "IRTnorm", truth = truth
  )
}
local_quiet_plot <- function(env = parent.frame()) {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off(), envir = env)
}
# Stand-in for the internal make_dummy_matrix(): one 0/1 column per category
# (codes 1..K), NA where the response is NA. Adapt to the real helper's coding.
mock_dummy <- function(info, K_item) {
  info <- as.matrix(info)
  do.call(cbind, lapply(seq_along(K_item), function(i)
    sapply(seq_len(K_item[i]), function(k) as.numeric(info[, i] == k))))
}

test_that("rejects objects that are not IRTnorm", {
  local_quiet_plot()
  expect_error(item_person_fit(list(), "item"), "estimate_norming_model")
  expect_error(item_person_fit(NULL,   "item"), "estimate_norming_model")
})

test_that("rejects a J that does not match the data", {
  local_quiet_plot()
  fit <- make_mock_ipf(3); fit$newdata$J <- 59
  expect_error(item_person_fit(fit, "item"), "data\\$J")
})

test_that("rejects an I that does not match the data", {
  local_quiet_plot()
  fit <- make_mock_ipf(3); fit$newdata$I <- 7
  expect_error(item_person_fit(fit, "item"), "data\\$I")
})

test_that("rejects itemD values outside 1..D", {
  local_quiet_plot()
  fit <- make_mock_ipf(9, D = 2); fit$newdata$itemD[1] <- 5
  expect_error(item_person_fit(fit, "item"), "itemD")
})

test_that("rejects an unknown parameter with a clear message", {
  local_quiet_plot()
  fit <- make_mock_ipf(3)
  expect_error(item_person_fit(fit, "foo"), "item.*person")
})

test_that("the input object is not modified", {
  local_quiet_plot()
  fit <- make_mock_ipf(3); before <- fit
  item_person_fit(fit, "person")
  expect_identical(fit, before)
})

test_that("returns a ggplot that can be built for all non-GRM models", {
  local_quiet_plot()
  skip_if_not_installed("ggplot2")
  for (m in c(1:3, 5:12)) for (par in c("item", "person")) {
    out <- item_person_fit(make_mock_ipf(m), par)
    expect_named(out, c("plot", "data"))
    expect_s3_class(out$plot, "ggplot")
    expect_no_error(ggplot2::ggplot_build(out$plot))
  }
})

test_that("titles identify the plot type", {
  local_quiet_plot()
  skip_if_not_installed("ggplot2", "3.5.0")
  fit <- make_mock_ipf(3)
  expect_match(ggplot2::get_labs(item_person_fit(fit, "item")$plot)$title,   "item fit")
  expect_match(ggplot2::get_labs(item_person_fit(fit, "person")$plot)$title, "person fit")
})

test_that("model 4 item fit builds as a plot", {
  local_quiet_plot()
  skip_if_not_installed("ggplot2")
  local_mocked_bindings(make_dummy_matrix = mock_dummy)
  out <- item_person_fit(make_mock_ipf(4, K_item = rep(3L, 6)), "item")
  expect_s3_class(out$plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(out$plot))
})

test_that("item fit matches an independent computation (models 1-3)", {
  local_quiet_plot()
  for (m in 1:3) {
    fit <- make_mock_ipf(m); tr <- attr(fit, "truth")
    d <- item_person_fit(fit, "item")$data
    expect_equal(nrow(d), 6L)
    expect_equal(d$observed, colMeans(tr$resp), ignore_attr = TRUE)
    expect_equal(d$irt_mod,  colMeans(tr$p),    ignore_attr = TRUE)
  }
})

test_that("multidimensional item fit uses each item's own domain ability", {
  local_quiet_plot()
  for (m in c(5, 9, 11)) {
    fit <- make_mock_ipf(m, D = 3, I = 6); tr <- attr(fit, "truth")
    d <- item_person_fit(fit, "item")$data
    expect_equal(d$observed, colMeans(tr$resp), ignore_attr = TRUE)
    expect_equal(d$irt_mod,  colMeans(tr$p),    ignore_attr = TRUE)
    expect_equal(d$domain, paste0("domain ", fit$newdata$itemD))
  }
})

test_that("missing responses are excluded from both observed and implied", {
  local_quiet_plot()
  fit <- make_mock_ipf(3); tr <- attr(fit, "truth")
  fit$raw_data[c(1, 5, 9), "item1"] <- NA
  fit$raw_data[c(2, 6),    "item3"] <- NA
  resp <- as.matrix(fit$raw_data[, fit$item_variables]); mask <- !is.na(resp)
  
  d <- item_person_fit(fit, "item")$data
  expect_equal(d$observed, colSums(resp, na.rm = TRUE) / colSums(mask), ignore_attr = TRUE)
  expect_equal(d$irt_mod,  colSums(tr$p * mask) / colSums(mask),        ignore_attr = TRUE)
})

test_that("items without any valid response are dropped with a warning", {
  local_quiet_plot()
  fit <- make_mock_ipf(3)
  fit$raw_data[, "item2"] <- NA
  expect_warning(p <- item_person_fit(fit, "item"), "no valid responses")
  expect_equal(nrow(p$data), 5L)
})

test_that("person fit equals raw score vs summed probabilities (models 1-3)", {
  local_quiet_plot()
  for (m in 1:3) {
    fit <- make_mock_ipf(m, J = 60, I = 6); tr <- attr(fit, "truth")
    d <- item_person_fit(fit, "person")$data
    expect_equal(nrow(d), 60L)
    expect_equal(d$observed, rowSums(tr$resp), ignore_attr = TRUE)
    expect_equal(d$irt_mod,  rowSums(tr$p),    ignore_attr = TRUE)
  }
})

test_that("person fit handles missing responses consistently", {
  local_quiet_plot()
  fit <- make_mock_ipf(3); tr <- attr(fit, "truth")
  fit$raw_data[c(1, 2), "item1"] <- NA
  resp <- as.matrix(fit$raw_data[, fit$item_variables]); mask <- !is.na(resp)
  d <- item_person_fit(fit, "person")$data
  expect_equal(d$observed, rowSums(resp, na.rm = TRUE), ignore_attr = TRUE)
  expect_equal(d$irt_mod,  rowSums(tr$p * mask),        ignore_attr = TRUE)
})

test_that("multidimensional person fit uses the items of each domain", {
  local_quiet_plot()
  fit <- make_mock_ipf(9, J = 60, I = 6, D = 2)
  tr <- attr(fit, "truth"); itemD <- fit$newdata$itemD
  d <- item_person_fit(fit, "person")$data
  
  exp_obs <- unlist(lapply(1:2, function(k) rowSums(tr$resp[, itemD == k, drop = FALSE])))
  exp_imp <- unlist(lapply(1:2, function(k) rowSums(tr$p[,    itemD == k, drop = FALSE])))
  expect_equal(nrow(d), 60L * 2)
  expect_equal(d$observed, exp_obs, ignore_attr = TRUE)
  expect_equal(d$irt_mod,  exp_imp, ignore_attr = TRUE)
  expect_equal(d$domain, rep(c("domain 1", "domain 2"), each = 60))
})


test_that("model 4 item fit: category proportions match", {
  local_quiet_plot()
  local_mocked_bindings(make_dummy_matrix = mock_dummy)
  fit <- make_mock_ipf(4, I = 6, K_item = rep(3L, 6)); tr <- attr(fit, "truth")
  K   <- fit$newdata$K_item
  d   <- item_person_fit(fit, "item")$data
  obs <- d$type == "observed"; imp <- d$type == "irt_mod"
  
  dummy <- mock_dummy(fit$raw_data[, fit$item_variables], K)
  expect_equal(nrow(d), 2 * sum(K))
  expect_equal(d$value[obs], colMeans(dummy),                ignore_attr = TRUE)
  expect_equal(d$value[imp], unlist(lapply(tr$P, colMeans)), ignore_attr = TRUE)
  expect_equal(d$item[obs],  rep(1:6, K))
  expect_equal(d$level[obs], sequence(K))
})

test_that("model 4 implied category probabilities sum to 1 per item", {
  local_quiet_plot()
  local_mocked_bindings(make_dummy_matrix = mock_dummy)
  fit <- make_mock_ipf(4, K_item = rep(3L, 6))
  d <- item_person_fit(fit, "item")$data
  d <- d[d$type == "irt_mod", ]
  expect_equal(as.vector(unname(tapply(d$value, d$item, sum))), rep(1, 6))
})

test_that("model 4 works with items that have different numbers of categories", {
  local_quiet_plot()
  local_mocked_bindings(make_dummy_matrix = mock_dummy)
  K <- c(2L, 3L, 4L, 3L, 2L, 3L)
  fit <- make_mock_ipf(4, I = 6, K_item = K); tr <- attr(fit, "truth")
  d <- item_person_fit(fit, "item")$data
  expect_equal(nrow(d), 2 * sum(K))
  expect_equal(d$value[d$type == "irt_mod"], unlist(lapply(tr$P, colMeans)),
               ignore_attr = TRUE)
})

test_that("model 4 person fit: raw score vs expected category score", {
  local_quiet_plot()
  local_mocked_bindings(make_dummy_matrix = mock_dummy)
  fit <- make_mock_ipf(4, K_item = rep(3L, 6)); tr <- attr(fit, "truth")
  d <- item_person_fit(fit, "person")$data
  
  exp_imp <- Reduce(`+`, lapply(tr$P, function(m) drop(m %*% seq_len(ncol(m)))))
  expect_equal(d$observed, rowSums(tr$resp), ignore_attr = TRUE)
  expect_equal(d$irt_mod,  exp_imp,          ignore_attr = TRUE)
})
