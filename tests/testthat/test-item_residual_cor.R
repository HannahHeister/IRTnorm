const_draws <- function(true, n = 20) matrix(rep(true, each = n), nrow = n)

make_mock_q3 <- function(irt_mod = 3, J = 80, I = 6, D = 2, seed = 1,
                         with_gamma = TRUE) {
  set.seed(seed)
  itemD <- if (irt_mod > 4) rep(seq_len(D), length.out = I) else NULL
  
  a <- runif(I, 0.8, 2); b <- rnorm(I); g <- runif(I, 0, 0.2)
  
  if (irt_mod > 4) {
    abil <- rnorm(J * D)
    info_sample <- data.frame(trait = rep(paste0("delta", seq_len(D)), each = J),
                              est_ability = abil)
    delta <- matrix(abil, ncol = D)          # J x D, column d = domain d
    theta_i <- sapply(seq_len(I), function(i) delta[, itemD[i]])
  } else {
    theta <- rnorm(J)
    info_sample <- list(theta_est = theta)
    theta_i <- sapply(seq_len(I), function(i) theta)
  }
  
  p <- sapply(seq_len(I), function(i)
    g[i] + (1 - g[i]) * plogis(a[i] * (theta_i[, i] - b[i])))
  resp <- matrix(rbinom(J * I, 1, p), nrow = J, ncol = I)
  colnames(resp) <- paste0("item", seq_len(I))
  
  norm_info <- list(alpha = const_draws(a), beta = const_draws(b))
  if (with_gamma) norm_info$gamma <- const_draws(g)
  
  structure(
    list(
      newdata = list(irt_mod = irt_mod, J = J, I = I, D = D, itemD = itemD,
                     K_item = rep(2L, I), K_max = 2L),
      raw_data = data.frame(resp, age = runif(J)),   # extra column on purpose
      item_variables = colnames(resp),
      norm_info = norm_info,
      model_info = list(info_sample = info_sample)
    ),
    class = "IRTnorm",
    truth = list(a = a, b = b, g = g, theta_i = theta_i, resp = resp)
  )
}

# independent expectation: Q3 = cor(observed - expected)
expected_q3 <- function(fit) {
  t <- attr(fit, "truth")
  E <- sapply(seq_along(t$a), function(i)
    t$g[i] + (1 - t$g[i]) * plogis(t$a[i] * (t$theta_i[, i] - t$b[i])))
  cor(t$resp - E, use = "pairwise.complete.obs")
}

off_diag <- function(m) m[row(m) != col(m)]

test_that("rejects objects that are not IRTnorm", {
  expect_error(item_residual_cor(list(), visual = FALSE), "estimate_norming_model")
  expect_error(item_residual_cor(NULL,   visual = FALSE), "estimate_norming_model")
})

test_that("rejects a domain that does not exist", {
  fit <- make_mock_q3(irt_mod = 5, D = 2)
  expect_error(item_residual_cor(fit, visual = FALSE, domain = 99), "domain")
})

test_that("domain is ignored for unidimensional models with warning", {
  fit <- make_mock_q3(irt_mod = 3)
  expect_warning(item_residual_cor(fit, visual = FALSE, domain = 1), "The domain argument is only for hierarchical norming models")
})

test_that("Q3 matches an independent computation (3PL)", {
  # (expected to FAIL: function returns NULL for irt_mod <= 4, bug 1)
  fit <- make_mock_q3(irt_mod = 3)
  out <- suppressMessages(item_residual_cor(fit, visual = FALSE))
  expect_equal(out$Q3, expected_q3(fit), ignore_attr = TRUE)
})

test_that("Q3 is a proper correlation matrix with item names", {
  # (expected to FAIL: bug 1)
  fit <- make_mock_q3(irt_mod = 3, I = 6)
  Q3  <- suppressMessages(item_residual_cor(fit, visual = FALSE)$Q3)
  expect_equal(dim(Q3), c(6L, 6L))
  expect_equal(dimnames(Q3), list(fit$item_variables, fit$item_variables))
  expect_true(isSymmetric(unname(Q3)))
  expect_equal(unname(diag(Q3)), rep(1, 6))
  expect_true(all(abs(Q3) <= 1 + 1e-12))
})

test_that("summary_info summarises the off-diagonal only", {
  # (expected to FAIL: bug 1)
  fit <- make_mock_q3(irt_mod = 3)
  out <- suppressMessages(item_residual_cor(fit, visual = FALSE))
  expect_equal(out$summary_info, summary(off_diag(out$Q3)))
})

test_that("summary is printed to the console", {
  fit <- make_mock_q3(irt_mod = 3)
  expect_output(item_residual_cor(fit, visual = FALSE), "Min.")
})

test_that("non-item columns in raw_data are ignored", {
  # (expected to FAIL: bug 1)
  fit <- make_mock_q3(irt_mod = 3, I = 5)
  out <- suppressMessages(item_residual_cor(fit, visual = FALSE))
  expect_equal(ncol(out$Q3), 5L)
  expect_false("age" %in% colnames(out$Q3))
})

test_that("missing responses are handled by pairwise deletion", {
  # (expected to FAIL: bug 1)
  fit <- make_mock_q3(irt_mod = 3)
  fit$raw_data[c(1, 5, 9), "item1"] <- NA
  fit$raw_data[c(2, 6),    "item3"] <- NA
  Q3 <- suppressMessages(item_residual_cor(fit, visual = FALSE)$Q3)
  expect_false(anyNA(Q3))
})

test_that("the input object is not modified", {
  fit    <- make_mock_q3(irt_mod = 3)
  before <- fit
  suppressMessages(item_residual_cor(fit, visual = FALSE))
  expect_identical(fit, before)
})


test_that("model 4 computes GRM expected scores from the right thresholds", {
  I <- 4; K_max <- 4
  K_item <- c(3L, 3L, 4L, 3L)
  fit <- make_mock_q3(irt_mod = 4, I = I, J = 80)
  fit$newdata$K_item <- K_item
  fit$newdata$K_max  <- K_max
  
  # true thresholds, deliberately different per item so mix-ups are visible
  thr_true <- list(c(-1, 0), c(-0.5, 0.5), c(-1.5, -0.2, 1.1), c(-0.8, 0.9))
  
  # draws matrix in the real Stan order: all items for k = 1, then k = 2, ...
  grid <- expand.grid(i = seq_len(I), k = seq_len(K_max - 1))
  grid <- grid[grid$k <= K_item[grid$i] - 1, ]
  vals <- mapply(function(i, k) thr_true[[i]][k], grid$i, grid$k)
  thr  <- const_draws(vals)
  colnames(thr) <- sprintf("b[%d,%d]", grid$i, grid$k)
  fit$norm_info$threshold <- thr
  
  # polytomous responses simulated from the GRM
  theta <- fit$model_info$info_sample$theta_est
  a     <- colMeans(fit$norm_info$alpha)
  set.seed(2)
  resp <- sapply(seq_len(I), function(i) {
    cum <- plogis(a[i] * outer(theta, thr_true[[i]], "-"))
    rowSums(matrix(runif(length(cum)), nrow(cum)) < cum)
  })
  colnames(resp) <- fit$item_variables
  fit$raw_data <- data.frame(resp, age = fit$raw_data$age)
  
  # independent expectation (does not use the package code)
  E <- sapply(seq_len(I), function(i) {
    cum <- plogis(a[i] * outer(theta, thr_true[[i]], "-"))
    p   <- cbind(1, cum) - cbind(cum, 0)
    drop(p %*% (seq_len(ncol(p)) - 1))
  })
  expected_Q3 <- cor(resp - E)
  
  fit$model_info$info_item <- list(threshold = thr_true, 
                                        alpha = a)
  out <- NULL
  invisible(capture.output(
    out <- item_residual_cor(fit, visual = FALSE)
  ))
  
  expect_equal(out$Q3, expected_Q3, ignore_attr = TRUE)
})


test_that("multidimensional models run without relying on a global `info`", {
  # (expected to FAIL: bugs 3 and 5)
  skip_if(exists("info", envir = globalenv()),
          "a global `info` masks the bug; remove it and rerun")
  local_mocked_bindings(
    compute_domain_averages = function(Q3, itemD) data.frame(domain = unique(itemD)),
    .package = "IRTnorm"
  )
  for (m in 5:12) {
    fit <- make_mock_q3(irt_mod = m, D = 2, I = 6)
    expect_no_error(item_residual_cor(fit, visual = FALSE))
  }
})

test_that("multidimensional Q3 uses the domain-specific ability per item", {
  # (expected to FAIL: bugs 3 and 5)
  local_mocked_bindings(
    compute_domain_averages = function(Q3, itemD) data.frame(domain = unique(itemD)),
    .package = "IRTnorm"
  )
  fit <- make_mock_q3(irt_mod = 9, D = 3, I = 6)
  out <- suppressMessages(item_residual_cor(fit, visual = FALSE))
  expect_equal(out$Q3, expected_q3(fit), ignore_attr = TRUE)
})

test_that("multidimensional result contains itemD and compare_domains", {
  # (expected to FAIL: bug 1, the first return is overwritten)
  local_mocked_bindings(
    compute_domain_averages = function(Q3, itemD) data.frame(domain = unique(itemD)),
    .package = "IRTnorm"
  )
  fit <- make_mock_q3(irt_mod = 9, D = 2, I = 6)
  out <- suppressMessages(item_residual_cor(fit, visual = FALSE))
  expect_named(out, c("Q3", "itemD", "compare_domains", "summary_info"),
               ignore.order = TRUE)
  expect_equal(out$itemD, fit$newdata$itemD)
})

test_that("compute_domain_averages receives the Q3 matrix and itemD", {
  seen <- NULL
  local_mocked_bindings(
    compute_domain_averages = function(Q3, itemD) {
      seen <<- list(Q3 = Q3, itemD = itemD); data.frame()
    },
    .package = "IRTnorm"
  )
  fit <- make_mock_q3(irt_mod = 5, D = 2, I = 6)
  try(item_residual_cor(fit, visual = FALSE), silent = TRUE)
  expect_equal(seen$itemD, fit$newdata$itemD)
  expect_equal(dim(seen$Q3), c(6L, 6L))
})

test_that("domain restricts Q3 to the items of that domain", {
  # (expected to FAIL: bug 6)
  local_mocked_bindings(
    compute_domain_averages = function(Q3, itemD) data.frame(),
    .package = "IRTnorm"
  )
  fit <- make_mock_q3(irt_mod = 5, D = 2, I = 6)   # itemD = 1 2 1 2 1 2
  out <- suppressMessages(item_residual_cor(fit, visual = FALSE, domain = 1))
  expect_equal(colnames(out$Q3), fit$item_variables[fit$newdata$itemD == 1])
  expect_false(anyNA(out$Q3))
})

test_that("multidimensional models use visual_multi", {
  # (expected to FAIL: bugs 3 and 5)
  called <- FALSE
  local_mocked_bindings(
    compute_domain_averages = function(Q3, itemD) data.frame(),
    visual_multi = function(Q3, itemD) { called <<- TRUE; ggplot2::ggplot() },
    .package = "IRTnorm"
  )
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  item_residual_cor(make_mock_q3(irt_mod = 5, D = 2), visual = TRUE)
  expect_true(called)
})

test_that("the result is returned invisibly", {
  # (expected to FAIL: bug 1)
  fit <- make_mock_q3(irt_mod = 3)
  expect_invisible(item_residual_cor(fit, visual = FALSE))
})