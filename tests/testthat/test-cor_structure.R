make_mock_irtnorm <- function(irt_mod = 4, D = 3, items_per_domain = 4,
                              n = 50, seed = 1) {
  set.seed(seed)
  itemD   <- rep(seq_len(D), each = items_per_domain)
  n_items <- length(itemD)
  
  items <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n, ncol = n_items)
  colnames(items) <- paste0("item", seq_len(n_items))
  raw_data <- data.frame(items, age = runif(n, 5, 15))
  
  structure(
    list(
      newdata = list(irt_mod = irt_mod, itemD = itemD, D = D, J = n),
      raw_data = raw_data,
      item_variables = colnames(items),
      age_variable = "age",
      model_info = list(info_relation = list(
        nu           = seq(0.5, 0.9, length.out = D),
        lambda       = runif(n, 0.5, 1),
        unique_sigma = runif(D, 0.1, 0.5)
      )),
      info_sample = data.frame(
        trait       = rep(paste0("delta", seq_len(D)), each = n),
        theta_sigma = runif(n * D)
      )
    ),
    class = "IRTnorm"
  )
}

pair_labels <- function(D) {
  m <- which(upper.tri(diag(D)), arr.ind = TRUE)
  paste0(m[, "row"], "&", m[, "col"])
}

test_that("rejects objects that are not IRTnorm", {
  expect_error(cor_structure(list()), "estimate_norming_model")
  expect_error(cor_structure(NULL),   "estimate_norming_model")
})

test_that("rejects unidimensional models", {
  for (m in 1:3) {
    expect_error(cor_structure(make_mock_irtnorm(irt_mod = m)),
                 "hierarchical")
  }
})

test_that("models 4-8 return the expected structure", {
  for (m in 4:8) {
    out <- cor_structure(make_mock_irtnorm(irt_mod = m, D = 3))
    expect_s3_class(out, "data.frame")
    expect_named(out, c("sumscore_cor", "model_cor", "domains"))
    expect_equal(nrow(out), choose(3, 2))
    expect_type(out$domains, "character")
  }
})

test_that("model_cor is the product of the domain loadings", {
  fit <- make_mock_irtnorm(irt_mod = 5, D = 4)
  nu  <- fit$model_info$info_relation$nu
  out <- cor_structure(fit)
  
  ij <- do.call(rbind, strsplit(out$domains, "&"))
  expected <- nu[as.integer(ij[, 1])] * nu[as.integer(ij[, 2])]
  expect_equal(out$model_cor, expected)
})

test_that("sumscore_cor matches the domain label in its own row", {
  # Compares by label, so it is robust to row ordering. Expected to FAIL for D >= 4.
  for (D in 2:5) {
    fit  <- make_mock_irtnorm(irt_mod = 4, D = D)
    out  <- cor_structure(fit)
    
    items <- as.matrix(fit$raw_data[, fit$item_variables])
    sums  <- sapply(seq_len(D), function(d)
      rowSums(items[, fit$newdata$itemD == d, drop = FALSE], na.rm = TRUE))
    cmat  <- cor(sums)
    
    ij <- do.call(rbind, strsplit(out$domains, "&"))
    expected <- cmat[cbind(as.integer(ij[, 1]), as.integer(ij[, 2]))]
    expect_equal(out$sumscore_cor, expected,
                 info = paste("D =", D))
  }
})

test_that("domain labels cover every pair exactly once", {
  for (D in 2:5) {
    out <- cor_structure(make_mock_irtnorm(irt_mod = 4, D = D))
    expect_setequal(out$domains, pair_labels(D))
    expect_false(anyDuplicated(out$domains) > 0)
  }
})

test_that("missing item responses are ignored in sum scores", {
  fit <- make_mock_irtnorm(irt_mod = 4, D = 3)
  fit$raw_data[1:5, fit$item_variables[1:2]] <- NA
  out <- cor_structure(fit)
  expect_false(anyNA(out$sumscore_cor))
})

test_that("models 9-12 return one column per domain pair plus age", {
  for (m in 9:12) for (D in 2:4) {
    fit <- make_mock_irtnorm(irt_mod = m, D = D)
    out <- cor_structure(fit, visual = FALSE)
    
    expect_s3_class(out, "data.frame")
    expect_equal(ncol(out), choose(D, 2) + 1, info = paste("m =", m, "D =", D))
    expect_equal(nrow(out), nrow(fit$raw_data))
    expect_equal(colnames(out)[ncol(out)], "age")
    expect_equal(out$age, fit$raw_data$age)
  }
})

test_that("column names use the same pairs as the sum-score output", {
  fit <- make_mock_irtnorm(irt_mod = 9, D = 3)
  out <- cor_structure(fit, visual = FALSE)
  expect_equal(length(setdiff(colnames(out), "age")), 3)
})

test_that("models 11/12: delta sigma is scaled per domain, not recycled", {
  # With J != D, recycling `delta_sigma * unique_sigma` multiplies the wrong
  # elements. Check one column by hand via the internal step.
  fit <- make_mock_irtnorm(irt_mod = 11, D = 3, n = 50)
  us  <- fit$model_info$info_relation$unique_sigma
  ds  <- sapply(1:3, function(d)
    fit$info_sample$theta_sigma[fit$info_sample$trait == paste0("delta", d)])
  correct <- sweep(ds, 2, us, "*")
  recycled <- ds * us
  expect_false(isTRUE(all.equal(correct, recycled)))  # documents the hazard
})

test_that("visual = FALSE does not draw and returns the same data", {
  fit <- make_mock_irtnorm(irt_mod = 9, D = 3)
  expect_equal(cor_structure(fit, visual = FALSE),
               suppressWarnings(cor_structure(fit, visual = TRUE)))
})

test_that("visual = TRUE runs without attaching ggplot2", {
  skip_if_not_installed("ggplot2")
  fit <- make_mock_irtnorm(irt_mod = 9, D = 3)
  pdf(NULL); on.exit(dev.off())
  expect_no_error(cor_structure(fit, visual = TRUE))   # fails: bare `ylab()`
})

test_that("plot labels line up with their values", {
  skip_if_not_installed("ggplot2")
  fit <- make_mock_irtnorm(irt_mod = 9, D = 3)
  pdf(NULL); on.exit(dev.off())
  out  <- cor_structure(fit, visual = TRUE)
  long <- ggplot2::last_plot()$data
  
  for (lab in setdiff(colnames(out), "age")) {
    expect_equal(long$cor_value[long$combi == lab], out[[lab]], ignore_attr = TRUE)
  }
})