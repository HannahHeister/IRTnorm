make_fit_info_known <- function(mus, sigmas, ages) {
  n_persons <- length(mus)
  stopifnot(length(sigmas) == n_persons, length(ages) == n_persons)
  list(model_info = list(
    info_sample = data.frame(
      irt_mod       = 1,
      age         = ages,
      theta_mu    = mus,
      theta_sigma = sigmas,
      theta_est   = mus   # use mu as a stand-in for point estimate
    )
  ))
}


test_that("quantile values exactly match qnorm for known mu and sigma", {
  mus    <- c(0, 1, -1,0, 1, -1,0, 1, -1)
  sigmas <- c(1, 0.5, 2,1, 0.5, 2,1, 0.5, 2)
  ages   <- c(7, 10, 14,7, 10, 14,7, 10, 14)
  probs <- c(0.05, 0.5, 0.95)
  nprobs <- length(probs)
  fit_info <- make_fit_info_known(mus, sigmas, ages)
  
  suppressWarnings(
    result <- quantile_curves(fit_info = fit_info, perspective = c("ability", "rawscore"), 
                                          probs =  probs)
  )
  

  
  for (i in seq_along(mus)) {
    expected <- qnorm(probs, mean = mus[i], sd = sigmas[i])
    # Rows for person i are at positions (i-1)*7 + 1 : i*7
    actual <- result$quan[((i - 1) * nprobs + 1):(i * nprobs)]
    
    expect_equal(actual, expected, tolerance = 1e-10,
                 info = paste("quantile values wrong for person", i,
                              "with mu =", mus[i], "sigma =", sigmas[i]))
  }
})

test_that("median quantile (p0.5) equals theta_mu for each person", {
  # qnorm(0.5, mean = mu, sd = sigma) == mu exactly
  mus    <- c(-2, 0, 2)
  sigmas <- c(1, 1, 1)
  ages   <- c(7, 10, 14)
  
  fit_info <- make_fit_info_known(mus, sigmas, ages)
  
  suppressWarnings(
    result <- quantile_curves(fit_info = fit_info, perspective = "ability")
  )
  
  median_rows <- result[result$qname == "p0.5", ]
  
  expect_equal(median_rows$quan, mus, tolerance = 1e-10,
               info = "p0.5 quantile should equal theta_mu for each person")
})


test_that("qname cycles through all 7 quantile labels per person, not sequentially", {
  # With 2 persons, qname should be:
  # p0.05, p0.1, p0.25, p0.5, p0.75, p0.9, p0.95,
  # p0.05, p0.1, p0.25, p0.5, p0.75, p0.9, p0.95
  # NOT: p0.05 repeated 14 times
  fit_info <- make_fit_info_known(
    mus    = c(0, 1),
    sigmas = c(1, 1),
    ages   = c(8, 12)
  )
  
  suppressWarnings(
    result <- quantile_curves(fit_info = fit_info, perspective = "ability")
  )
  
  expected_qnames <- rep(paste0("p", c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)), 2)
  
  expect_equal(result$qname, expected_qnames,
               info = "qname should cycle through all 7 labels per person")
})

test_that("multidim: quantile values match qnorm for known mu and sigma per dimension", {
  # Give each dimension a distinct mu so we can verify dimension-level correctness
  n_persons <- 3
  ages      <- c(7, 10, 14)
  
  info_sample <- rbind(
    data.frame(irt_mod = 5, age = ages, trait = "delta1",
               theta_mu = c(0, 0, 0), theta_sigma = c(1, 1, 1),
               est_ability = c(0, 0, 0)),
    data.frame(irt_mod = 5, age = ages, trait = "delta2",
               theta_mu = c(2, 2, 2), theta_sigma = c(1, 1, 1),
               est_ability = c(2, 2, 2))
  )
  fit_info <- list(model_info = list(info_sample = info_sample), 
                   newdata = list(D = 2, 
                                  J = 3))
  
  suppressWarnings(
    result <- quantile_curves(fit_info = fit_info, perspective = "ability")
  )
  
  probs <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
  
  # delta1 rows: mu = 0, so median should be 0
  delta1_median <- result$data$quan[
    result$data$domain == "delta1" & result$data$qname == "p0.5"
  ]
  expect_true(all(abs(delta1_median - 0) < 1e-10),
              info = "delta1 median quantile should equal 0")
  
  # delta2 rows: mu = 2, so median should be 2
  delta2_median <- result$data$quan[
    result$data$domain == "delta2" & result$data$qname == "p0.5"
  ]
  expect_true(all(abs(delta2_median - 2) < 1e-10),
              info = "delta2 median quantile should equal 2")
})

