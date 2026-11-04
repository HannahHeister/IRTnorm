test_that("models 3 calculates the unidimensional log posterior", {
  theta <- 0.5
  mu <- 0
  sigma <- 1
  
  a_vec <- c(1.2, 0.8, 1.5)
  b_vec <- c(-0.5, 0.2, 1)
  c_vec <- c(0, 0, 0)
  y_vec <- c(1, 0, 1)
  
 irt_mod <- 3
    result <- log_post_draw_theta(
      data = list(),
      theta = theta,
      a_vec = a_vec,
      b_vec = b_vec,
      c_vec = c_vec,
      y_vec = y_vec,
      mu = mu,
      sigma = sigma,
      irt_mod = irt_mod
    )
    
    p <- stats::plogis(a_vec * (theta - b_vec))
    
    expected <- stats::dnorm(
      theta, mean = mu, sd = sigma, log = TRUE
    ) +
      sum(stats::dbinom(y_vec, size = 1, prob = p, log = TRUE))
    
    expect_equal(result, expected)
    
    
})

test_that("NULL and default item parameters give equivalent results", {
  theta <- 0.3
  y_vec <- c(1, 0, 1)
  b_vec <- c(-1, 0, 1)
  
  common <- list(
    data = list(),
    theta = theta,
    b_vec = b_vec,
    y_vec = y_vec,
    mu = 0,
    sigma = 1,
    irt_mod = 2
  )
  
  baseline <- do.call(
    log_post_draw_theta,
    c(common, list(
      a_vec = rep(1, 3),
      c_vec = rep(0, 3)
    ))
  )
  
  null_parameters <- do.call(
    log_post_draw_theta,
    c(common, list(
      a_vec = NULL,
      c_vec = NULL
    ))
  )
  
  zero_guessing <- do.call(
    log_post_draw_theta,
    c(common, list(
      a_vec = rep(1, 3),
      c_vec = rep(0, 3)
    ))
  )
  
  expect_equal(null_parameters, baseline)
  expect_equal(zero_guessing, baseline)
})
test_that("nonzero guessing probabilities affect the likelihood", {
  args <- list(
    data = list(),
    theta = 0.5,
    a_vec = c(1, 1),
    b_vec = c(0, 0),
    y_vec = c(0, 1),
    mu = 0,
    sigma = 1,
    irt_mod = 2
  )
  
  without_guessing <- do.call(
    log_post_draw_theta,
    c(args, list(c_vec = c(0, 0)))
  )
  
  with_guessing <- do.call(
    log_post_draw_theta,
    c(args, list(c_vec = c(0.2, 0.2)))
  )
  
  expect_false(isTRUE(all.equal(
    without_guessing, with_guessing
  )))
})

test_that("model 4 calculates the graded response log posterior", {
  theta <- 0.4
  
  data <- list(
    K_item = c(3, 2), 
    K_max = 3
  )
  
  a_vec <- c(1.2, 0.8)
  b_vec <- c(-1, 0.5,0.2, 0)
  
  y_vec <- c(2, 1)
  
  result <- log_post_draw_theta(
    data = data,
    theta = theta,
    a_vec = a_vec,
    b_vec = b_vec,
    c_vec = NULL,
    obs_K_item = c(3, 2),
    y_vec = y_vec,
    mu = 0,
    sigma = 1,
    irt_mod = 4
  )
  
  y_mat <- make_dummy_matrix(
    info = as.matrix(t(y_vec)),
    K_item = data$K_item
  )
  
  all_p <- c(
    grm_probs(theta, a_vec[1], b_vec[1:2]),
    grm_probs(theta, a_vec[2], b_vec[3])
  )
  
  expected <- stats::dnorm(
    theta, mean = 0, sd = 1, log = TRUE
  ) + sum(y_mat * log(all_p))
  
  expect_equal(result, expected)
  
  
})

test_that("models 5 to 8 calculate the multidimensional log posterior", {
  theta <- c(0.3)
  unique <- c(0.2, 0.5)
  
  mu <- c(0.1, -0.2)
  sigma <- c(1.1, 0.9)
  nu_vec <- c(0.4, -0.3)
  
  a_vec <- c(1.2, 0.8, 1.5)
  b_vec <- c(-0.5, 0.2, 0.7)
  c_vec <- c(0, 0, 0)
  y_vec <- c(1, 0, 1)
  itemD <- c(1, 2, 1)
  
  data <- list(D = 2)
  
  for (irt_mod in 5:8) {
    result <- log_post_draw_theta(
      data = data,
      theta = theta,
      a_vec = a_vec,
      b_vec = b_vec,
      c_vec = c_vec,
      y_vec = y_vec,
      mu = mu,
      sigma = sigma,
      unique = unique,
      nu_vec = nu_vec,
      irt_mod = irt_mod,
      itemD = itemD
    )
    
    delta <- mu +
      nu_vec * sigma * theta +
      sqrt(1 - nu_vec^2) * sigma * unique
    
    p <- stats::plogis(
      a_vec * (delta[itemD] - b_vec)
    )
    
    expected <- sum(stats::dnorm(
      theta, mean = 0, sd = 1, log = TRUE
    )) +
      sum(stats::dnorm(
        unique, mean = 0, sd = 1, log = TRUE
      )) +
      sum(stats::dbinom(
        y_vec, size = 1, prob = p, log = TRUE
      ))
    
    expect_equal(result, expected)
    
    
  }
})

test_that("models 9 and 10 use lambda_vec correctly", {
  theta <- c(0.3)
  unique <- c(0.2, 0.5)
  
  mu <- c(0.1, -0.2)
  sigma <- c(1.1, 0.9)
  nu_vec <- c(0.4, -0.3)
  lambda_vec <- c(0.8, 0.6)
  
  a_vec <- c(1.2, 0.8, 1.5)
  b_vec <- c(-0.5, 0.2, 0.7)
  y_vec <- c(1, 0, 1)
  itemD <- c(1, 2, 1)
  
  for (irt_mod in 9:10) {
    result <- log_post_draw_theta(
      data = list(D = 2),
      theta = theta,
      a_vec = a_vec,
      b_vec = b_vec,
      c_vec = NULL,
      y_vec = y_vec,
      mu = mu,
      sigma = sigma,
      unique = unique,
      nu_vec = nu_vec,
      lambda_vec = lambda_vec,
      irt_mod = irt_mod,
      itemD = itemD
    )
    
    rho <- nu_vec * lambda_vec
    
    delta <- mu +
      rho * sigma * theta +
      sqrt(1 - rho^2) * sigma * unique
    
    p <- stats::plogis(
      a_vec * (delta[itemD] - b_vec)
    )
    
    expected <- sum(stats::dnorm(theta, log = TRUE)) +
      sum(stats::dnorm(unique, log = TRUE)) +
      sum(stats::dbinom(y_vec, size = 1, prob = p, log = TRUE))
    
    expect_equal(result, expected)
    
  }
})

test_that("models 11 and 12 use unique_vec correctly", {
  theta <- c(0.3)
  unique <- c(0.2, 0.5)
  
  mu <- c(0.1, -0.2)
  sigma <- c(1.1, 0.9)
  unique_vec <- c(0.7, 0.4)
  
  a_vec <- c(1.2, 0.8, 1.5)
  b_vec <- c(-0.5, 0.2, 0.7)
  y_vec <- c(1, 0, 1)
  itemD <- c(1, 2, 1)
  
  for (irt_mod in 11:12) {
    result <- log_post_draw_theta(
      data = list(D = 2),
      theta = theta,
      a_vec = a_vec,
      b_vec = b_vec,
      c_vec = NULL,
      y_vec = y_vec,
      mu = mu,
      sigma = sigma,
      unique = unique,
      unique_vec = unique_vec,
      irt_mod = irt_mod,
      itemD = itemD
    )
    
    delta <- mu + sigma * theta + unique_vec * unique
    
    p <- stats::plogis(
      a_vec * (delta[itemD] - b_vec)
    )
    
    expected <- sum(stats::dnorm(theta, log = TRUE)) +
      sum(stats::dnorm(unique, log = TRUE)) +
      sum(stats::dbinom(y_vec, size = 1, prob = p, log = TRUE))
    
    expect_equal(result, expected)
    
  }
})

test_that("extreme binary responses do not produce NaN", {
  result <- log_post_draw_theta(
    data = list(),
    theta = 10,
    a_vec = c(2, 2),
    b_vec = c(-5, 5),
    c_vec = c(0, 0),
    y_vec = c(1, 0),
    mu = 0,
    sigma = 1,
    irt_mod = 2
  )
  
  
  expect_false(is.nan(result))
  expect_false(is.infinite(result) && result > 0)
})