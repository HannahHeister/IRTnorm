
make_age_data <- function(
    age_mod = 1,
    S = 2,
    W = 2,
    D = 1,
    M = 1
) {
  list(
    age_mod = age_mod,
    S = S,
    W = W,
    D = D,
    W = W,
    M = M
  )
}

make_post_draw <- function(
    p1 = matrix(c(1, 2, 3, 4), nrow = 2),
    p2 = matrix(c(0, 0.1, 0.2, 0.3), nrow = 2),
    p3 = NULL
) {
  list(
    p1_theta = p1,
    p2_theta = p2,
    p3_theta = p3
  )
}

test_that("age_mod = 1 computes correct unidimensional effects", {
  age_scale <- c(0, 1, 2)
  
  post_draw <- list(
    p1_theta = matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      byrow = TRUE
    ),
    p2_theta = matrix(
      c(0, 0, log(2), 0),
      nrow = 2,
      byrow = TRUE
    )
  )
  
  data <- make_age_data(S = 2, W = 2)
  
  result <- compute_age_effects(
    post_draw = post_draw,
    data = data,
    age_scale = age_scale,
    nperson = length(age_scale),
    irt_mod = 5,
    multidim = FALSE
  )
  
  expected_mu <- matrix(
    c(1, 3, 5, 3, 7, 11),
    nrow = 2,
    byrow = TRUE
  )
  
  expected_sigma <- exp(
    post_draw$p2_theta %*%
      t(model.matrix(
        ~ poly(age, 1, raw = TRUE),
        data = data.frame(age = age_scale)
      ))
  )
  
  expect_equal(unname(result$mu_post), expected_mu)
  expect_equal(result$sigma_post, expected_sigma)
  expect_null(result$lambda_post)
  expect_null(result$unique_post)
})

test_that("single-coefficient effects are constant over age", {
  age_scale <- c(10, 20, 30, 40)
  
  post_draw <- list(
    p1_theta = matrix(c(2, 5), nrow = 2),
    p2_theta = matrix(c(0, log(2)), nrow = 2)
  )
  
  data <- make_age_data(S = 1, W = 1)
  
  result <- compute_age_effects(
    post_draw = post_draw,
    data = data,
    age_scale = age_scale,
    nperson = length(age_scale),
    irt_mod = 5,
    multidim = FALSE
  )
  
  expect_equal(
    result$mu_post,
    matrix(
      c(2, 2, 2, 2, 5, 5, 5, 5),
      nrow = 2,
      byrow = TRUE
    )
  )
  
  expect_equal(
    result$sigma_post,
    matrix(
      c(1, 1, 1, 1, 2, 2, 2, 2),
      nrow = 2,
      byrow = TRUE
    )
  )
})

test_that("models 5 and 6 repeat effects across dimensions", {
  age_scale <- c(0, 1, 2)
  D <- 3
  
  post_draw <- list(
    p1_theta = matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE),
    p2_theta = matrix(c(0, 0.1, 0.2, 0.3), nrow = 2, byrow = TRUE)
  )
  
  data <- make_age_data(S = 2, W = 2, D = D)
  
  result <- compute_age_effects(
    post_draw = post_draw,
    data = data,
    age_scale = age_scale,
    nperson = length(age_scale),
    irt_mod = 5,
    multidim = TRUE
  )
  
  expect_equal(dim(result$mu_post), c(2, 3, D))
  expect_equal(dim(result$sigma_post), c(2, 3, D))
  
  for (d in seq_len(D)) {
    expect_equal(
      result$mu_post[, , d],
      result$mu_post[, , 1]
    )
    
    expect_equal(
      result$sigma_post[, , d],
      result$sigma_post[, , 1]
    )
  }
  
  expect_null(result$lambda_post)
  expect_null(result$unique_post)
})

test_that("models 9 and 10 return lambda_post", {
  age_scale <- c(0, 1, 2)
  
  post_draw <- list(
    p1_theta = matrix(c(1, 2), nrow = 2),
    p2_theta = matrix(c(0, 0.1), nrow = 2),
    p3_theta = matrix(c(0.5, 1), nrow = 2)  )
  
  data <- make_age_data(S = 1, W = 1, M = 1)
  
  result <- compute_age_effects(
    post_draw = post_draw,
    data = data,
    age_scale = age_scale,
    nperson = length(age_scale),
    irt_mod = 9,
    multidim = TRUE
  )
  
  expect_type(result, "list")
  expect_named(
    result,
    c("mu_post", "sigma_post", "lambda_post", "unique_post")
  )
  
  expect_equal(
    result$lambda_post,
    matrix(plogis(c(0.5, 0.5, 0.5, 1, 1, 1)), nrow = 2, byrow = TRUE)
  )
  
  expect_null(result$unique_post)
})

test_that("models 11 and 12 return positive unique effects", {
  age_scale <- c(0, 1, 2)
  D <- 2
  
  post_draw <- list(
    p1_theta = matrix(c(1, 2, 3, 4), nrow = 2),
    p2_theta = matrix(c(0, 0.1, 0.2, 0.3), nrow = 2),
    p3_theta = matrix(
      c(0.5, 1, 1.5, 2),
      nrow = 2,
      byrow = TRUE
    )
  )
  
  data <- make_age_data(
    S = 1, W = 1, M = 1, D = D
  )
  
  result <- compute_age_effects(
    post_draw = post_draw,
    data = data,
    age_scale = age_scale,
    nperson = length(age_scale),
    irt_mod = 11,
    multidim = TRUE
  )
  
  expect_null(result$lambda_post)
  expect_equal(dim(result$unique_post), c(2, 3, D))
  expect_true(all(result$unique_post > 0))
  
  for (d in seq_len(D)) {
    expected <- matrix(
      log1p(exp(post_draw$p3_theta[, d])),
      nrow = 2,
      ncol = length(age_scale)
    )
    
    expect_equal(result$unique_post[, , d], expected)
  }
})
test_that("output dimensions match the requested age scale", {
  age_scale <- c(20, 30, 40, 50, 60)
  
  post_draw <- list(
    p1_theta = matrix(c(1, 2), nrow = 2),
    p2_theta = matrix(c(0, 0.1), nrow = 2)
  )
  
  data <- make_age_data(S = 1, W = 1)
  
  result <- compute_age_effects(
    post_draw = post_draw,
    data = data,
    age_scale = age_scale,
    nperson = length(age_scale),
    irt_mod = 5,
    multidim = FALSE
  )
  
  expect_equal(dim(result$mu_post), c(2, 5))
  expect_equal(dim(result$sigma_post), c(2, 5))
})