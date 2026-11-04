
vspec <- function(irt = "1PLnorm", age = "poly", ...) {
  c(list(irt_model = irt, age_model = age, poly_mean = 2, poly_sd = 3), list(...))
}

test_that("valid input returns TRUE invisibly", {
  x <- list(A = vspec(), B = vspec("2PLnorm"))
  expect_invisible(validate_model_specifications_list(x))
  expect_true(validate_model_specifications_list(x))
})

test_that("extra fields in a specification are allowed", {
  x <- list(A = c(vspec(), poly_mean = 1), B = vspec("2PLnorm"))
  expect_true(validate_model_specifications_list(x))
})

test_that("non-list input errors", {
  expect_error(validate_model_specifications_list("A"), "must be a list")
})

test_that("fewer than two models errors and reports the count", {
  expect_error(validate_model_specifications_list(list()),
               "at least two.*has 0")
  expect_error(validate_model_specifications_list(list(A = vspec())),
               "at least two.*has 1")
})

test_that("a single flat specification passed by mistake is caught", {
  # length 2, but its elements are character strings, not lists
  expect_error(validate_model_specifications_list(vspec()),
               "must itself be a list.*irt_model, age_model")
})

test_that("non-list element is reported by name or by position", {
  expect_error(
    validate_model_specifications_list(list(A = vspec(), B = "oops")),
    "Offending element\\(s\\): B"
  )
  expect_error(
    validate_model_specifications_list(list(vspec(), "oops")),
    "Offending element\\(s\\): 2"
  )
})

test_that("missing irt_model or age_model is reported per entry", {
  expect_error(
    validate_model_specifications_list(
      list(A = vspec(), B = list(irt_model = "2PLnorm"))),
    "B: missing age_model"
  )
  expect_error(
    validate_model_specifications_list(
      list(A = list(age_model = "poly"), B = vspec())),
    "A: missing irt_model"
  )
  expect_error(
    validate_model_specifications_list(list(A = vspec(), B = list())),
    "B: missing irt_model, age_model"
  )
})

test_that("explicit NULL values count as missing", {
  x <- list(A = vspec(), B = list(irt_model = NULL, age_model = "poly"))
  # list(irt_model = NULL) keeps the entry, but is.null() is still TRUE
  expect_error(validate_model_specifications_list(x), "B: missing irt_model")
})

test_that("unnamed entries are labelled by position", {
  expect_error(
    validate_model_specifications_list(list(vspec(), list(irt_model = "2PLnorm"))),
    "\\[2\\]: missing age_model"
  )
})

test_that("partially named lists use position for the unnamed entry", {
  x <- list(A = vspec(), list(irt_model = "2PLnorm"))
  expect_error(validate_model_specifications_list(x), "\\[2\\]: missing age_model")
})

test_that("all problems are reported together", {
  x <- list(A = list(), B = list(irt_model = "2PLnorm"))
  expect_error(validate_model_specifications_list(x),
               "A: missing irt_model, age_model\nB: missing age_model")
})