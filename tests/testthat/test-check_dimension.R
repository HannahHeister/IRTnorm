
uni   <- function(m = "2PLnorm")  list(irt_model = m, age_model = "poly")
multi <- function(m = "HO1PLnorm-C") list(irt_model = m, age_model = "poly")

test_that("all unidimensional models pass", {
  expect_invisible(check_dimension(list(uni("1PLnorm"), uni("2PLnorm"), uni("GRMnorm"))))
  expect_true(check_dimension(list(uni(), uni("3PLnorm"))))
})

test_that("all multidimensional models pass", {
  expect_true(check_dimension(list(multi(), multi("HO1PLnorm-C"))))
})

# EXPECTED TO FAIL with current code: `unidim <- ...` overwrites the whole
# vector with a scalar each iteration, so only the LAST model is inspected
# and any(unidim) & any(!unidim) can never both be TRUE.
# Fix: `unidim[n] <- model_specifications[[n]]$irt_model %in% c(...)`
test_that("mixing unidimensional and multidimensional models errors", {
  expect_error(check_dimension(list(uni(), multi())), "not imply the same dimensionality")
  expect_error(check_dimension(list(multi(), uni())), "not imply the same dimensionality")
  expect_error(check_dimension(list(uni(), uni("1PLnorm"), multi())),
               "not imply the same dimensionality")
})