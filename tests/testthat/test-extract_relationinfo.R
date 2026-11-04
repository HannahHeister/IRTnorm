
test_that("extract_relationinfo returns the right components per irt_mod", {
  local_mocked_bindings(
    get_value = function(fit, var, type) paste(var, type),   # recognisable marker
    .package = "IRTnorm"
  )
  
  expected <- list(            # which components should be filled
    `1`  = character(), `2`  = character(), `3`  = character(), `4` = character(),
    `5`  = "nu", `6`  = "nu", `7`  = "nu", `8` = "nu",
    `9`  = c("nu", "lambda"), `10` = c("nu", "lambda"),
    `11` = "unique_sigma", `12` = "unique_sigma"
  )
  
  for (m in names(expected)) {
    res <- extract_relationinfo(fit = NULL, irt_mod = as.integer(m), type = "mean")
    expect_named(res, c("nu", "lambda", "unique_sigma"), info = m)   # always all three slots
    filled <- names(Filter(Negate(is.null), res))
    expect_setequal(filled, expected[[m]])
  }
})

test_that("type is passed through to get_value and validated", {
  local_mocked_bindings(get_value = function(fit, var, type) type, .package = "IRTnorm")
  expect_equal(extract_relationinfo(NULL, 9, "draws")$lambda, "draws")
  expect_equal(extract_relationinfo(NULL, 9)$nu, "mean")   # default is "mean"
  expect_error(extract_relationinfo(NULL, 9, "bogus"))
})