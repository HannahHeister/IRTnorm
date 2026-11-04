
test_that("extract_info wires each helper into its own slot with type = 'mean'", {
  seen <- list()
  local_mocked_bindings(
    extract_iteminfo = function(fit, data, type) {
      seen$item <<- type; "ITEM"
    },
    extract_coefinfo = function(fit, age_mod, irt_mod, type) {
      seen$coef <<- list(age_mod = age_mod, irt_mod = irt_mod, type = type); "COEF"
    },
    extract_personinfo = function(fit, data) "PERSON",
    extract_relationinfo = function(fit, irt_mod, type) {
      seen$rel <<- list(irt_mod = irt_mod, type = type); "REL"
    },
    .package = "IRTnorm"
  )
  
  y_rep_mock <- matrix(1:6, 2, 3)
  fit <- list(draws = function(variables, format) {
    seen$draws <<- list(variables = variables, format = format)
    y_rep_mock
  })
  
  res <- extract_info(fit, list(age_mod = 2, irt_mod = 9))
  
  expect_named(res, c("info_sample", "info_item", "y_rep", "info_coef", "info_relation"))
  expect_equal(res$info_sample, "PERSON")
  expect_equal(res$info_item, "ITEM")
  expect_equal(res$info_coef, "COEF")
  expect_equal(res$info_relation, "REL")
  expect_identical(res$y_rep, y_rep_mock)
  
  expect_equal(seen$item, "mean")
  expect_equal(seen$coef, list(age_mod = 2, irt_mod = 9, type = "mean"))
  expect_equal(seen$rel, list(irt_mod = 9, type = "mean"))
  expect_equal(seen$draws, list(variables = "y_rep", format = "matrix"))
})