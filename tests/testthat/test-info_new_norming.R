
test_that("info_new_norming forwards args with type = 'draws' and returns blocks in order", {
  seen <- list()
  local_mocked_bindings(
    extract_coefinfo = function(fit, age_mod, irt_mod, type) {
      seen$coef <<- list(age_mod = age_mod, irt_mod = irt_mod, type = type)
      list(coef_a = "COEF")                       # replace with your real names
    },
    extract_iteminfo = function(fit, data, type) {
      seen$item <<- type
      list(alpha = "ITEM", beta = "ITEM", gamma = "ITEM")
    },
    extract_relationinfo = function(fit, irt_mod, type) {
      seen$rel <<- list(irt_mod = irt_mod, type = type)
      list(nu = "REL", lambda = NULL, unique_sigma = NULL)
    },
    .package = "IRTnorm"
  )
  
  res <- info_new_norming(fit = "fake_fit", data = list(age_mod = 2, irt_mod = 9))
  
  # arguments forwarded, always in "draws" mode
  expect_equal(seen$coef, list(age_mod = 2, irt_mod = 9, type = "draws"))
  expect_equal(seen$item, "draws")
  expect_equal(seen$rel, list(irt_mod = 9, type = "draws"))
  
  # flat list, blocks ordered coef -> item -> relation, no duplicate names
  expect_identical(names(res),
                   c("coef_a", "alpha", "beta", "gamma", "nu", "lambda", "unique_sigma"))
  expect_false(anyDuplicated(names(res)) > 0)
  
  # each value comes from the right helper
  expect_equal(res$coef_a, "COEF")
  expect_equal(res$alpha, "ITEM")
  expect_equal(res$nu, "REL")
  expect_null(res$lambda)
})