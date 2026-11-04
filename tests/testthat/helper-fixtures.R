ITEMS   <- paste0("i", 1:4)
ITEM_D4 <- c(1, 1, 2, 2)          # item -> domain, for the multidimensional models

# 6 persons x 4 items, binary, with missing values
make_raw <- function() data.frame(
  age = c(5, 6, 7, 8, 9, 10),
  i1 = c(1, 0, 1, NA, 1, 0), i2 = c(0, 1, NA, 1, 0, 1),
  i3 = c(1, 1, 0, 0, NA, 1), i4 = c(0, NA, 1, 1, 1, 0))

# same layout, 3 response categories coded 1..3 (for the GRM)
make_raw_grm <- function() data.frame(
  age = c(5, 6, 7, 8, 9, 10),
  i1 = c(1, 2, 3, NA, 1, 2), i2 = c(2, 3, NA, 1, 2, 3),
  i3 = c(1, 1, 2, 3, NA, 3), i4 = c(3, NA, 1, 2, 2, 1))

# Stand-in for build_spline(): right structure, dummy values, no mgcv.
# Use with local_mocked_bindings(build_spline = fake_spline, .package = PKG)
fake_spline <- function(age_df, k, fixed_knots) {
  J <- nrow(age_df)
  list(S_fixed = 2,
       K       = k - 2,
       X       = matrix(1, J, 2),
       Z       = matrix(2, J, k - 2),
       R_Z     = matrix(3, k, k - 2),
       R_X     = matrix(4, k, 2),
       sm      = paste0("sm_k", k))     # marker: shows which k produced it
}

# data_prep() with sensible defaults per model; extra args override (NULL removes)
prep_data <- function(m, age_model = "poly", raw = make_raw(), ...) {
  args <- list(raw_data = raw, age_variable = "age", item_variables = ITEMS,
               irt_model = m, age_model = age_model)
  if (age_model == "poly")  args <- c(args, list(poly_mean = 2, poly_sd = 1))
  if (m %in% MULTI_MODELS)  args$itemD  <- ITEM_D4
  if (m == "GRMnorm")       args$K_item <- rep(2, length(ITEMS))
  args <- modifyList(args, list(...))
  utils::capture.output(res <- do.call(data_prep, args))   # silences print()
  res
}


raw <- data.frame(age = c(1, 2, 3, 4), item1 = c(0, 1, 0, 1), item2 = c(1, 1, 0, 0))
