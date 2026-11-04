test_that("generated Stan code passes the syntax check for every model, age model and scenario", {
  skip_on_cran(); skip_if_not(cmdstan_ok(), "cmdstan not available")
  for (x in combos())
    expect_true(syntax_ok(build_code(x$m, x$a, x$scen)), info = paste(x$m, x$a, x$s))
})

test_that("every data variable the program declares is supplied by data_prep()", {
  skip_on_cran(); skip_if_not(cmdstan_ok())
  local_mocked_bindings(build_spline = fake_spline, .package = "IRTnorm")
  for (m in ALL_MODELS) for (a in AGE_MODELS) {
    declared <- stan_vars(build_code(m, a, list()))$data
    expect_identical(setdiff(declared, names(suppressWarnings(prep_data(m, a)))), character(0),
                     info = paste(m, a))
  }
})

test_that("the R extractors can find the Stan variables they read", {
  skip_on_cran(); skip_if_not(cmdstan_ok())
  for (m in ALL_MODELS) for (a in AGE_MODELS) {
    v <- stan_vars(build_code(m, a, list())); i <- match(m, ALL_MODELS)
    info <- paste(m, a)
    # extract_personinfo / extract_info
    expect_true("z_theta" %in% v$parameters, info = info)
    expect_true(all(c("theta_mean", "theta_sigma") %in% v$transformed_parameters), info = info)
    expect_true(all(c("y_rep", "log_lik") %in% v$generated_quantities), info = info)
    expect_true(if (i <= 4) "theta" %in% v$transformed_parameters
                else "substraits" %in% v$transformed_parameters, info = info)
    # extract_iteminfo: estimated exactly where it does not default
    expect_equal("alpha" %in% v$parameters, !i %in% c(1, 5, 7, 9, 11), info = info)
    expect_equal("gamma" %in% v$parameters, i == 3,                    info = info)
    expect_equal("b"     %in% v$parameters, i == 4,                    info = info)
    expect_equal("beta"  %in% v$parameters, i != 4,                    info = info)
    # extract_relationinfo
    expect_equal("nu" %in% v$parameters, i %in% 5:10, info = info)
    expect_equal("lambda" %in% c(v$parameters, v$transformed_parameters), i %in% 9:10, info = info)
    expect_equal("unique_sigma" %in% v$transformed_parameters, i %in% 11:12, info = info)
  }
})

test_that("every sampled parameter gets a prior, and every model has a likelihood", {
  skip_on_cran(); skip_if_not(cmdstan_ok())
  for (x in combos()) {
    code <- build_code(x$m, x$a, x$scen)
    info <- paste(x$m, x$a, x$s)
    expect_identical(setdiff(stan_vars(code)$parameters, model_lhs(code)),
                     character(0), info = info)
    expect_true(grepl("\\by(\\[n\\])?\\s*~|target\\s*\\+=", model_section(code)), info = info)
  }
})

test_that("supplied priors are used, and parameters supplied as data are never sampled", {
  skip_on_cran(); skip_if_not(cmdstan_ok())
  for (x in combos()) {
    code <- build_code(x$m, x$a, x$scen)
    v <- stan_vars(code); info <- paste(x$m, x$a, x$s)
    prior_vars <- grep("_prior_(mu|sd)$", v$data, value = TRUE)
    for (p in prior_vars)
      expect_true(grepl(p, model_section(code), fixed = TRUE), info = paste(info, p))
    expect_identical(setdiff(intersect(v$data, model_lhs(code)), "y"),
                     character(0), info = info)
  }
})

test_that("every prior requested in prior_knowledge is declared as data", {
  skip_on_cran(); skip_if_not(cmdstan_ok())
  for (m in ALL_MODELS) for (a in AGE_MODELS) {
    asked <- scenarios(m)$prior_all$prior_knowledge
    v <- stan_vars(build_code(m, a, list(prior_knowledge = asked)))$data
    for (p in asked)
      expect_true(paste0(p, "_prior_mu") %in% v, info = paste(m, a, p))
  }
})

test_that("spline random effect u_k is scaled by its own tau_k", {
  for (m in ALL_MODELS) {
    ms <- regmatches(build_code(m, "splines", list()),
                     gregexpr("\\bu([123])(\\[d\\])?\\s*~\\s*normal\\(0,\\s*tau([123])", 
                              build_code(m, "splines", list())))[[1]]
    expect_gt(length(ms), 0)
    expect_identical(sub("^u([123]).*", "\\1", ms), sub(".*tau([123])$", "\\1", ms), info = m)
  }
})

test_that("only the 3PL model gets the custom functions block", {
  for (m in ALL_MODELS)
    expect_equal(grepl("logit_3pl_lpmf", build_code(m, "poly", list())), m == "3PLnorm", info = m)
})

