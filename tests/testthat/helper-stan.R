cmdstan_ok <- function()
  requireNamespace("cmdstanr", quietly = TRUE) &&
  !is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL))

stan_model_obj <- function(code)
  cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), compile = FALSE)

syntax_ok <- function(code)
  isTRUE(tryCatch(stan_model_obj(code)$check_syntax(quiet = TRUE),
                  error = function(e) FALSE))

# list(data, parameters, transformed_parameters, generated_quantities) -> variable names
# (needs a recent cmdstanr with $variables() on uncompiled models)
stan_vars <- function(code) lapply(stan_model_obj(code)$variables(), names)

# the model { ... } section of the generated program
model_section <- function(code)
  regmatches(code, regexpr("(?s)\\bmodel\\s*\\{.*?\\}\\s*generated quantities",
                           code, perl = TRUE))

# names on the left of `~` in the model block
model_lhs <- function(code) {
  sec <- model_section(code)
  x <- regmatches(sec, gregexpr("[A-Za-z_][A-Za-z0-9_]*(\\[[^]]*\\])?\\s*~", sec))[[1]]
  unique(sub("\\[.*|\\s*~", "", x))
}

# Item / relation parameters a model samples (matches the flags in create_stan_data)
used_params <- function(m) c(
  if (m == "GRMnorm") c("b", "alpha") else "beta",
  if (m %in% c("2PLnorm", "3PLnorm", "Testlet2PLnorm", "HO2PLnorm-C",
               "HO2PLnorm-P", "HO2PLnorm-F")) "alpha",
  if (m == "3PLnorm") "gamma",
  if (m %in% ALL_MODELS[5:10]) "nu")

# prior_knowledge / parameter_fixed as fit_IRTnorm builds them (character vectors)
scenarios <- function(m) {
  p     <- used_params(m)
  items <- intersect(p, c("alpha", "beta", "gamma", "b"))
  first <- if (m == "GRMnorm") "b" else "beta"
  list(default   = list(),
       fixed_one = list(parameter_fixed = first),
       fixed_all = list(parameter_fixed = items),
       prior_all = list(prior_knowledge = p),
       mixed     = list(parameter_fixed = first,
                        prior_knowledge = setdiff(p, first)))
}

build_code <- function(m, a, scen = list())
  do.call(create_stan_data, c(list(irt_model = m, age_model = a), scen))

# every model x age model x scenario, so all tests loop the same way
combos <- function() {
  out <- list()
  for (m in ALL_MODELS) for (a in AGE_MODELS) for (s in names(scenarios(m)))
    out[[length(out) + 1]] <- list(m = m, a = a, s = s, scen = scenarios(m)[[s]])
  out
}