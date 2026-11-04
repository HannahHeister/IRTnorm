

.DIAG_FUNS <- list(
  rhat     = posterior::rhat,
  ess_bulk = posterior::ess_bulk,
  ess_tail = posterior::ess_tail
)

.SUMMARY_NAMES <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")


# ---- Helpers ---------------------------------------------------------------
assert_cmdstanr <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package 'cmdstanr' is required but not installed.\n",
      "Install it with: install.packages('cmdstanr', ",
      "repos = c('https://stan-dev.r-universe.dev', getOption('repos')))",
      call. = FALSE
    )
  }
}

validate_info_inputs <- function(fit, data) {
  if (!is.list(data)) {
    stop("`data` must be a list.", call. = FALSE)
  }
  if (!is.function(fit$draws)) {
    stop("`fit` must have a `$draws()` method (e.g. a CmdStanMCMC object).",
         call. = FALSE)
  }
  check_choice <- function(x, name, allowed) {
    if (!is.numeric(x) || length(x) != 1L || !x %in% allowed) {
      stop(sprintf("`data$%s` must be a single number in {%s}.",
                   name, paste(allowed, collapse = ", ")),
           call. = FALSE)
    }
  }
  check_choice(data$age_mod, "age_mod", 1:2)
  check_choice(data$irt_mod, "irt_mod", 1:12)
  invisible(TRUE)
}

# Pure function: which parameters should be diagnosed for this model?
get_info_params <- function(data) {
  has_third <- data$irt_mod > 8
  
  param1 <- switch(
    as.character(data$age_mod),
    "1" = c("p1_theta", "p2_theta", if (has_third) "p3_theta"),
    "2" = c("l1", "l2", "u1", "u2", if (has_third) c("l3", "u3"))
  )
  
  irt_mod_params <- list(
    "1"  = c("z_theta", param1, "beta"),
    "2"  = c("z_theta", param1, "alpha", "beta"),
    "3"  = c("z_theta", param1, "alpha", "beta", "gamma"),
    "4"  = c("z_theta", param1, "alpha", "b"),
    "5"  = c("z_theta", "unique", "nu", param1, "beta"),
    "6"  = c("z_theta", "unique", "nu", param1, "alpha", "beta"),
    "7"  = c("z_theta", "unique", "nu", param1, "beta"),
    "8"  = c("z_theta", "unique", "nu", param1, "alpha", "beta"),
    "9"  = c("z_theta", "unique", "nu", "lambda", param1, "alpha", "beta"),
    "10" = c("z_theta", "unique", "nu", "lambda", param1, "alpha", "beta"),
    "11" = c("z_theta", "unique", "unique_sigma", param1, "alpha", "beta"),
    "12" = c("z_theta", "unique", "unique_sigma", param1, "alpha", "beta")
  )
  
  param <- irt_mod_params[[as.character(data$irt_mod)]]
  no_fit <- intersect(names(data), c("alpha","beta","gamma","b"))
  
  setdiff(param, no_fit)
}

# NA-safe six-number summary (always returns exactly 6 values)
summarise_values <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(setNames(rep(NA_real_, 6L), .SUMMARY_NAMES))
  }
  q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
  setNames(c(min(x), q[1], q[2], mean(x), q[3], max(x)), .SUMMARY_NAMES)
}

check_params_in_draws <- function(draws, param) {
  base_names <- unique(sub("\\[.*$", "", posterior::variables(draws)))
  missing <- setdiff(param, base_names)
  if (length(missing) > 0L) {
    stop("Parameters not found in fit draws: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
}

# One parameter -> long data frame (one row per diagnostic)
diagnose_param <- function(draws, param) {
  info <- posterior::summarise_draws(
    posterior::subset_draws(draws, variable = param),
    .DIAG_FUNS
  )
  stats <- t(vapply(names(.DIAG_FUNS),
                    function(d) summarise_values(info[[d]]),
                    numeric(6L)))
  data.frame(
    diagnostic = names(.DIAG_FUNS),
    variable   = param,
    stats,
    check.names      = FALSE,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
}


# ---- Main function ---------------------------------------------------------

fit_info_param2 <- function(fit) {
  assert_cmdstanr()
  data <- fit$newdata 
  fit <- fit$model_stan_obj
  
  
  validate_info_inputs(fit, data)
  
  param <- get_info_params(data)
  if (length(param) == 0L) {
    stop("No free parameters to diagnose: all parameters were supplied as data.",
         call. = FALSE)
  }
  
  draws <- fit$draws()
  check_params_in_draws(draws, param)
  
  res <- do.call(rbind, lapply(param, diagnose_param, draws = draws))
  
  # Keep the original ordering: diagnostic first, then parameter
  res <- res[order(match(res$diagnostic, names(.DIAG_FUNS)),
                   match(res$variable, param)), ]
  rownames(res) <- NULL
  res
}
