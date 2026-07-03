#'@title Fit an IRT-based continuous norming model using Stan
#'@name fit_IRTnorm
#'@description
#' 

#' The \code{fit_IRTnorm} function fits an IRT-based continuous norming model to
#' the data prepared with \code{\link{data_prep}} using Bayesian estimation via \pkg{cmdstanr}.
#' The specific model is specified with data_prep().
#' Depending on the specified model in the prepared data object, the
#' The specifics of the IRT-based continuous norming model are specified in
#' \code{data} using \code{\link{data_prep}}.
#' The function supports fitting the norming model with either polynomial or spline-based
#' age effects.
#'
#'
#' @param data a named list as returned by \code{\link{data_prep}} containing
#'   all data and indexing information required for fitting the IRT norming
#'   model and evaluating its fit.
#'
#' @param seed an optional \code{integer} specifying the random seed for
#'   reproducibility. If \code{NULL} (default), a random seed is generated
#'   internally.
#'
#' @param iter_warmup an \code{integer} specifying the number of warm-up
#'   iterations per chain used by the MCMC sampler.
#'
#' @param iter_sampling an \code{integer} specifying the number of post
#'   warm-up (sampling) iterations per chain.
#'
#' @param chains an \code{integer} specifying the number of Markov chains to
#'   run in parallel. The default is 4 chains.
#'
#' @param parallel_chains an optional \code{integer} specifying the number of
#'   CPU cores used for parallel sampling. If \code{NULL} (default), all
#'   available cores detected by \code{\link[parallel]{detectCores}} are used.
#' @param ... Additional arguments passed directly to the
#'   \code{sample()} method of the \pkg{cmdstanr} model object (e.g.,
#'   \code{adapt_delta}, \code{max_treedepth}).
#'
#' @return A \pkg{cmdstanr} \code{CmdStanMCMC} object containing posterior
#'   draws, sampler diagnostics, and metadata for the fitted IRT norming
#'   model.
#' @seealso [extract_info()],[data_prep()]
#'
#' @examples
#' \dontrun{
#' ## Simulate example data
#' set.seed(1234)
#' ## load data
#' data(response_data)
#'
#
#' ## Prepare data for the 2PL norming model
#' prep_data <- data_prep(
#'   raw_data = response_data,
#'   age_variable = "age",
#'   int_variables = 1:50,
#'   irt_mod = "2PLnorm",
#'   age_mod = "polynom",
#'   poly_mean = 2,
#'   poly_sd = 1
#' )
#'
#' ## Fit the IRT norming model
#' fit <- fit_IRTnorm(
#'   data = prep_data,
#'   iter_warmup = 250,
#'   iter_sampling = 250,
#'   chains = 2
#' )
#'}
#'@export
#'@keywords internal
fit_IRTnorm <- function(data, seed = NULL,
                        iter_warmup = 500, iter_sampling = 500,
                        chains = 4, parallel_chains = NULL, ...){
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package 'cmdstanr' is required but not installed.\n",
      "Install it with: install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))",
      call. = FALSE
    )
  }
  cmdstanr::check_cmdstan_toolchain()
  age_model <- c("poly","splines")
  IRTmodel <- c("1PLnorm","2PLnorm","3PLnorm", "GRMnorm",
                "Testlet1PLnorm","Testlet2PLnorm",
                "HO1PLnorm-C","HO2PLnorm-C",
                "HO1PLnorm-P","HO2PLnorm-P",
                "HO1PLnorm-F","HO2PLnorm-F")

  if(data$age_mod == 2){
    data$sm1 <- NULL
    data$sm2 <- NULL
    data$sm3 <- NULL
    
  }
  # any parameter_fixed
  if(any(duplicated(c(names(data),"alpha","beta","gamma", "b")))){
    param_fixed <-intersect(names(data), c("alpha","beta","gamma", "b"))
  }else{
    param_fixed <- NULL
  }
  if(any(duplicated(c(names(data), "alpha_prior_mu","alpha_prior_sd",
                      "beta_prior_mu","beta_prior_sd",
                      "gamma_prior_mu","gamma_prior_sd",
                      "b_prior_mu","b_prior_sd",
                      "nu_prior_mu","nu_prior_sd",
                      "lambda_prior_mu","lambda_prior_sd",
                      "unique_sigma_prior_mu", "unique_sigma_prior_sd")))){
    prior_known_temp <-intersect(names(data), c("alpha_prior_mu","alpha_prior_sd",
                                           "beta_prior_mu","beta_prior_sd",
                                           "gamma_prior_mu","gamma_prior_sd",
                                           "b_prior_mu","b_prior_sd",
                                           "nu_prior_mu","nu_prior_sd",
                                           "lambda_prior_mu","lambda_prior_sd",
                                           "unique_sigma_prior_mu", "unique_sigma_prior_sd"))
    prior_known <- NULL
    if(any(substr(prior_known_temp,1,5) == "alpha")){
      prior_known <- c(prior_known,"alpha")
    }
    if(any(substr(prior_known_temp,1,4) == "beta")){
      prior_known <- c(prior_known,"beta")
    }
    if(any(substr(prior_known_temp,1,5) == "gamma")){
      prior_known <- c(prior_known,"gamma")
    }
    if(any(substr(prior_known_temp,1,2) == "b_")){
      prior_known <- c(prior_known,"b")
    }
    if(any(substr(prior_known_temp,1,2) == "nu")){
      prior_known <- c(prior_known,"nu")
    }

    }else{
    prior_known <- NULL
  }

  mod <- get_compiled_model(age_model = age_model[data$age_mod],
                            irt_model = IRTmodel[data$irt_mod],
                            prior_knowledge = prior_known,
                            parameter_fixed  = param_fixed)

  if(length(parallel_chains) == 0){
    parallel_chains <- parallel::detectCores()
  }

  if(is.null(seed)){
    seed <- sample.int(.Machine$integer.max, 1)
  }

  # Capture extra arguments
  #  dots <- list(...)
  if(data$irt_mod != 4){
    fit <- mod$sample(
      data = data,
      seed = seed,
      chains = chains,
      parallel_chains = parallel_chains,
      refresh = (iter_sampling + iter_warmup)/10 ,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling, ...)
  }
  if(data$irt_mod == 4){
    fit <- mod$sample(
      data = data,
      seed = seed,
      chains = chains,
      parallel_chains = parallel_chains,
      refresh = (iter_sampling + iter_warmup)/10 ,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,init = 1, ...)
  }
  # Base arguments shared by both cases
  #base_args <- list(
  #  data             = data,
  #  seed             = seed,
  #  chains           = chains,
  #  parallel_chains  = parallel_chains,
  #  refresh          = (iter_sampling + iter_warmup) / 10,
  # iter_warmup      = iter_warmup,
  #  iter_sampling    = iter_sampling
  #)

  #if(data$irt_mod == 8){
  # Use init = 1 as default, but allow user to override via ...
  #  if(!"init" %in% names(dots)){
  #  base_args$init <- 1
  # }
  #}

  # Merge base args with dots — dots override base_args if names clash
  #final_args <- c(base_args, dots)

  #fit <- do.call(mod$sample, final_args)

  return(fit)
}




# Helper functions for extract info 
extract_iteminfo <- function(fit, data,
                             type = c("mean", "draws")) {
  
  type <- match.arg(type)
  nitem <- data$I
  ndraws <-  fit$metadata()$iter_sampling * fit$num_chains()
  # if item information was fixed use that information
  for (var in c("alpha", "beta", "gamma", "b")) {
    if (var %in% names(data)) {
      if(type == "mean"){
        assign(var, data[[var]])
      }
      if(type == "draws"){
        assign(var, rep(data[[var]], each = ndraws))
      }
    }
  }
  
  # beta
  if (!("beta" %in% names(data)) && data$irt_mod != 4) {
    beta <- get_value(fit = fit, var = "beta", type = type)
  }
  
  # alpha
  if (!("alpha" %in% names(data))) {
    if (!(data$irt_mod %in% c(1, 5, 7, 9, 11))) {
      alpha <- get_value(fit = fit, var = "alpha", type = type)
    } else {
      if(type == "mean"){
        alpha <- 1
      }
      if(type == "draws"){
        alpha <- matrix(1,ndraws, nitem)
      }
    }
  }
  
  # gamma
  if (!("gamma" %in% names(data))) {
    if (data$irt_mod == 3) {
      gamma <- get_value(fit = fit, var = "gamma", type = type)
    } else {
      if(type == "mean"){
        gamma <- 0
      }
      if(type == "draws"){
        gamma <- matrix(0,ndraws, nitem)
      }
    }
  }
  
  # thresholds
  if (!("b" %in% names(data)) && data$irt_mod == 4) {
    b <- get_value(fit = fit, var = "b", type = type)
  }
  
  if (data$irt_mod != 4) {
    if (type == "mean") {
      info_item <- data.frame(
        alpha = alpha,
        beta  = beta,
        gamma = gamma)
    } else {
      info_item <- list(
        alpha = alpha,
        beta  = beta,
        gamma = gamma)
    }
  } else {
    info_item <- list(
      alpha     = alpha,
      threshold = b
    )
  }
  
  return(info_item)
}

get_value <- function(fit, var, type) {
  if (type == "mean") {
    fit$summary(variables = var, "mean")$mean
  } else {
    fit$draws(variables = var, format = "matrix")
  }
}

extract_coefinfo <- function(fit, age_mod, irt_mod, type = c("mean", "draws")) {
  
  type <- match.arg(type)
  
  vars <- switch(
    as.character(age_mod),
    "1" = c("p1_theta", "p2_theta", if (irt_mod > 8) "p3_theta"),
    "2" = c("l1", "u1", "l2", "u2", if (irt_mod > 8) c("l3", "u3"))
  )
  
  stats::setNames(lapply(vars, function(x){get_value(fit = fit, type = type, var = x)}), vars)
}
extract_personinfo <- function(fit, data){
  if(data$irt_mod <= 4){
    
    info_sample <- data.frame(age = data$age_variab)
    
    normscore <- fit$summary(variables = c("z_theta"))
    
    info_sample$normscore_est <- normscore$mean
    info_sample$normscore_est_upper <- normscore$q95
    info_sample$normscore_est_lower <- normscore$q5
    info_sample$theta_est <- as.data.frame(fit$summary(variables = c("theta"), "mean"))$mean
    info_sample$theta_mu <- as.data.frame(fit$summary(variables = c("theta_mean"), "mean"))$mean
    info_sample$theta_sigma <- as.data.frame(fit$summary(variables = c("theta_sigma"), "mean"))$mean
    info_sample$irt_mod <- data$irt_mod
  }
  # HO and Testlet model
  if(data$irt_mod > 4){
    nitem <- data$I
    nperson <- data$J
    D <- data$D
    itemD <- data$itemD
    theta_est <- fit$summary(variables = "z_theta")
    delta_est_all <- fit$draws("substraits",format =  "draws_matrix")
    
    int_delta <- which(data$ii %in% which(!duplicated(itemD)))
    delta_est <-  delta_est_all[,int_delta]
    
    mu_est <- fit$draws( "theta_mean",format =  "draws_matrix")
    sd_est <- fit$draws(  "theta_sigma",format =  "draws_matrix")
    
    if(data$irt_mod == 5 |data$irt_mod == 6){
      mu_est <-   do.call(cbind, replicate(D, mu_est, simplify = FALSE))
      sd_est <-   do.call(cbind, replicate(D, sd_est, simplify = FALSE))
      
    }
    
    delta_normed <- (delta_est - mu_est)/(sd_est)
    
    q_delta <- apply(delta_est, 2, stats::quantile, probs = c(0.05, 0.95))
    q_deltanorm <- apply(delta_normed, 2, stats::quantile, probs = c(0.05, 0.95))
    
    m_mu_est <- colMeans(mu_est)
    m_sd_est <- colMeans(sd_est)
    
    
    info_sample <- data.frame(est_ability = c(colMeans(delta_est),
                                              colMeans(delta_normed),
                                              theta_est$mean),
                              ability_est_upper = c(q_delta[2,], q_deltanorm[2,], theta_est$q95),
                              ability_est_lower = c(q_delta[1,], q_deltanorm[1,], theta_est$q5),
                              trait = c(rep(paste0("delta",1:D), each = nperson),
                                        rep(paste0("normeddelta",1:D), each = nperson),
                                        rep(c("theta"),nperson)),
                              theta_mu = c(m_mu_est, rep(NA, nperson*(D+1))),
                              theta_sigma = c(m_sd_est, rep(NA, nperson*(D +1))),
                              age = rep(data$age_variab, D*2 + 1),
                              irt_mod = data$irt_mod)
    
  }
  return(info_sample)
}

extract_relationinfo <- function(fit, irt_mod,
                                 type = c("mean", "draws")) {
  
  type <- match.arg(type)
  
 #if (irt_mod <= 4) {
  #  return(NULL)
  #}
  
  relation <- list(
    nu = NULL,
    lambda = NULL,
    unique_sigma = NULL
  )
  
  if (irt_mod %in% c(5, 6, 7, 8, 9, 10)) {
    relation[["nu"]] <- get_value(fit = fit, var = "nu", type = type)
  }
  
  if (irt_mod %in% c(9, 10)) {
    relation[["lambda"]] <- get_value(fit = fit, "lambda", type = type)
  }
  
  if (irt_mod %in% c(11, 12)) {
    relation[["unique_sigma"]] <- get_value(fit = fit, "unique_sigma", type = type)
  }
  
  return(relation)
}
#'@title Extract person-, item-, and norming-level information from a
#'fitted IRT-based continuous norming model.
#'@name extract_info
#'@description
#' The \code{extract_info} function extracts key posterior summaries and
#' posterior draws from a fitted Stan-based IRT norming model obtained with
#' \code{\link{fit_IRTnorm}}. The function returns person-level norm scores and
#' latent trait estimates, item parameter estimates, full posterior draws, and
#' replicated response data for posterior predictive checks.
#' @param fit A \pkg{cmdstanr} \code{CmdStanMCMC} object as returned by
#'   \code{\link{fit_IRTnorm}}.
#'
#' @param data A named list as returned by \code{\link{data_prep}} containing
#'   the original data and indexing information used for model fitting.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{info_sample}{A list containing the original prepared data augmented
#'     with person-level posterior summaries, including estimated norm scores,
#'     latent trait estimates, and age-dependent mean and standard deviation of
#'     the latent trait.}
#'
#'   \item{info_item}{A \code{data.frame} containing posterior mean estimates of
#'     item parameters.}
#'
#'   \item{y_rep}{A matrix containing posterior predictive replicated item
#'     responses (\code{y_rep}), useful for posterior predictive checks.}
#'
#'   \item{info_coef}{A \code{list}) containing the posterior mean estimates of
#'   the regression parameters for the age effects.}
#'
#'   \item{nu_est}{For multidimensinal models a \code{vector} containing posterior mean estimates of
#'     item parameters, and for all models holding the value 1.}
#' }
#'@details
#' Person-level norm scores are extracted from the posterior distribution of
#' \code{z_theta}. The function reports posterior means as point estimates and
#' the 5% and 95% quantiles as uncertainty bounds.
#'
#' Item parameters are summarized using posterior means only.
#'
#' The returned object is intended to be used for further norming analyses,
#' reporting, and model diagnostics.
#'
#'@seealso \code{\link{data_prep}}, \code{\link{fit_IRTnorm}},  \code{\link{norm_new_individual}}
#'
#' @examples
#' \dontrun{
#' # Step 1: load data
#' data(response_data)
#'
#' # Step 2: Prepare data
#' prep_data <- data_prep(
#'   raw_data = response_data,
#'   age_variable = "age",
#'   int_variables = 1:50,
#'   model = "2PLnorm",
#'
#' # Step 3: Fit the IRT norming model
#' mod_sp <- fit_IRTnorm(
#'   data = prep_data,
#'   iter_warmup = 500,
#'   iter_sampling = 500,
#'   chains = 4)
#'
#' # Step 4: Extract information
#' info <- extract_info(fit = mod_sp, data = prep_data)
#' }
#'@export
#'@keywords internal

extract_info <- function(fit, data){
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package 'cmdstanr' is required but not installed.\n",
      "Install it with: install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))",
      call. = FALSE
    )
  }
  info_item <- extract_iteminfo(fit = fit, data = data, type = "mean")

  info_coef <- extract_coefinfo(fit = fit, age_mod = data$age_mod, irt_mod = data$irt_mod,  type = "mean")

  info_sample <- extract_personinfo(fit = fit, data = data) 
  
  info_relation <- extract_relationinfo(fit = fit, irt_mod  = data$irt_mod, type = "mean")
  
  y_rep <- fit$draws("y_rep", format = "matrix")

  info <- list(info_sample = info_sample,
               info_item = info_item,
               y_rep = y_rep,
               info_coef = info_coef,
               info_relation = info_relation)
  return(info)
}

#'#'@title Extract all posterior draws of the norming model to enable norming 
#' of new individuals. 
#'@name info_new_norming
#'@description
#' The \code{info_new_norming} function extracts the posterior draws of all 
#' item parameter and regression coeficents from a fitted Stan-based IRT norming
#' model obtained with \code{\link{fit_IRTnorm}}. 
#' @param fit A \pkg{cmdstanr} \code{CmdStanMCMC} object as returned by
#'   \code{\link{fit_IRTnorm}}.
#'
#' @param data A named list as returned by \code{\link{data_prep}} containing
#'   the original data and indexing information used for model fitting.
#'
#'@return a list containing posterior draws relevant for norming new individuals. 
#'
#'#'@seealso \code{\link{data_prep}}, \code{\link{fit_IRTnorm}}
#'
#' @examples
#' \dontrun{
#' # Step 1: load data
#' data(response_data)
#'
#' # Step 2: Prepare data
#' prep_data <- data_prep(
#'   raw_data = response_data,
#'   age_variable = "age",
#'   int_variables = 1:50,
#'   model = "2PLnorm",
#'   age_model = "splines")
#'#' # Step 3: Fit the IRT norming model
#' mod_sp <- fit_IRTnorm(
#'   data = prep_data,
#'   iter_warmup = 500,
#'   iter_sampling = 500,
#'   chains = 4)
#'
#' # Step 4: Extract information
#' info <- info_new_norming(fit = mod_sp, data = prep_data)
#' }
#'@export
#'@keywords internal
#'
#'

info_new_norming  <- function(fit, data){
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package 'cmdstanr' is required but not installed.\n",
      "Install it with: install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))",
      call. = FALSE
    )
  }
  info_coef <- extract_coefinfo(fit = fit, age_mod = data$age_mod, 
                                irt_mod = data$irt_mod,  type = "draws")

  info_item <- extract_iteminfo(fit = fit, data = data, type = "draws")
  
  info_relation <- extract_relationinfo(fit = fit, irt_mod  = data$irt_mod, type = "draws")
  return(c(info_coef,info_item,info_relation))
}
