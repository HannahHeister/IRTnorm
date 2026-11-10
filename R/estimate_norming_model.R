#'@title Convergence check of IRTnorm model parameters
#'@name fit_info_param
#'@description Computes the convergence diagnostics \strong{R-hat},
#' \strong{Effective Sample Size (bulk)} and \strong{Effective Sample Size (tail)}
#' for selected groups of IRTnorm model parameters from a fitted Bayesian IRTnorm model.
#' For each parameter group, the function summarizes the distribution of
#' convergence measures across all indexed parameters (e.g., all estimated
#' beta parameters \code{beta[i]}).
#'@details
#' The function operates on posterior draws extracted via
#' \code{fit$draws()}, and uses the \pkg{posterior} package to compute
#' convergence diagnostics. Parameters are first grouped by name
#' (e.g., all \code{beta[*]} parameters), then diagnostics are computed
#' for each element, and finally summary statistics are taken across elements
#' within each group.
#'
#' The output is ordered by diagnostic, with all parameter groups reported
#' first for \code{rhat}, followed by \code{ess_bulk}, and then \code{ess_tail}.
#'
#' This function is intended for IRTnorm model diagnostics and quality control.
#' Interpretation guidelines commonly used in practice include:
#' \itemize{
#'   \item R-hat values close to 1 (e.g., < 1.01) indicate good mixing,
#'   \item larger ESS values indicate more reliable posterior summaries.
#' }
#'@param fit  the data output from `fit_IRTnorm()`.
#'@param data named list as returned by \code{\link{data_prep}} containing
#'   all data and indexing information required for fitting the IRT norming
#'   model and evaluating its fit.
#'@seealso
#' \code{\link[posterior]{summarise_draws}},
#' \code{\link[posterior]{rhat}},
#' \code{\link[posterior]{ess_bulk}},
#' \code{\link[posterior]{ess_tail}}
#'@return A data.frame summarizing Rhat, Ess bulk and ess tail per variable group.
#'@examples
#' \dontrun{
#' # Fit an IRTnorm model
#' data(response_data)
#' mod_poly <- fit_IRTnorm(data = response_data)
#'
#' # Compute convergence diagnostics
#' converg_tab <- fit_info_param(mod_poly)
#' converg_tab
#' # Inspect R-hat summaries
#' subset(converg_tab, diagnostic == "rhat")
#'
#' # Check for low effective sample sizes
#' subset(converg_tab, diagnostic == "ess_bulk" & Mean < 100)}
#'@export
#'@keywords internal

fit_info_param <- function(fit, data){
  cmdstan_available()
  #fixed irt_mod parameters
  no_fit <- intersect(names(data), c("alpha","beta","gamma","b"))
  
  if(data$age_mod == 1){
    param1 <- c("p1_theta","p2_theta")
    if(data$irt_mod > 8){
      param1 <- c(param1, "p3_theta")
    }
  }
  if(data$age_mod == 2){
    param1 <- c("l1","l2","u1","u2")
    if(data$irt_mod > 8){
      param1 <- c(param1, c("l3","u3"))
    }
  }
  draws <- fit$draws()
  
  irt_mod_params <- list(
    "1"  = c("z_theta", param1, "beta"),
    "2"  = c("z_theta", param1, "alpha", "beta"),
    "3"  = c("z_theta", param1, "alpha", "beta", "gamma"),
    "4"  = c("z_theta", param1, "alpha", "b"),
    "5"  = c("z_theta", "unique", "nu", param1, "beta"),
    "6"  = c("z_theta", "unique", "nu", param1, "alpha", "beta"),
    "7"  = c("z_theta", "unique", "nu", param1, "beta"),
    "8"  = c("z_theta", "unique", "nu", param1, "alpha", "beta"),
    "9"  = c("z_theta", "unique", "nu", "lambda", param1, "beta"),
    "10" = c("z_theta", "unique", "nu", "lambda", param1, "alpha", "beta"),
    "11" = c("z_theta", "unique", "unique_sigma", param1, "beta"),
    "12" = c("z_theta", "unique", "unique_sigma", param1, "alpha", "beta")
  )
  
  param <- irt_mod_params[[as.character(data$irt_mod)]]
  param <- setdiff(param, no_fit)
  
  diag_funs <- list(
    rhat     = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  
  out <- vector("list", length(diag_funs))
  d <- 1
  
  for(diag in names(diag_funs)){
    
    res_param <- vector("list", length(param))
    
    for(i in seq_along(param)){
      
      info <- draws |>
        posterior::subset_draws(variable = param[i]) |>
        posterior::summarise_draws(diag_funs[[diag]])
      
      stats <- as.data.frame(
        t(vapply(info[,-1], summary, numeric(6)))
      )
      
      res_param[[i]] <- cbind(
        diagnostic = rep(diag, nrow(stats)),
        variable   = rep(param[i], nrow(stats)),
        stats
      )
    }
    
    out[[d]] <- do.call(rbind, res_param)
    d <- d + 1
  }
  
  final_info <- do.call(rbind, out)
  rownames(final_info) <- NULL
  return(final_info)
}


#' Estimate an IRT-based norming model
#'
#' Fits a Bayesian item response theory (IRT) norming model that jointly
#' estimates the item parameters and the age-dependent mean and standard
#' deviation of the latent trait. The returned object contains everything
#' needed to evaluate the model and to norm new individuals.
#'
#' @param raw_data A \code{data.frame} or \code{matrix} with individual-level
#'   data, one row per individual. It must contain the age variable and the
#'   item response variables. Item responses may contain missing values
#'   (\code{NA}).
#' @param age_variable Column of \code{raw_data} containing age, given as a
#'   numeric index or a column name. Age is treated as continuous.
#' @param item_variables Columns of \code{raw_data} containing the item
#'   responses, given as a numeric or character vector of indices or names.
#' @param irt_model Character string naming the IRT norming model. One of:
#'   \describe{
#'     \item{\code{"1PLnorm"}}{One-parameter logistic (Rasch) model with
#'       age-dependent mean and variance.}
#'     \item{\code{"2PLnorm"}}{Two-parameter logistic model with age-dependent
#'       mean and variance.}
#'     \item{\code{"3PLnorm"}}{Three-parameter logistic model with
#'       age-dependent mean and variance.}
#'     \item{\code{"GRMnorm"}}{Graded response model with age-dependent mean
#'       and variance.}
#'     \item{\code{"Testlet1PLnorm"}, \code{"Testlet2PLnorm"}}{Higher-order 1PL
#'       or 2PL model with age-dependent mean and variance that are constant
#'       across domains.}
#'     \item{\code{"HO1PLnorm-C"}, \code{"HO2PLnorm-C"}}{Higher-order 1PL or
#'       2PL model with age-dependent mean and variance and a constant
#'       correlation matrix across age.}
#'     \item{\code{"HO1PLnorm-P"}, \code{"HO2PLnorm-P"}}{Higher-order 1PL or
#'       2PL model with age-dependent mean and variance and a correlation
#'       matrix that is proportional across age.}
#'     \item{\code{"HO1PLnorm-F"}, \code{"HO2PLnorm-F"}}{Higher-order 1PL or
#'       2PL model with age-dependent mean and variance and a freely varying
#'       correlation matrix across age.}
#'   }
#' @param age_model Character string specifying how the effect of age on the
#'   latent trait is modeled. One of:
#'   \describe{
#'     \item{\code{"poly"}}{Polynomials of standardized age for the mean and
#'       standard deviation of the latent trait.}
#'     \item{\code{"splines"}}{Penalized B-splines of standardized age for the
#'       mean and standard deviation of the latent trait.}
#'   }
#' @param itemD Only required for multidimensional models. Integer vector of
#'   length equal to the number of items, indicating the domain each item
#'   loads on. Domains must be numbered consecutively starting at 1 (e.g.,
#'   \code{1, 2, 3, 4} for four domains).
#' @param K_item Only required for \code{irt_model = "GRMnorm"}. Integer
#'   vector of length equal to the number of items, giving the number of
#'   response categories of each item (minimum 2, for binary items).
#' @param poly_mean Only required if \code{age_model = "poly"}. Integer giving
#'   the polynomial degree for the age-dependent mean of the latent trait.
#'   For higher-order models with domain-specific trajectories, either a
#'   vector with one degree per domain, or a single integer used for all
#'   domains.
#' @param poly_sd Only required if \code{age_model = "poly"}. Integer giving
#'   the polynomial degree for the age-dependent standard deviation of the
#'   latent trait. Domain-specific degrees are specified as for
#'   \code{poly_mean}.
#' @param poly_sd_unique Only required if \code{age_model = "poly"} and
#'   \code{irt_model} is a higher-order model with a age varing correlation matrix
#'   (\code{-P}) or (\code{-F}).
#'   For HO-P a \code{integer} specifying the degree of the polynomial used to model the
#'   age-dependent lambda of the latent trait.
#'   For HO-F either a \code{vector} indicating the polynomial degree per domain or a \code{integer}
#'   then the same polynomial degree is used for all dimensions.
#' @param k_mean,k_sd Only used if \code{age_model = "splines"}. Number of
#'   knots (basis functions) for the spline of the latent trait mean
#'   (\code{k_mean}, default 10) and standard deviation (\code{k_sd},
#'   default 8).
#' @param k_sd2 Only used if \code{age_model = "splines"} and \code{irt_model}
#'   is a higher-order model with a age-dependent correlation
#'   matrix (\code{-P}) or (\code{-F}).
#'    A \code{num} indicating the number of knots used for modelling the lambda for HO-P 
#'    and the loading of the unique age effect for HO-F. Default is 8.
#' @param fixed_knots_mean,fixed_knots_sd,fixed_knots_sd2 Only used if
#'   \code{age_model = "splines"}. By default, knot positions are chosen from
#'   the data, as in \code{\link[mgcv]{s}}. Use \code{"equidistant"} for
#'   equally spaced knots, or a numeric vector to specify knot positions
#'   directly. They apply to the splines for the mean, the standard
#'   deviation, and the unique age effect (see \code{k_sd2}), respectively.
#' @param prior_knowledge Optional. A named list with one element per
#'   parameter group for which prior knowledge should be incorporated. Each
#'   element is a list of prior means and prior standard deviations. For
#'   parameter \code{b} in the GRM, means and standard deviations must each
#'   be supplied as a list. Default is \code{NULL} (no prior knowledge).
#' @param parameter_fixed Optional. A named list in which each name is a
#'   parameter group to be fixed and each element gives the value(s) to which
#'   it is fixed. Default is \code{NULL} (no fixed parameters).
#' @param iter_warmup Integer. Number of warm-up iterations per chain. Default
#'   is 500.
#' @param iter_sampling Integer. Number of post-warm-up (sampling) iterations
#'   per chain. Default is 500.
#' @param chains Integer. Number of Markov chains. Default is 4.
#' @param parallel_chains Optional integer. Number of chains run in parallel.
#'   If \code{NULL} (default), as many chains as \code{chains} are run in
#'   parallel, limited by the number of available cores
#' @param seed Optional integer random seed for reproducibility. If
#'   \code{NULL} (default), a seed is drawn at random.
#' @param extract_stan_object Logical. If \code{TRUE} (default), the raw Stan
#'   fit object is included in the result, giving direct access to posterior
#'   draws and sampler diagnostics. Set to \code{FALSE} to reduce memory use.
#'
#' @return An object of class \code{"IRTnorm"}, a list with components:
#' \describe{
#'   \item{\code{model_info}}{Model-level summary, including parameter
#'     estimates, convergence diagnostics, and model fit indices.}
#'   \item{\code{norm_info}}{Information required to norm new individuals.}
#'   \item{\code{conv_info}}{Table summarizing the effective sample size and Rhat for
#'     each parameter group.}
#'   \item{\code{newdata}}{The prepared data used for fitting, including the
#'     design matrices.}
#'   \item{\code{raw_data}}{The data passed to \code{raw_data}.}
#'   \item{\code{age_variable}}{The age variable, as supplied.}
#'   \item{\code{item_variables}}{The item response variables, as supplied.}
#'   \item{\code{model_stan_obj}}{The raw Stan fit object. Only present if
#'     \code{extract_stan_object = TRUE}.}
#' }
#' @details
#' The function prepares the data, fits the model with a Bayesian sampler,
#' and extracts the information needed to assess convergence, summarize the
#' model, and norm new individuals. To choose the best model specification, see
#' \code{\link{compare_age_models}} and \code{\link{compare_IRTnorm_models}}.

#' @seealso \code{\link{compare_age_models}},
#'   \code{\link{compare_IRTnorm_models}}
#'
#' @examples
#' \dontrun{
#' data(my_data)
#' 
#' # 1PL norming model with a spline age trajectory
#' result <- estimate_norming_model(
#'   raw_data = my_data,
#'   age_variable = "age",
#'   item_variables = 1:25,
#'   irt_model = "1PLnorm",
#'   age_model = "splines" )
#'
#' # Model summary
#' result$model_info
#'
#'# 1PL norming model with prior knowledge on the beta values 
#' prior_info = list(beta = list(beta_prior_mu = c(-0.7, -0.4, -0.6, -0.8, -0.5, -2.4,
#'                                                 -2.0, -1.0, -2.2, 2.4, -2.5, 2.9, 
#'                                                 -0.3, -2.1, 2.8, 1.2, 0.7, -1.0, 2.7,
#'                                                 -1.7, -2.0, -1.0, -2.5, -1.6, -1.6),
#'                                                 beta_prior_sd = rep(0.3, 25)))
#'
#'result2 <- estimate_norming_model(raw_data = my_data,
#'                                 age_variable = "age",
#'                                 item_variables = 1:25,
#'                                 irt_model = "1PLnorm",
#'                                  age_model = "splines", 
#'                                  prior_knowledge = prior_info)
#' # Model summary
#' result2$model_info
#' }
#' @export
estimate_norming_model <- function(raw_data, age_variable, item_variables, 
                                   irt_model,
                                   age_model,
                                   itemD = NULL, K_item = NULL,
                                   poly_mean = NULL, poly_sd = NULL, poly_sd_unique = NULL,
                                   k_mean = 10, k_sd = 8,k_sd2 = 8, fixed_knots_mean = NULL, fixed_knots_sd = NULL,fixed_knots_sd2 = NULL,
                                   prior_knowledge = NULL, parameter_fixed = NULL, 
                                   iter_warmup = 500, iter_sampling = 500,
                                   chains = 4, parallel_chains = NULL, seed = NULL,
                                   extract_stan_object = TRUE){
  
  if(!is.logical(extract_stan_object)){
    stop("extract_stan_object has to be logical, e.g either TRUE or FALSE.")
  }
  check_rawinput(raw_data = raw_data, item_variables = item_variables, 
                 age_variable = age_variable)
# data preperation 
  newdata <- data_prep(raw_data = raw_data,
                       age_variable = age_variable,
                       item_variables = item_variables, 
                       irt_model = irt_model, 
                       age_model = age_model,
                       itemD = itemD, K_item = K_item,
                       poly_mean = poly_mean, poly_sd = poly_sd, 
                       poly_sd_unique = poly_sd_unique,
                       k_mean = k_mean, k_sd = k_sd,
                       k_sd2 = k_sd2,
                       fixed_knots_mean = fixed_knots_mean,
                       fixed_knots_sd = fixed_knots_sd,
                       fixed_knots_sd2 = fixed_knots_sd2,
                      prior_knowledge = prior_knowledge,
                      parameter_fixed = parameter_fixed)
  
  # fit IRTnorm model 
  fit_model <- fit_IRTnorm(data = newdata,
                          seed = seed, iter_warmup = iter_warmup,
                          iter_sampling = iter_sampling, chains = chains, 
                          parallel_chains = parallel_chains)
    
  # extract information from the norming model 
  model_info <- extract_info(fit =  fit_model, data = newdata)
  
  # extract information for norming new individual
  norm_info <- info_new_norming(fit =  fit_model, data = newdata)
  
  
  # convergence information 
  conv_info <-  fit_info_param(fit = fit_model, data = newdata)
  
  
  # combin all information 
  all_info <- list(model_info = model_info,
                   norm_info = norm_info,
                   conv_info = conv_info, 
                   newdata = newdata,
                   raw_data = raw_data, 
                   age_variable = age_variable,
                   item_variables = item_variables)
  
    
  if(extract_stan_object == TRUE){
    all_info[["model_stan_obj"]] <- fit_model
  }
  class(all_info) <- "IRTnorm"
    return(all_info)
}


  
