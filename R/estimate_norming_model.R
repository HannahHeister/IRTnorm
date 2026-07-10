
#' Estimate a IRT-based norming model
#' @name estimate_norming_model
#' @description
#' Fits a Bayesian Item Response Theory (IRT) norming model that jointly estimates
#' item parameters and age-based regression models for the mean and standard
#' deviation of the latent trait. The function prepares the data, fits the Stan
#' model, and extracts all information needed for model evaluation, and 
#' subsequently norm of new individuals.
#'
#'@param raw_data a \code{data.frame} or \code{matrix} containing individual-level
#'   data. Each row represents one individual. The object must contain at least
#'   an age variable and item response variables. Item responses may
#'   contain missing values (\code{NA}).
#'
#'@param age_variable either a \code{numeric} index or a \code{character} string specifying
#'   the column in \code{raw_data} that contains the individuals' age information.
#'   Age is assumed to be continuous.
#'
#'@param item_variables a \code{numeric} or \code{character} vector specifying the columns in
#'   \code{raw_data} that contain the item response variables.
#'
#'@param irt_model a \code{character} scalar indicating for which IRT norming model the
#'   data should be prepared. The supported options are:
#'   \describe{
#'     \item{\code{"1PLnorm"}}{One-parameter logistic IRT model / Rasch model with
#'     age-dependent mean and variance.}
#'      \item{\code{"2PLnorm"}}{Two-parameter logistic IRT model with
#'     age-dependent mean and variance.}
#'      \item{\code{"3PLnorm"}}{Three-parameter logistic IRT model with
#'     age-dependent mean and variance.}
#'      \item{\code{"GRMnorm"}}{Graded Response model with
#'     age-dependent mean and variance.}
#'      \item{\code{"Testlet1PLnorm"}}{Higher order 1PL IRT model with
#'     age-dependent mean and variance constant across domains.}
#'      \item{\code{"Testlet2PLnorm"}}{Higher order 2PL IRT model with
#'     age-dependent mean and variance constant across domains.}
#'      \item{\code{"HO1PLnorm-C"}}{Higher order 1PL IRT model with
#'     age-dependent mean and variance and constant correlation matrix.}
#'      \item{\code{"HO2PLnorm-C"}}{Higher order 2PL IRT model with
#'     age-dependent mean and variance and constant correlation matrix.}
#'    \item{\code{"HO1PLnorm-P"}}{Higher order 1PL IRT model with
#'     age-dependent mean and variance and proportional correlation matrix across age.}
#'      \item{\code{"HO2PLnorm-P"}}{Higher order 2PL IRT model with
#'     age-dependent mean and variance and proportional correlation matrix across age.}
#'     \item{\code{"HO1PLnorm-F"}}{Higher order 1PL IRT model with
#'     age-dependent mean and variance and free correlation matrix across age.}
#'      \item{\code{"HO2PLnorm-F"}}{Higher order 2PL IRT model with
#'     age-dependent mean and variance and free correlation matrix across age.}
#'   }
#'
#' @param age_model a \code{character} scalar indicating how the effect of age on the
#'   latent trait should be modeled. Available options are:
#'   \describe{
#'     \item{\code{"polynom"}}{Polynomial functions of standardized age are used
#'     to model the mean and standard deviation of the latent trait.}
#'     \item{\code{"spline"}}{Penalized B-spline basis functions of standardized age are
#'     used to model the mean and standard deviation of the latent trait.}
#'   }

#' @param poly_mean only required if \code{age_model = "polynom"}.
#'   For all models except HO2PLnorm an \code{integer} specifying the degree of the polynomial used to model the
#'   age-dependent mean of the latent trait. For HO2PLnorm
#'   either a \code{vector} indicating the polynomial degree per domain or a \code{integer}
#'   then the same polynomial degree is used for all dimensions.
#'
#' @param poly_sd only required if \code{age_model = "polynom"}.
#'   For all models except HO2PLnorm an \code{integer} specifying the degree of the polynomial used to model the
#'   age-dependent standard deviation of the latent trait. For HO2PLnorm
#'   either a \code{vector} indicating the polynomial degree per domain or a \code{integer}
#'   then the same polynomial degree is used for all dimensions.
#'
#' @param poly_sd_unique only required if \code{age_model = "polynom"} and \code{irt_model} a HO-P or HO-C model.
#'  For HO-P a\code{integer} specifying the degree of the polynomial used to model the
#'   age-dependent lambda of the latent trait. For HO-F
#'   either a \code{vector} indicating the polynomial degree per domain or a \code{integer}
#'   then the same polynomial degree is used for all dimensions.
#' 
#'@param itemD only required if \code{model} is a multidimensional model.
#' A \code{vector} indicating on which domain each item is loading on.
#' It is important that the first domain has the value 1 and the following domains
#' have the subsequent integers. E.g., for 4 domains the vector should only contain the
#' values 1,2,3,4. itemD has the length of the number of tested items.
#'
#' @param K_item only required if \code{model = "GRMnorm"}.
#' A \code{vector} indicating the number of categories of each item.
#' The minimum number in K_item is 2, since binary items have two response categories.
#' K_item has the length of the number of tested items.
#'
#' @param k_mean only required if \code{age_model = "splines"}.
#' A \code{num} indicating the number of knots used for modelling the
#' latent trait mean.
#'
#' @param k_sd only required if \code{age_model = "splines"}.
#' A \code{num} indicating the number of knots used for modelling the
#' latent trait standard deviation
#' 
#' @param k_sd2 only required if \code{age_model = "splines"}  and \code{irt_model} a HO-P or HO-C model.
#' A \code{num} indicating the number of knots used for modelling the lambda for HO-P 
#' and the loading of the unique age effect for HO-F. 
#'
#'@param fixed_knots_mean only required if knots should not be chosen based on
#' data such as the default in \code{s} from \code{mgcv}. It can be choose between
#' equidistant knots by setting fixed_knots_sd =  "equidistant" or fully user
#' specified which a numerical vector of knot values
#'
#'@param fixed_knots_sd only required if knots should not be chosen based on
#' data such as the default in \code{s} from \code{mgcv}. It can be choose between
#' equidistant knots by setting fixed_knots_sd =  "equidistant" or fully user
#' specified which a numerical vector of knot values.
#'
#'@param fixed_knots_sd2 only required if knots should not be chosen based on
#' data such as the default in \code{s} from \code{mgcv}. It can be choose between
#' equidistant knots by setting fixed_knots_sd =  "equidistant" or fully user
#' specified which a numerical vector of knot values.
#'
#' @param prior_knowledge only required if prior knowledge should be incoperated
#' in the model. A list of the parameter groups for which prior_knowledge
#' should be incorporated for each parameter group with prior knowledge a list of
#' prior means and prior standard deviations has to be indicated. For b in the GRM model
#' this has to be a list per means and standard deviations.
#' Below in the example a example list is used.
#'
#' @param parameter_fixed only required if model parameters should be fixed.
#' A list where the parameter groups that should be fixed are indicated with the list names
#' and the list elements indicating the values to which the parameters should be fixed.
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
#'   
#' @param extract_stan_object A logical indicating whether to include the raw
#'   Stan/CmdStanR fit object in the returned list. Default is \code{TRUE}, such that
#'   posterior draws, diagnostics, and other Stan-level output can be access directly.
#'   Set to \code{FALSE} if working memory is an issue.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{\code{model_info}}{A list of summarized model-level information
#'       extracted by \code{extract_info()}, including parameter estimates,
#'       convergence diagnostics, and model fit indices.}
#'     \item{\code{norm_info}}{A list of information required to norm new
#'       individuals, as extracted by \code{info_new_norming()}.}
#'        \item{\code{conv_info}}{Table with a summary of the effective sample 
#'        size and Rhat for every parameter group.}
#'     \item{\code{newdata}}{The prepared data object returned by
#'       \code{data_prep()}, containing the design matrices and Stan input
#'       list.}
#'     \item{\code{raw_data}}{The original data frame passed to
#'       \code{raw_data}.}
#'     \item{\code{age_variable}}{The name of the age variable, as supplied.}
#'     \item{\code{item_variables}}{The interaction variable names, as supplied.}
#'     \item{\code{model_stan_obj}}{The raw Stan fit object. Only present when
#'       \code{extract_stan_object = TRUE}.}
#'   }
#'
#' @details
#' The function implements a three-stage workflow:
#' \enumerate{
#'   \item \strong{Data preparation}: \code{data_prep()} constructs the Stan
#'     input list, including the IRT design matrix, spline or polynomial basis
#'     matrices for age, and any interaction terms.
#'   \item \strong{Model fitting}: The Stan model is called via
#'     \code{fit_IRTnorm()}, which runs a Bayesian HMC sampler with the
#'     specified number of warmup and sampling iterations.
#'   \item \strong{Post-processing}: \code{fit_info_param} extracts information 
#'   for convergence checks, \code{extract_info()} summarizes the
#'     posterior and \code{info_new_norming()} computes the quantities needed
#'     to derive norm scores for new respondents.
#' }
#'
#' @seealso
#' \code{\link{data_prep}}, \code{\link{extract_info}},
#' \code{\link{info_new_norming}}
#'
#' @examples
#' \dontrun{
#' # Minimal example with a 2PLnorm model 
#' result <- estimate_norming_model(
#'   raw_data      = my_data,
#'   age_variable  = "age",
#'   item_variables = 1:50,
#'   irt_model     = "2PLnorm",
#'   age_model     = "splines",
#'   iter_warmup   = 500,
#'   iter_sampling = 1000,
#'   chains        = 4,
#'   parallel_chains = 4
#' )
#'
#' # Access model summary
#' result$model_info
#'}
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
  return(all_info)
}


  
