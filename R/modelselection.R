loo_extr <- function(data, fit){
  # info for loo comparison 
  log_lik <- fit$draws("log_lik", format = "matrix")
  loo_loglike <- loo::loo(log_lik)
  return(loo_loglike)
}

loo_comp <- function(loo_info){
  best <-  max(loo_info[[1]]$ELPD)
  reference <-  which(loo_info[[1]]$ELPD == max(loo_info[[1]]$ELPD))
  
  best_list <- order(loo_info[[1]]$ELPD,decreasing = TRUE)
  loo_order <- loo_info[[1]][best_list,]
  loo_order$ELPD_diff <- loo_order$ELPD - best 
  
  d <- loo_info[[2]][,reference[1]] - loo_info[[2]][,best_list]
  loo_order$SE <- sqrt(nrow(d) * apply(d, 2, stats::var))
  loo_order$rank <- seq_len(nrow(loo_order))
  n_ref <- length(reference)
  if(n_ref > 1){
    loo_order$rank <- c(rep(1,n_ref), 
                        2:(nrow(loo_order) - n_ref +1))
  }
  return(loo_order)
}

fit_one_model <- function(model_specifications, raw_data, age_variable, item_variables,
                          iter_warmup = 500, iter_sampling = 500, traj_info = FALSE, 
                          chains = 4, parallel_chains = NULL, seed = NULL) {
  
  defaults <- list(
    irt_model        = NULL,
    age_model        = NULL,
    itemD            = NULL,
    K_item           = NULL,
    poly_mean        = NULL,
    poly_sd          = NULL,
    poly_sd_unique   = NULL,
    k_mean           = 10,
    k_sd             = 8,
    k_sd2            = 8,
    fixed_knots_mean = NULL,
    fixed_knots_sd   = NULL,
    fixed_knots_sd2  = NULL,
    prior_knowledge  = NULL,
    parameter_fixed  = NULL
  )
  model_specifications <- utils::modifyList(defaults, model_specifications)
  
  newdata <- do.call(data_prep, c(
    list(raw_data = raw_data, age_variable = age_variable, item_variables = item_variables),
    model_specifications
  ))
  
  fit <- fit_IRTnorm(
    data            = newdata,
    seed            = seed,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    chains          = chains,
    parallel_chains = parallel_chains
  )
  
  loo_info <- loo_extr(data = newdata, fit = fit)
  
  # generic model label, works for polynomial grid (mu_sd) or any other spec
  modelcombi <- if (!is.null(model_specifications$poly_mean) || !is.null(model_specifications$poly_sd)) {
    paste0(model_specifications$poly_mean, "_", model_specifications$poly_sd)
  } else {
    paste0(model_specifications$irt_model, "_", model_specifications$age_model)
  }
  if(is.null(model_specifications$poly_mean)){
    model_specifications$poly_mean <- NA 
  }
  if(is.null(model_specifications$poly_sd)){
    model_specifications$poly_sd <- NA 
  }
  return_info <-   list(
    elpd_pointwise = loo_info$pointwise[, "elpd_loo"],
    elpd_summary   = data.frame(
      ELPD       = loo_info$estimates[[1]],
      irt_model  = model_specifications$irt_model,
      age_model  = model_specifications$age_model,
      poly_mean  = model_specifications$poly_mean,
      poly_sd    = model_specifications$poly_sd,
      modelcombi = modelcombi
    ), 
    paretok =  loo_info$diagnostics$pareto_k, 
    neff =  loo_info$diagnostics$n_eff,
    reff =  loo_info$diagnostics$r_eff)
  if(traj_info){
    age_mean <- fit$summary(variables = "theta_mean")$mean
    age_sd   <- fit$summary(variables = "theta_sigma")$mean
    age      <- newdata$age_variable
    nperson  <- newdata$J
    trajectory <- data.frame(
      est        = c(age_mean, age_sd),
      age        = rep(age, 2),
      moment     = rep(c("mean", "sd"), each = nperson),
      modelcombi = modelcombi
    )
    return_info[["trajectory"]] <- trajectory
  }
  
  return(return_info)
  
}

fit_one_row <- function(row, raw_data, irt_model, age_variable, item_variables,
                        prior_knowledge  = NULL,parameter_fixed  = NULL, 
                        seed = NULL, iter_warmup = 500, iter_sampling = 500,
                        chains = 2, parallel_chains = NULL) {
  
  fit_one_model(
    model_specifications = list(
      irt_model = irt_model,
      age_model = "poly",
      poly_mean = row$mu,
      poly_sd   = row$sd, 
      prior_knowledge =  prior_knowledge, 
      parameter_fixed = parameter_fixed
    ),
    raw_data        = raw_data,
    age_variable    = age_variable,
    item_variables   = item_variables,
    traj_info = TRUE,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    chains          = chains,
    parallel_chains = parallel_chains,
    seed            = seed
  )
}

validate_model_specifications_list <- function(model_specifications_list){
  
  # must be a list
  if (!is.list(model_specifications_list)){
    stop("model_specifications_list must be a list of model_specifications entries (e.g. list(A = list(...), B = list(...))).")
  }
  
  # must contain at least two model entries
  if (length(model_specifications_list) < 2){
    stop("model_specifications_list must contain at least two model entries, but has ",
         length(model_specifications_list), ".")
  }
  
  # must be a list of lists (not a single flat model_specifications passed by mistake)
  is_sublist <- vapply(model_specifications_list, is.list, logical(1))
  if (!all(is_sublist)){
    bad <- names(model_specifications_list)[!is_sublist]
    if (length(bad) == 0) bad <- which(!is_sublist)
    stop("Each element of model_specifications_list must itself be a list. ",
         "Offending element(s): ", paste(bad, collapse = ", "))
  }
  
  # each sub-list must have non-NULL irt_model and age_model
  required <- c("irt_model", "age_model")
  required_poly <- c("poly_mean", "poly_sd")
  
  problems <- lapply(seq_along(model_specifications_list), function(i){
    mi <- model_specifications_list[[i]]
    nm <- names(model_specifications_list)[i]
    if (is.null(nm) || nm == "") nm <- paste0("[", i, "]")
    missing_fields <- required[vapply(required, function(f) is.null(mi[[f]]), logical(1))]
    
    if (identical(mi$age_model, "poly")){
      missing_fields <- c(
        missing_fields,
        required_poly[vapply(required_poly, function(f) is.null(mi[[f]]), logical(1))]
      )
    }
    
    if (length(missing_fields) > 0){
      paste0(nm, ": missing ", paste(missing_fields, collapse = ", "))
    } else {
      NULL
    }
  })

  
  problems <- unlist(problems)
  if (length(problems) > 0){
    stop("model_specifications_list has invalid entries:\n",
         paste(problems, collapse = "\n"),call. = FALSE)
  }
  
  invisible(TRUE)
}


check_dimension <- function(model_specifications){
  nmodels <- length(model_specifications)
  unidim <- vector(length = nmodels)
  for(n in seq_len(nmodels)){
    unidim[n] <- model_specifications[[n]]$irt_model %in% c("1PLnorm", "2PLnorm","3PLnorm", "GRMnorm")
  }
  if(any(unidim == TRUE) & any(unidim == FALSE)){
    stop("models in model_specifications do not imply the same dimensionality and can therefore not be compared.",call. = FALSE)
  }
    invisible(TRUE)
}

#' Compare age trajectory models for IRT norming models
#'
#' Fits a series of unidimensional IRT norming models that differ in how the
#' age trajectory of the latent trait is modeled, and compares them by their
#' expected log predictive density (ELPD), estimated with approximate
#' leave-one-out cross-validation (PSIS-LOO-CV).
#'
#' Age trajectories are modeled as polynomials of age, with separate degrees
#' for the mean (\code{mu}) and the standard deviation (\code{sd}) of the latent
#' trait. Candidate \code{(mu, sd)} combinations can be supplied as a grid
#' (\code{method = "grid"}) or found by forward selection
#' (\code{method = "forward"}). Optionally, a penalized-spline model is fitted
#' and compared against the polynomial models.
#'
#' @param method Character string, either \code{"grid"} (default) or
#'   \code{"forward"}.
#'   \itemize{
#'     \item \code{"grid"}: fits every \code{(mu, sd)} combination in
#'       \code{grid}.
#'     \item \code{"forward"}: greedy forward selection. See Details.
#'   }
#' @param grid A \code{data.frame} with columns \code{mu} and \code{sd}, giving
#'   the polynomial degrees to evaluate for the mean and the standard deviation
#'   of the age trajectory (e.g., \code{expand.grid(mu = 1:3, sd = 0:2)}).
#'   Required when \code{method = "grid"}; ignored when
#'   \code{method = "forward"}.
#'@param maxeval Integer. Upper limit on the number of polynomial models
#'   fitted during forward selection. Default is 10. Only used when
#'   \code{method = "forward"}. The initial model counts as one fit and each
#'   round adds two, and a new round is started only if both of its fits stay
#'   within the limit. At least \code{maxeval = 4} is therefore needed to
#'   evaluate any candidate beyond the starting model. The spline model (see
#'   \code{include_splines}) is not counted.
#' @param irt_model Character string naming the unidimensional IRT norming
#'   model to fit. One of:
#'   \describe{
#'     \item{\code{"1PLnorm"}}{One-parameter logistic (Rasch) model with
#'       age-dependent mean and variance. This is the default.}
#'     \item{\code{"2PLnorm"}}{Two-parameter logistic model with age-dependent
#'       mean and variance.}
#'     \item{\code{"3PLnorm"}}{Three-parameter logistic model with
#'       age-dependent mean and variance.}
#'     \item{\code{"GRMnorm"}}{Graded response model with age-dependent mean
#'       and variance.}
#'   }
#'   The hierarchical and testlet models (\code{"Testlet1PLnorm"},
#'   \code{"Testlet2PLnorm"}, \code{"HO1PLnorm-C"}, \code{"HO2PLnorm-C"},
#'   \code{"HO1PLnorm-P"}, \code{"HO2PLnorm-P"}, \code{"HO1PLnorm-F"},
#'   \code{"HO2PLnorm-F"}) are not supported by this function and cause an
#'   error. For these models, select the optimal polynomial degrees separately
#'   for each domain with a unidimensional model.
#' @param raw_data A \code{data.frame} or \code{matrix} with individual-level
#'   data, one row per individual. It must contain the age variable and the
#'   item response variables. Item responses may contain missing values
#'   (\code{NA}).
#' @param age_variable Column of \code{raw_data} containing age, given as a
#'   numeric index or a column name. Age is treated as continuous.
#' @param item_variables Columns of \code{raw_data} containing the item
#'   responses, given as a numeric or character vector of indices or names.
#' @param prior_knowledge Optional. A named list with prior means and prior
#'   standard deviations for each parameter group for which prior knowledge
#'   should be incorporated. For parameter \code{b} in the GRM, means and
#'   standard deviations must each be supplied as a list. Default is
#'   \code{NULL} (no prior knowledge). The same priors are used for all
#'   candidate models.
#' @param parameter_fixed Optional. A named list in which each name is a
#'   parameter group to be fixed and each element gives the value(s) to which
#'   it is fixed. Default is \code{NULL} (no fixed parameters). The same
#'   constraints are applied to all candidate models.
#' @param seed Optional integer random seed for reproducibility. If
#'   \code{NULL} (default), a seed is drawn at random. The same seed is used
#'   for every model fit.
#' @param iter_warmup Integer. Number of warm-up iterations per chain.
#'   Default is 500, which is usually sufficient for well-fitting data but may
#'   need to be increased for complex models.
#' @param iter_sampling Integer. Number of post-warm-up (sampling) iterations
#'   per chain. Default is 500; increase it for complex models or for more
#'   precise ELPD estimates.
#' @param chains Integer. Number of Markov chains. Default is 4.
#' @param parallel_chains Optional integer. Number of chains run in parallel.
#'   If \code{NULL} (default), as many chains as \code{chains} are run in
#'   parallel, limited by the number of available cores
#'   (see \code{\link[parallel]{detectCores}}).
#' @param include_splines Logical. If \code{TRUE} (default), a penalized-spline
#'   age model is fitted in addition to the polynomial models and included in
#'   the comparison.
#'
#' @details
#' \strong{Grid search.} Every row of \code{grid} is fitted once, and all
#' models are compared on equal footing.
#'
#' \strong{Forward selection.} The search starts at the simplest model
#' (\code{mu = 1}, \code{sd = 0}). In each round, the two neighboring models
#' \code{(mu + 1, sd)} and \code{(mu, sd + 1)} are fitted. If the better of the
#' two has a strictly higher ELPD than the current best model, it becomes the
#' new current best and the search continues. Otherwise the search stops, as
#' it also does when the \code{maxeval} limit is reached. Because the search is
#' greedy, it can end in a local optimum and may miss combinations that a grid
#' search would find.
#'
#' \strong{Model comparison.} All fitted models, including the spline model if
#' requested, are ranked by ELPD. The model with the highest ELPD receives rank
#' 1 and is flagged as \code{best} in the returned trajectory data. The ranking
#' is based on the point estimate of ELPD and does not account for uncertainty
#' in the differences (\code{ELPD_diff}) or for parsimony. Inspect
#' \code{ELPD_diff} before choosing the top-ranked model, and check
#' \code{paretok} for unreliable LOO estimates (commonly, values above 0.7).
#'
#' @return An object of class \code{"compare_info"}, a list with components:
#' \describe{
#'   \item{\code{best_model}}{A model comparison table with one row per fitted
#'     model, including \code{ELPD}, \code{ELPD_diff}, \code{modelcombi}
#'     (model label), the polynomial degrees \code{mu} and \code{sd}, and
#'     \code{rank}.}
#'   \item{\code{trajectory}}{A data frame with the estimated age trajectories
#'     (mean and SD of \code{theta} across the age range) for all fitted
#'     models, merged with each model's \code{rank}. The logical column
#'     \code{best} is \code{TRUE} for rows of the top-ranked model.}
#'   \item{\code{paretok}}{A matrix of Pareto-\eqn{\hat{k}} diagnostics from
#'     the LOO approximation, one column per fitted model.}
#'   \item{\code{best_configuration}}{A list describing the top-ranked model
#'     (\code{irt_model}, \code{age_model}, \code{poly_mean}, \code{poly_sd},
#'     \code{prior_knowledge}, \code{parameter_fixed}). \code{age_model} is
#'     \code{"poly"}, or \code{"splines"} if the spline model ranks first, in
#'     which case \code{poly_mean} and \code{poly_sd} are \code{NA}. It can be
#'     passed on to further model fitting.}
#' }
#'
#' @note
#' \code{grid} must be supplied when \code{method = "grid"}, otherwise the
#' function stops with an error. Fitting many Bayesian models can be
#' time-consuming; consider a small grid or \code{method = "forward"} for
#' initial exploration.
#'
#' @seealso \code{\link[loo]{loo}} for the underlying cross-validation method.
#'
#' @examples
#' \dontrun{
#' data(my_data)
#' # Grid search over a fixed set of polynomial degrees
#' my_grid <- expand.grid(mu = 1:3, sd = 0:2)
#' res_grid <- compare_age_models(
#'   method = "grid", grid = my_grid,
#'   raw_data = my_data, irt_model = "1PLnorm",
#'   age_variable = "age", item_variables = 1:25
#' )
#' res_grid$best_model
#'
#' # Forward selection, without the spline comparison
#' res_fwd <- compare_age_models(
#'   method = "forward",
#'   raw_data = my_data, irt_model = "1PLnorm",
#'   age_variable = "age", item_variables = 1:25,
#'   maxeval = 12, include_splines = FALSE
#' )
#' res_fwd$best_configuration
#' }
#'
#' @export
compare_age_models <- function(method = c("grid", "forward"),
                               grid = NULL,
                               maxeval = 10,
                               irt_model = c("1PLnorm","2PLnorm","3PLnorm","GRMnorm"),
                               raw_data,
                               age_variable,
                               item_variables,
                               prior_knowledge = NULL, 
                               parameter_fixed = NULL,
                               seed = NULL, 
                               iter_warmup = 500,
                               iter_sampling = 500,
                               chains = 4,
                               parallel_chains = NULL,
                               include_splines = TRUE){
  
  method <- match.arg(method)
  if (is.null(irt_model)) {
    stop("`irt_model` must be a single model name.")
  }
  check_rawinput(raw_data = raw_data, item_variables = item_variables, 
                 age_variable = age_variable)
  if(any(is.na(irt_model))){
    stop("The IRT model selected is NA, you have to select one of the 
         followingirt_model: \"1PLnorm\", \"2PLnorm\", \"3PLnorm\", \"GRMnorm\"")
  }
  if(any(irt_model == c("Testlet1PLnorm","Testlet2PLnorm",
                   "HO1PLnorm-C", "HO2PLnorm-C",
                   "HO1PLnorm-P", "HO2PLnorm-P",
                   "HO1PLnorm-F", "HO2PLnorm-F"))){
    stop("compare_age_models is for unidimensional IRT models. If you want to fit a 
         hierarchical IRT model, you should search for each domains optimal polynomial degree 
         seperatly.")
  }
  irt_model <- match.arg(irt_model)
  
  if (length(parallel_chains) == 0){
    parallel_chains <- min(chains, parallel::detectCores())
  }
  if (is.null(seed)){
    seed <- sample.int(.Machine$integer.max, 1)
  }
  if (method == "grid" && is.null(grid)){
    stop("`grid` must be supplied when method = 'grid'.")
  }
  
  # --- shared helper: fit a set of candidate (mu, sd) rows ---
  fit_grid <- function(grid){
    results <- lapply(seq_len(nrow(grid)), function(i){
      fit_one_row(
        row             = grid[i, ,drop = FALSE],
        raw_data        = raw_data,
        irt_model       = irt_model,
        age_variable    = age_variable,
        item_variables   = item_variables,
        prior_knowledge = prior_knowledge,
        parameter_fixed = parameter_fixed,
        seed            = seed,
        iter_warmup     = iter_warmup,
        iter_sampling   = iter_sampling,
        chains          = chains,
        parallel_chains = parallel_chains
      )
    })
    info <- list(
      loo_ELPD   = do.call(cbind, lapply(results, function(x) x$elpd_pointwise)),
      elpd_value = do.call(rbind, lapply(results, function(x) x$elpd_summary)),
      trajectory = do.call(rbind, lapply(results, function(x) x$trajectory)),
      paretok = do.call(cbind, lapply(results, function(x) x$paretok))
    )
    return(info)
  }
  
  # --- shared helper: finalize results (loo comparison, ranks, best flag) ---
  finalize <- function(loo_ELPD, elpd_value, trajectory, paretok, #best_configuration = NULL, 
                       irt_model,include_splines){
    if(include_splines){
      temp <- fit_one_model(
        model_specifications = list(
          irt_model = irt_model,
          age_model = "splines",
          prior_knowledge = prior_knowledge, 
          parameter_fixed = parameter_fixed
        ),
        raw_data        = raw_data,
        age_variable    = age_variable,
        item_variables   = item_variables,
        traj_info = TRUE,
        iter_warmup     = iter_warmup,
        iter_sampling   = iter_sampling,
        chains          = chains,
        parallel_chains = parallel_chains,
        seed            = seed
      )
        loo_ELPD <- cbind(loo_ELPD,temp$elpd_pointwise )
      elpd_value <- rbind(elpd_value,temp$elpd_summary )
      trajectory <- rbind(trajectory,temp$trajectory )
      paretok <- cbind(paretok, temp$paretok)
      
    }
    
    loo_info   <- list(elpd_value = elpd_value, loo_ELPD = loo_ELPD)
    best_model <- loo_comp(loo_info = loo_info)
    
    trajectory <- merge(
      trajectory,
      best_model[, c("modelcombi", "rank")],
      by = "modelcombi",
      all.x = TRUE
    )
    trajectory$best <- trajectory$rank == 1
    
    #if (is.null(best_configuration)){
      best_row <- best_model[1, ]
      
     # if(best_row$modelcombi == "splines"){
    #    best_agemod <- "splines"
    #  }
     # else{
    #    best_agemod <- "poly"
     # }
      best_configuration <- list(
        irt_model       = irt_model,
        age_model       = "poly",
        poly_mean       = best_row$mu,
        poly_sd         = best_row$sd,
        prior_knowledge = prior_knowledge,
        parameter_fixed = parameter_fixed
      )
      if(include_splines){
        best_configuration$age_model[is.na(best_configuration$poly_mean)] <- "splines"
      }
    #}
    
      structure(list(best_model = best_model, trajectory = trajectory,  paretok = paretok,
         best_configuration =  best_configuration),class = "compare_info")
  }
  
  # ============================================================
  # METHOD 1: evaluate a fixed grid of (mu, sd) combinations
  # ============================================================
  if (method == "grid"){
    out <- fit_grid(grid)
    #out$loo_ELPD
    return(finalize(out$loo_ELPD, out$elpd_value, out$trajectory, out$paretok,
                    irt_model = irt_model, include_splines = include_splines))
  }
  
  # ============================================================
  # METHOD 2: greedy forward selection
  # ============================================================
  # step 0: starting model (lowest polynomial order, mu = 1, sd = 0)
  grid0 <- data.frame(mu = 1, sd = 0)
  out0  <- fit_grid(grid0)
  
  loo_ELPD   <- out0$loo_ELPD
  elpd_value <- out0$elpd_value
  trajectory <- out0$trajectory
  paretok <- out0$paretok
  
  best_row  <- elpd_value[nrow(elpd_value), ]
  best_elpd <- best_row$ELPD
  best_mu   <- best_row$poly_mean
  best_sd   <- best_row$poly_sd
  
  n_eval   <- 1
  improved <- TRUE
  
  while (improved && n_eval + 2 < maxeval){
    
    candidates <- data.frame(mu = c(best_mu + 1, best_mu),
                             sd = c(best_sd,     best_sd + 1))
    
    out <- fit_grid(candidates)
    
    loo_ELPD   <- cbind(loo_ELPD, out$loo_ELPD)
    elpd_value <- rbind(elpd_value, out$elpd_value)
    trajectory <- rbind(trajectory, out$trajectory)
    paretok    <- cbind(paretok, out$paretok )
    n_eval     <- n_eval + nrow(candidates)
    
    new_rows <- utils::tail(elpd_value, nrow(candidates))
    best_new <- new_rows[which.max(new_rows$ELPD), ]
    
    if (best_new$ELPD > best_elpd){
      best_elpd <- best_new$ELPD
      best_mu   <- best_new$poly_mean
      best_sd   <- best_new$poly_sd
      improved  <- TRUE
    } else {
      improved <- FALSE
    }
  }
  
  #best_configuration <- list(irt_model       = irt_model,
  #                           age_model       = "poly",
  #                           poly_mean       = best_mu ,
  #                           poly_sd         = best_sd,
  #                           prior_knowledge = prior_knowledge,
  #                           parameter_fixed = parameter_fixed)
  
  finalize(loo_ELPD, elpd_value, trajectory,paretok, #best_configuration,
           include_splines = include_splines, irt_model = irt_model)
}





#' Visualize results of an age model comparison
#'
#' Creates a diagnostic plot from the output of
#' \code{\link{compare_age_models}}. The results can be inspected either from
#' the perspective of model fit (ELPD differences relative to the best model,
#' across polynomial degree combinations) or from the perspective of the
#' estimated age trajectories themselves.
#'
#' @param compare_info An object of class \code{"compare_info"}, as returned by
#'   \code{\link{compare_age_models}}. Any other input causes an error.
#' @param perspective Character string, either \code{"ELPD"} (default) or
#'   \code{"trajectory"}.
#'   \itemize{
#'     \item \code{"ELPD"}: plots the ratio of the ELPD difference to its
#'       standard error for each fitted model, by polynomial degree of the
#'       mean. See Details.
#'     \item \code{"trajectory"}: plots the estimated age trajectories of all
#'       fitted models, with the top-ranked model(s) highlighted. See Details.
#'   }
#' @param highlight_ranks Numeric value of length 1. Models with a rank less
#'   than or equal to this value are highlighted: labeled with their rank
#'   under \code{perspective = "ELPD"}, and drawn as thick colored lines
#'   (colored by rank) under \code{perspective = "trajectory"}. Default is
#'   1, which highlights only the best model. If the value is not a single
#'   number, a warning is issued and the default is used.
#'
#' @details
#' \strong{perspective = "ELPD".} Each point is one fitted model. The x-axis
#' is the polynomial degree for the mean (\code{mu}), the point color encodes
#' the polynomial degree for the standard deviation (\code{sd}), and the
#' y-axis shows the ELPD difference from the best model divided by its
#' standard error (\code{ELPD_diff / SE}).
#' Values close to zero indicate models that are hard to distinguish from the
#' best model, which is plotted at 0; strongly negative values indicate
#' clearly worse fit. Models within the top
#' \code{highlight_ranks} ranks are labeled with their rank.
#' If the spline model was included in
#' \code{\link{compare_age_models}}, it has no polynomial degrees and is
#' plotted at the x position and in the color category \code{"splines"}.
#'
#' \strong{perspective = "trajectory".} Each line is the estimated age
#' trajectory of one fitted model, shown separately for the mean and the
#' standard deviation of the latent trait. Models within the top
#' \code{highlight_ranks} ranks are highlighted in color, while all other
#' models are drawn faintly in the background. This shows whether top-ranked
#' models agree on the shape of the trajectory or diverge.
#' 
#' @return A \code{ggplot} object, which can be printed directly or further
#'   customized with additional \pkg{ggplot2} layers.
#'
#' @examples
#' \dontrun{
#' data(my_data)
#' res <- compare_age_models(
#'   method = "grid", grid = expand.grid(mu = 1:3, sd = 0:2),
#'   raw_data = my_data, irt_model = "1PLnorm",
#'   age_variable = "age", item_variables = 1:25
#' )
#'
#' # ELPD comparison across (mu, sd) combinations
#' visual_age_models(res, perspective = "ELPD")
#'
#' # Estimated age trajectories, highlighting the three best models
#' visual_age_models(res, perspective = "trajectory", highlight_ranks = 3)
#' }
#'
#' @seealso \code{\link{compare_age_models}}
#'
#' @export
visual_age_models <- function(compare_info, perspective = c("ELPD", "trajectory"), 
                                        highlight_ranks = 1){
  perspective <-  match.arg(perspective)
  if(!is.numeric(highlight_ranks) | length(highlight_ranks) != 1){
    warning("highlight_ranks has to be a numeric value of length 1. The default 
            of 1 was used instead.")
    highlight_ranks <- 1
  }
  if (!inherits(compare_info, "compare_info")) {
    stop("`compare_info` must be created by `compare_age_models()`.")
  }
  if(perspective == "ELPD"){
    best_model <- compare_info$best_model
    best_model$ranked <- as.factor(best_model$rank)
    if(any(is.na(best_model$poly_mean))){
      best_model$poly_mean[is.na(best_model$poly_mean)] <- "splines"
      best_model$poly_sd[is.na(best_model$poly_sd)] <- "splines"
    }
    best_model$poly_sd <- as.factor(best_model$poly_sd)
    best_model$SE[best_model$SE == 0] <- 1
    best_model$ratio <- best_model$ELPD_diff/best_model$SE
    p <- ggplot2::ggplot(best_model) +  ggplot2::theme_bw() +
      ggplot2::geom_point( ggplot2::aes(.data$poly_mean, .data$ratio, color = .data$poly_sd,size =1.2  )) +
      ggplot2::geom_text(
        data = subset(best_model,rank <highlight_ranks + 1),   # only these points get labels
        ggplot2::aes(.data$poly_mean, .data$ratio, label = .data$rank),              # the number to display
        hjust = 2                    # nudge label above the point
      )+
      ggplot2::xlab("polynomial degree for mean") +
      ggplot2::theme(text =  ggplot2::element_text(size=13,family = "serif"), 
                     axis.line =  ggplot2::element_line(color='black'),
                     plot.background =  ggplot2::element_blank(),
                     panel.grid.minor =  ggplot2::element_blank(),
                     panel.grid.major =  ggplot2::element_blank(), 
                     legend.position = "bottom") + 
      ggplot2::ylab("ELPD difference/SE") +  ggplot2::guides(size = "none") +
      ggplot2::labs(color ="degree for sd" )
    
  }
  if(perspective == "trajectory"){
    trajectory <- compare_info$trajectory
    trajectory$ranked <- as.factor(trajectory$rank)
    p <-  ggplot2::ggplot(trajectory) +  ggplot2::theme_bw()+
      ggplot2::geom_line( ggplot2::aes(.data$age, .data$est, group = .data$modelcombi), alpha = 0.2) +
      ggplot2::geom_line( ggplot2::aes(.data$age, .data$est, color = .data$ranked), 
                          data = trajectory[trajectory$rank  <= highlight_ranks,],linewidth = 1.2) +
      ggplot2::facet_wrap(~.data$moment, scale = "free") +  ggplot2::ylab("value") +
      ggplot2::theme(text =  ggplot2::element_text(size=13,family = "serif"), 
                     axis.line =  ggplot2::element_line(color='black'),
                     plot.background =  ggplot2::element_blank(),
                     panel.grid.minor =  ggplot2::element_blank(),
                     panel.grid.major =  ggplot2::element_blank()) 
  }
  return(p)
}


#' Compare IRT norming models via leave-one-out cross-validation
#'
#' Fits a set of user-specified IRT norming models, which may differ in the
#' IRT model, the age trajectory specification, or both, and compares them by
#' their expected log predictive density (ELPD), estimated with approximate
#' leave-one-out cross-validation (LOO-CV). Unlike
#' \code{\link{compare_age_models}}, which searches over polynomial degrees
#' for a single IRT model, this function compares arbitrary, independently
#' specified models directly.
#'
#' @param model_specifications A list of model specifications. Each element is
#'   itself a list describing one model, with fields such as \code{irt_model},
#'   \code{age_model}, \code{poly_mean}, \code{poly_sd},
#'   \code{prior_knowledge}, and \code{parameter_fixed}. See
#'   \code{\link{estimate_norming_model}} for all accepted fields. All models
#'   must have the same dimensionality (all unidimensional or all
#'   multidimensional), because ELPD values are only comparable across models
#'   for the same data and outcome structure.
#' @param raw_data A \code{data.frame} with individual-level data, one row per
#'   individual. It must contain the age variable and the item response
#'   variables. Item responses may contain missing values (\code{NA}).
#' @param age_variable Character string giving the name of the column in
#'   \code{raw_data} that contains age. Age is treated as continuous.
#' @param item_variables Character vector giving the names of the columns in
#'   \code{raw_data} that contain the item responses.
#' @param seed Optional integer random seed for reproducibility. If
#'   \code{NULL} (default), a seed is drawn at random. The same seed is used
#'   for every model.
#' @param iter_warmup Integer. Number of warm-up iterations per chain. Default
#'   is 500, which is usually sufficient for well-fitting data but may need to
#'   be increased for complex models.
#' @param iter_sampling Integer. Number of post-warm-up (sampling) iterations
#'   per chain. Default is 500; increase it for complex models or for more
#'   precise ELPD estimates.
#' @param chains Integer. Number of Markov chains per model. Default is 2.
#' @param parallel_chains Optional integer. Number of chains run in parallel.
#'   If \code{NULL} (default), as many chains as \code{chains} are run in
#'   parallel, limited by the number of available cores
#'   (see \code{\link[parallel]{detectCores}}).
#'
#' @details
#' All models are fitted to the same data with the same \code{seed} and
#' sampler settings; only the model specification differs. Models are
#' labeled by their position in \code{model_specifications} and their IRT
#' model (e.g., \code{"1_1PLnorm"}), so list names are not used.
#'
#' If a model fails to fit, a warning is issued and the model is left out of
#' the comparison. The function stops with an error if fewer than two models
#' are fitted successfully.
#'
#' Models are ranked by ELPD, with rank 1 given to the model with the highest
#' ELPD. The ranking is based on the point estimate and does not account for
#' uncertainty in the differences, so inspect \code{ELPD_diff} and its
#' standard error before choosing the top-ranked model. Also check the Pareto
#' \eqn{\hat{k}} diagnostics (values above 0.7 indicate unreliable LOO
#' estimates).
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{best_model}}{A data frame with one row per successfully
#'     fitted model, containing the ELPD (\code{ELPD}), the difference to the
#'     best model (\code{ELPD_diff}), and the \code{rank}, along with columns
#'     identifying the model (e.g., \code{modelcombi}).}
#'   \item{\code{diagnostics}}{A data frame with the model label
#'     (\code{model}), \code{irt_model}, \code{age_model}, polynomial degrees
#'     (\code{poly_mean}, \code{poly_sd}; \code{NA} for non-polynomial age
#'     models), and the sampling diagnostics \code{pareto_k}, \code{neff},
#'     and \code{reff}. There is one row per observation and model, because
#'     \code{pareto_k} is computed for each observation.}
#' }
#'
#' @examples
#' \dontrun{
#' data(my_data)
#' model_specifications <- list(
#'   list(irt_model = "1PLnorm", age_model = "poly",
#'        poly_mean = 2, poly_sd = 1),
#'   list(irt_model = "2PLnorm", age_model = "poly",
#'        poly_mean = 2, poly_sd = 1)
#' )
#'
#' result <- compare_IRTnorm_models(
#'   model_specifications = model_specifications,
#'   raw_data = my_data,
#'   age_variable = "age",
#'   item_variables = paste0("item", 1:25)
#' )
#' result$best_model
#' }
#'
#' @seealso \code{\link{compare_age_models}}
#'
#' @export
compare_IRTnorm_models <- function(model_specifications, raw_data, age_variable,
                                   item_variables, seed = NULL,
                                   iter_warmup = 500, iter_sampling = 500,
                                   chains = 2, parallel_chains = NULL) {
  
  validate_model_specifications_list(model_specifications)
  check_dimension(model_specifications)
  check_rawinput(raw_data = raw_data, item_variables = item_variables, 
                 age_variable = age_variable)
  if (is.null(parallel_chains)) {
    parallel_chains <- min(chains, parallel::detectCores())
  }
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1)
  }
  
  # Unique, informative model labels
  model_ids <- vapply(seq_along(model_specifications), function(n) {
    paste0(n, "_", model_specifications[[n]]$irt_model)
  }, character(1))
  
  # Fit all models; a failed fit returns NULL instead of aborting everything
  fits <- lapply(seq_along(model_specifications), function(n) {
    tryCatch(
      fit_one_model(
        model_specifications = model_specifications[[n]],
        raw_data = raw_data, age_variable = age_variable,
        item_variables = item_variables, seed = seed,
        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
        chains = chains, parallel_chains = parallel_chains
      ),
      error = function(e) {
        warning("Model ", model_ids[n], " failed: ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )
  })
  
  # Keep only successful fits, with specs and labels aligned
  ok       <- !vapply(fits, is.null, logical(1))
  if (sum(ok) < 2) stop("Fewer than two models fitted successfully; cannot compare.")
  fits     <- fits[ok]
  specs_ok <- model_specifications[ok]
  ids_ok   <- model_ids[ok]
  names(fits) <- ids_ok
  
  # Pointwise ELPD matrix (one named column per model)
  loo_ELPD <- do.call(cbind, lapply(fits, `[[`, "elpd_pointwise"))
  colnames(loo_ELPD) <- ids_ok
  
  # ELPD summaries (one named row per model)
  elpd_value <- do.call(rbind, lapply(fits, `[[`, "elpd_summary"))
  rownames(elpd_value) <- ids_ok
  
  # Diagnostics, including age-model information
  diagn <- do.call(rbind, lapply(seq_along(fits), function(n) {
    spec <- specs_ok[[n]]
    
    if (identical(spec$age_model, "poly")) {
      poly_mean_temp <- spec$poly_mean
      poly_sd_temp   <- spec$poly_sd
    } else {
      poly_mean_temp <- NA
      poly_sd_temp   <- NA
    }
    
    data.frame(model     = ids_ok[n],
               irt_model = spec$irt_model,
               age_model = spec$age_model,
               poly_mean = poly_mean_temp,
               poly_sd   = poly_sd_temp,
               pareto_k  = fits[[n]]$paretok,
               neff      = fits[[n]]$neff,
               reff      = fits[[n]]$reff)
  }))
  
  best_model <- loo_comp(loo_info = list(elpd_value = elpd_value,
                                         loo_ELPD = loo_ELPD))
  
  list(best_model = best_model, diagnostics = diagn)
}
