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
  
  d <- loo_info[[2]][,reference] - loo_info[[2]][,best_list]
  loo_order$SE <- sqrt(nrow(d) * apply(d, 2, stats::var))
  loo_order$rank <- seq_len(nrow(loo_order))
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
  
  problems <- lapply(seq_along(model_specifications_list), function(i){
    mi <- model_specifications_list[[i]]
    nm <- names(model_specifications_list)[i]
    if (is.null(nm) || nm == "") nm <- paste0("[", i, "]")
    
    missing_fields <- required[vapply(required, function(f) is.null(mi[[f]]), logical(1))]
    if (length(missing_fields) > 0){
      paste0(nm, ": missing ", paste(missing_fields, collapse = ", "))
    } else {
      NULL
    }
  })
  
  problems <- unlist(problems)
  if (length(problems) > 0){
    stop("model_specifications_list has invalid entries:\n",
         paste(problems, collapse = "\n"))
  }
  
  invisible(TRUE)
}


check_dimension <- function(model_specifications){
  nmodels <- length(model_specifications)
  unidim <- vector(length = nmodels)
  for(n in seq_len(nmodels)){
    unidim <- model_specifications[[n]]$irt_model %in% c("1PLnorm", "2PLnorm","3PLnorm", "GRMnorm")
  }
  if(any(unidim == TRUE) & any(unidim == FALSE)){
    stop("models in model_specifications do not imply the same dimensionallity and can therefore not be compared.")
  }
    invisible(TRUE)
}

#' Compare IRTnorm models modeled with different polynomial degrees 
#' via Grid Search or via Forward Selection
#'
#' Fits a series of IRTnorm with polynomial age trajectories (varying the
#' polynomial degree for the mean, \code{mu}, and standard deviation, \code{sd},
#' of the trajectory) and compares them using leave-one-out cross-validation
#' (LOO-ELPD). Models can be evaluated either over a user-supplied grid of
#' \code{(mu, sd)} combinations, or via a forward-selection search that
#' starts from the simplest model and incrementally increases polynomial
#' complexity as long as model fit improves.
#'
#' @param method Character string, either \code{"grid"} or \code{"forward"}.
#'   \code{"grid"} evaluates every \code{(mu, sd)} combination supplied in
#'   \code{grid}. \code{"forward"} performs greedy forward selection, starting
#'   from \code{mu = 1, sd = 0} and at each step trying to increase either
#'   \code{mu} or \code{sd} by one, keeping whichever candidate improves the
#'   ELPD, and stopping when neither candidate improves on the current best
#'   model or \code{maxeval} model fits have been reached.
#' @param grid A data frame with columns \code{mu} and \code{sd} giving the
#'   polynomial degree combinations to evaluate. Required (and only used) when
#'   \code{method = "grid"}; ignored otherwise.
#'@param raw_data A \code{data.frame} or \code{matrix} containing individual-level
#'   data. Each row represents one individual. The object must contain at least
#'   an age variable and item response variables. Item responses may
#'   contain missing values (\code{NA}).
#'
#'@param age_variable Either a \code{numeric} index or a \code{character} string specifying
#'   the column in \code{raw_data} that contains the individuals' age information.
#'   Age is assumed to be continuous.
#'
#'@param item_variables A \code{numeric} or \code{character} vector specifying the columns in
#'   \code{raw_data} that contain the item response variables.
#'
#'@param irt_model A \code{character} scalar indicating for which IRT norming model the
#'   data should be prepared. The currently supported options are:
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
#' @param prior_knowledge Only required if prior knowledge should be incorporated
#' in the model. A list of prior means and prior standard deviations has to be indicated,
#' for each parameter group for which prior knowledge is to be incorporated. For b in the GRM model
#' this has to be a list per means and standard deviations.
#' Below in the example an example list is used.
#' @param parameter_fixed Only required if model parameters should be fixed.
#' A list which indicates the parameter groups that should be fixed, with the list names
#' and the list elements indicating the values to which the parameters should be fixed.
#' @param maxeval Integer. Maximum number of individual model fits to perform
#'   during forward selection before stopping, regardless of whether
#'   improvement has plateaued. Only used when \code{method = "forward"}.
#'   Defaults to 10.

#' @param seed An optional \code{integer} specifying the random seed for
#'   reproducibility. If \code{NULL} (default), a random seed is generated
#'   internally.
#'
#' @param iter_warmup An \code{integer} specifying the number of warm-up
#'   iterations per chain used by the MCMC sampler.  The default value of 500 is 
#'   sufficient for well-fitting data, but it may need to be increased for 
#'   complex models or to improve precision.  
#'
#' @param iter_sampling An \code{integer} specifying the number of post
#'   warm-up (sampling) iterations per chain.  The default value of 500 is 
#'   sufficient for well-fitting data, but it may need to be increased for
#'    complex models or to improve precision.  
#'    
#' @param chains An \code{integer} specifying the number of Markov chains to
#'   run in parallel. The default is 4 chains.
#'
#' @param parallel_chains An optional \code{integer} specifying the number of
#'   CPU cores used for parallel sampling. If \code{NULL} (default), all
#'   available cores detected by \code{\link[parallel]{detectCores}} are used.
#'
#'@param include_splines Logical, if penalized splines should be compared to the 
#'polynomial selection or not. 
#' @details
#' Returns LOO-ELPD values (\code{loo_extr}) along with the
#' estimated age trajectory (mean and SD of \code{theta} across the age
#' range).
#'
#' Under \code{method = "grid"}, every row of \code{grid} is fit once and all
#' models are compared on equal footing.
#'
#' Under \code{method = "forward"}, the search begins at the simplest model
#' (\code{mu = 1, sd = 0}) and proceeds greedily: at each iteration, the two
#' neighboring candidates (\code{mu + 1, sd} and \code{mu, sd + 1}) are fit
#' and compared against the current best ELPD. The search advances to
#' whichever candidate scores higher, or stops if neither candidate improves
#' on the current best model. The search also stops once \code{maxeval} total
#' model fits have been performed, even if improvement has not yet
#' plateaued. Note that \code{maxeval} counts individual model fits (the
#' initial model plus two per round), not the number of forward-selection
#' rounds.
#'
#' After fitting, all candidate models are compared based on their ELPD 
#' differences and ranks. The model with rank 1 is flagged as
#' \code{best} in the returned \code{trajectory} data frame.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{best_model}}{A model comparison
#'     table including \code{ELPD}, \code{ELPD_diff}, \code{modelcombi}, and
#'     \code{rank} for every fitted model.}
#'   \item{\code{trajectory}}{A data frame combining the estimated age
#'     trajectories (mean and SD) for all fitted models, merged with model
#'     rank, and with a logical \code{best} column flagging rows belonging to
#'     the top-ranked (rank 1) model.}
#'   \item{\code{best_configuration}}{A list describing the winning model
#'     configuration (\code{irt_model}, \code{age_model}, \code{poly_mean},
#'     \code{poly_sd}, \code{prior_knowledge}, \code{parameter_fixed}).
#'     \strong{Only populated when \code{method = "forward"}}; under
#'     \code{method = "grid"} this element is \code{NULL}, since no single
#'     configuration is privileged as a search outcome.}
#' }
#'
#' @note
#' When \code{method = "grid"}, \code{grid} must be supplied or the function
#' stops with an error. When \code{method = "forward"}, \code{grid} is
#' ignored if supplied.
#'
#' @examples
#' \dontrun{
#' # Grid search over a fixed set of polynomial degrees
#' my_grid <- expand.grid(mu = 1:3, sd = 0:2)
#' res_grid <- compare_age_models(
#'   method = "grid", grid = my_grid,
#'   raw_data = my_data, irt_model = "1PLnorm",
#'   age_variable = "age", item_variables = 1:20
#' )
#'
#' # forward selection
#' res_fwd <- compare_age_models(
#'   method = "forward",
#'   raw_data = my_data, irt_model = "1PLnorm",
#'   age_variable = "age", item_variables = 1:20,
#'   maxeval = 12
#' )
#' res_fwd$best_configuration
#' }
#'
#' @export
compare_age_models <- function(method = c("grid", "forward"),
                                       grid = NULL,
                                       raw_data, irt_model,
                                       age_variable, item_variables,
                                       maxeval = 10,
                                       prior_knowledge = NULL, parameter_fixed = NULL,
                                       seed = NULL, iter_warmup = 500,
                                       iter_sampling = 500,
                                       chains = 2, parallel_chains = NULL,
                                       include_splines = TRUE){
  
  method <- match.arg(method)
  
  if(any(irt_model == c("Testlet1PLnorm","Testlet2PLnorm",
                   "HO1PLnorm-C", "HO2PLnorm-C",
                   "HO1PLnorm-P", "HO2PLnorm-P",
                   "HO1PLnorm-F", "HO2PLnorm-F"))){
    stop("compare_age_models is for unidimensional IRT models. If you want to fit a 
         hierarchical IRT model, you should search for each domains optimal polynomial degree 
         seperatly.")
  }
  if (length(parallel_chains) == 0){
    parallel_chains <- parallel::detectCores()
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
  finalize <- function(loo_ELPD, elpd_value, trajectory, paretok, best_configuration = NULL, irt_model,include_splines){
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
    
    if (is.null(best_configuration)){
      best_row <- best_model[1, ]
      
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
    }
    
    
    list(best_model = best_model, trajectory = trajectory,  paretok = paretok,
         best_configuration =  best_configuration)
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
  
  while (improved && n_eval < maxeval){
    
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
  
  best_configuration <- list(irt_model       = irt_model,
                             age_model       = "poly",
                             poly_mean       = best_mu ,
                             poly_sd         = best_sd,
                             prior_knowledge = prior_knowledge,
                             parameter_fixed = parameter_fixed)
  
  finalize(loo_ELPD, elpd_value, trajectory,paretok, best_configuration,
           include_splines = include_splines, irt_model = irt_model)
}





#' Visualize Polynomial Model Selection Results
#'
#' Creates a diagnostic plot from the output of
#' \code{\link{compare_age_models}}, allowing inspection of model
#' comparison results either from the perspective of model fit (ELPD
#' differences across polynomial degree combinations) or from the
#' perspective of the estimated age trajectories themselves.
#'
#' @param compare_info A list as returned by
#'   \code{\link{compare_age_models}}
#' @param perspective Character string, either \code{"ELPD"} or
#'   \code{"trajectory"}. \code{"ELPD"} plots ELPD differences for each
#'   \code{(mu, sd)} combination, faceted by \code{sd}, with the
#'   top-ranked model(s) highlighted. \code{"trajectory"} plots the
#'   estimated age trajectories for every fitted model, faceted by moment
#'   (\code{mean}/\code{sd}), with the top-ranked model(s) highlighted.
#'   Defaults to \code{"ELPD"}. Partial matching is supported via
#'   \code{\link[base]{match.arg}}.
#'@param  highlight_ranks Numeric value indicating the top highlight_ranks to be 
#' highlighted in color. 
#' @details
#' Under \code{perspective = "ELPD"}, each point represents one fitted
#' model's ELPD difference from the best-performing model, plotted against
#' \code{mu} and faceted by \code{sd}. Models within the top
#' \code{highlight_ranks} ranks are colored by rank; all other models are
#' shown as faint unfilled points for context.
#'
#' Under \code{perspective = "trajectory"}, each line represents one fitted
#' model's estimated age trajectory (separately for the mean and SD moments,
#' via \code{facet_wrap(~moment)}). All trajectories are drawn faintly in
#' the background, with trajectories from models within the top
#' \code{highlight_ranks} ranks overlaid in color and at greater line
#' thickness, making it easy to see whether top-ranked models agree on the
#' shape of the trajectory or diverge.
#'
#' @return A \code{ggplot} object, which can be further customized or
#'   printed directly.
#'
#'
#' @examples
#' \dontrun{
#' res <- compare_age_models(
#'   method = "grid", grid = expand.grid(mu = 1:3, sd = 0:2),
#'   raw_data = my_data, irt_model = "1PLnorm",
#'   age_variable = "age", item_variables = c("sex")
#' )
#'
#' # ELPD comparison across (mu, sd) combinations
#' visual_polynomial_selection(res, perspective = "ELPD")
#'
#' # Estimated age trajectories, with best model(s) highlighted
#' visual_polynomial_selection(res, perspective = "trajectory")
#' }
#'
#' @seealso \code{\link{compare_age_models}}
#'
#' @export
visual_polynomial_selection <- function(compare_info, perspective = c("ELPD", "trajectory"), 
                                        highlight_ranks = 1){
  perspective <-  match.arg(perspective)
  if(perspective == "ELPD"){
    best_model <- compare_info$best_model
    best_model$ranked <- as.factor(best_model$rank)
    if(any(is.na(best_model$poly_mean))){
      best_model$poly_mean[is.na(best_model$poly_mean)] <- "splines"
      best_model$poly_sd[is.na(best_model$poly_sd)] <- "splines"
    }
    best_model$poly_sd <- as.factor(best_model$poly_sd)
    best_model$SE[best_model$SE == 0] <- 1
    p <- ggplot2::ggplot(best_model) +  ggplot2::theme_bw() +
      ggplot2::geom_point( ggplot2::aes(poly_mean, ELPD_diff/SE, color = poly_sd,size =1.2  )) +
      ggplot2::geom_text(
        data = subset(best_model,rank <highlight_ranks + 1),   # only these points get labels
        ggplot2::aes(poly_mean, ELPD_diff/SE, label = rank),              # the number to display
        hjust = 2                    # nudge label above the point
      )+
      #ggplot2:: geom_point( ggplot2::aes(poly_mean, ELPD_diff/SE, shape = ranked, color = poly_sd, size = 1.2), 
      #                     data = best_model[best_model$rank <= highlight_ranks,]) +
      ggplot2::xlab("polynomial degree for mean") +
      ggplot2::theme(text =  ggplot2::element_text(size=13,family = "serif"), 
                     axis.line =  ggplot2::element_line(color='black'),
                     plot.background =  ggplot2::element_blank(),
                     panel.grid.minor =  ggplot2::element_blank(),
                     panel.grid.major =  ggplot2::element_blank(), 
                     legend.position = "bottom") + 
      ggplot2::ylab("ELPD difference/SE") +  ggplot2::guides(size = "none") +
      ggplot2::labs(color ="degree for sd", shape = "rank"
      )
    
  }
  if(perspective == "trajectory"){
    trajectory <- compare_info$trajectory
    trajectory$ranked <- as.factor(trajectory$rank)
    p <-  ggplot2::ggplot(trajectory) +  ggplot2::theme_bw()+
      ggplot2::geom_line( ggplot2::aes(age, est, group = modelcombi), alpha = 0.2) +
      ggplot2::geom_line( ggplot2::aes(age, est, color = ranked), 
                          data = trajectory[trajectory$rank  <= highlight_ranks,],linewidth = 1.2) +
      ggplot2::facet_wrap(~moment, scale = "free") +  ggplot2::ylab("value") +
      ggplot2::theme(text =  ggplot2::element_text(size=13,family = "serif"), 
                     axis.line =  ggplot2::element_line(color='black'),
                     plot.background =  ggplot2::element_blank(),
                     panel.grid.minor =  ggplot2::element_blank(),
                     panel.grid.major =  ggplot2::element_blank()) 
  }
  return(p)
}


#' Compare Multiple IRTnorm models via Leave-One-Out Cross-Validation
#'
#' Fits a set of user-specified IRTnorm models (e.g. differing in IRT model
#' type, or age trajectory specification) and
#' compares them using leave-one-out cross-validation (LOO-ELPD). Unlike
#' \code{\link{compare_age_models}}, which searches over polynomial
#' degree combinations for a single model specification,
#' \code{compare_IRTnorm_models} allows arbitrary, independently specified
#' model configurations to be compared directly against one another.
#'
#' @param model_specifications A list of model specifications, where each element is
#'   itself a list describing one model's properties (e.g. \code{irt_model},
#'   \code{age_model}, \code{poly_mean}, \code{poly_sd}, and any other fields
#'   accepted by \code{\link{estimate_norming_model}}).
#'   \strong{All models in \code{model_specifications} must be of the same
#'   dimensionality} — i.e. all unidimensional or all multidimensional.
#'   Comparing models of different dimensionality is not supported, since
#'   LOO-ELPD comparisons assume a common outcome/likelihood structure
#'   across models.
#'   
#'@param raw_data A \code{data.frame} or \code{matrix} containing individual-level
#'   data. Each row represents one individual. The object must contain at least
#'   an age variable and item response variables. Item responses may
#'   contain missing values (\code{NA}).
#'
#'@param age_variable Either a \code{numeric} index or a \code{character} string specifying
#'   the column in \code{raw_data} that contains the individuals' age information.
#'   Age is assumed to be continuous.
#'
#'@param item_variables A \code{numeric} or \code{character} vector specifying the columns in
#'   \code{raw_data} that contain the item response variables.
#'
#' @param seed An optional \code{integer} specifying the random seed for
#'   reproducibility. If \code{NULL} (default), a random seed is generated
#'   internally.
#'
#' @param iter_warmup An \code{integer} specifying the number of warm-up
#'   iterations per chain used by the MCMC sampler.  The default value of 500 is 
#'   sufficient for well-fitting data, but it may need to be increased for 
#'   complex models or to improve precision.  
#'
#' @param iter_sampling An \code{integer} specifying the number of post
#'   warm-up (sampling) iterations per chain.  The default value of 500 is 
#'   sufficient for well-fitting data, but it may need to be increased for
#'    complex models or to improve precision.  
#'
#' @param chains An \code{integer} specifying the number of Markov chains to
#'   run in parallel. The default is 4 chains.
#'
#' @param parallel_chains An optional \code{integer} specifying the number of
#'   CPU cores used for parallel sampling. If \code{NULL} (default), all
#'   available cores detected by \code{\link[parallel]{detectCores}} are used.
#'
#' @details
#' For every model the data is prepared, the model is fitted, 
#' and LOO-ELPD values are extracted and reported. All fitted models share
#' the same \code{raw_data}, \code{age_variable}, \code{item_variables}, and
#' MCMC settings (\code{seed}, \code{iter_warmup}, \code{iter_sampling},
#' \code{chains}, \code{parallel_chains}); only the model specification
#' itself varies across entries in \code{model_specifications}.
#' For direct model comparison the ELPD differences and the resulting ranks
#' across the full set of compared models are reported.
#'
#' @return A data frame giving the LOO-ELPD,
#'   \code{ELPD_diff}, and \code{rank} for each fitted model, along with
#'   model-identifying columns (e.g. \code{modelcombi}, \code{irt_model},
#'   \code{age_model}) carried over from each model's \code{elpd_summary}.
#'
#' @examples
#' \dontrun{
#' model_specifications <- list(
#'   A = list(irt_model = "1PLnorm", age_model = "poly",
#'            poly_mean = 2, poly_sd = 1),
#'   B = list(irt_model = "1PLnorm", age_model = "splines"),
#'   C = list(irt_model = "2PLnorm", age_model = "poly",
#'            poly_mean = 2, poly_sd = 1)
#' )
#'
#' result <- compare_IRTnorm_models(
#'   model_specifications   = model_specifications,
#'   raw_data     = my_data,
#'   age_variable = "age",
#'   item_variables = 1:50
#' )
#' }
#'
#' @seealso \code{\link{compare_age_models}}
#'
#' @export
compare_IRTnorm_models <- function(model_specifications, raw_data, age_variable, item_variables, 
                                   seed = NULL, iter_warmup = 500,
                                   iter_sampling = 500,
                                   chains = 2, parallel_chains = NULL){

  validate_model_specifications_list(model_specifications)
  check_dimension(model_specifications)
  
  if (length(parallel_chains) == 0){
    parallel_chains <- parallel::detectCores()
  }
  if (is.null(seed)){
    seed <- sample.int(.Machine$integer.max, 1)
  }
  
  loo_ELPD <- NULL
  elpd_value <- NULL
  diagn <- NULL
  for(n in seq_len(length(model_specifications))){
    info_model <- fit_one_model(model_specifications = model_specifications[[n]], 
                                raw_data = raw_data, age_variable =age_variable, item_variables =item_variables,
                                seed = seed, iter_warmup = iter_warmup, 
                                iter_sampling = iter_sampling, 
                                chains = chains, parallel_chains = parallel_chains)
    loo_ELPD <- cbind(loo_ELPD, info_model$elpd_pointwise)
    elpd_value <- rbind(elpd_value,info_model$elpd_summary)
    diagn <- rbind(diagn, data.frame(pareto_k = info_model$paretok, 
                                     neff = info_model$neff, 
                                     reff = info_model$reff, 
                                     irt_model = model_specifications[[n]]$irt_model))
  } 
  loo_info   <- list(elpd_value = elpd_value, loo_ELPD = loo_ELPD)
  best_model <- loo_comp(loo_info = loo_info)
  
  info <- list(best_model, diagn)
  return(info)
}



