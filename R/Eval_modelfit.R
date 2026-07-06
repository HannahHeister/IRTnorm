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
#'@param fit  the data output from `fit_IRTnorm()`
#'@param data named list as returned by \code{\link{data_prep}} containing
#'   all data and indexing information required for fitting the IRT norming
#'   model and evaluating its fit.
#'@seealso
#' \code{\link[posterior]{summarise_draws}},
#' \code{\link[posterior]{rhat}},
#' \code{\link[posterior]{ess_bulk}},
#' \code{\link[posterior]{ess_tail}}
#'@return a data.frame summarizing Rhat, Ess bulk and ess tail per variable group
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
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package 'cmdstanr' is required but not installed.\n",
      "Install it with: install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))",
      call. = FALSE
    )
  }
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
    "9"  = c("z_theta", "unique", "nu", "lambda", param1, "alpha", "beta"),
    "10" = c("z_theta", "unique", "nu", "lambda", param1, "alpha", "beta"),
    "11" = c("z_theta", "unique", "unique_sigma", param1, "alpha", "beta"),
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


#' Compute Observed and IRTnorm model-implied Raw Score Distributions
#'
#' @title Posterior Predictive Raw Score Distribution Analysis
#' @name NC
#'
#' @description
#' The function \code{NC} computes both the observed raw score distribution from data
#' and the posterior predictive (IRTnorm model-implied) raw score distributions from a fitted
#' Bayesian IRT-based norming model. The function can perform the analysis either
#' for the entire sample or separately for multiple age groups.
#'
#' @param fit_info a named \code{list} as returned by \code{\link{estimate_norming_model}}
#' @param age_range Numeric \code{vector} defining age group boundaries. If \code{NULL}
#'   (default), no age grouping is performed. If specified, must contain at least
#'   two values defining the boundaries of age groups (e.g., \code{c(0, 30, 60, 90)}
#'   creates three groups: (0,30], (30,60], (60,90]).
#' @return
#' If \code{age_range = NULL}, returns a list containing:
#' \itemize{
#'   \item \code{df}: Data frame with columns:
#'     \itemize{
#'       \item \code{score}: Raw scores from 0 to \code{nitem}
#'       \item \code{obs}: Observed frequency for each score
#'       \item \code{Q5}: 5th percentile of posterior predictive frequencies
#'       \item \code{Q50}: Median (50th percentile) of posterior predictive frequencies
#'       \item \code{Q95}: 95th percentile of posterior predictive frequencies
#'     }
#'   \item \code{replicated}: Matrix of replicated frequencies (draws × scores)
#'   \item \code{nitem}: Number of items
#' }
#'
#' If \code{age_range} is specified, returns a list of lists, one for each age group,
#' where each element contains the same structure as above plus:
#' \itemize{
#'   \item \code{age}: Character string describing the age range (e.g., "age: (30,60]")
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Computes raw scores by summing item responses for each individual
#'   \item Calculates the observed frequency distribution of raw scores
#'   \item Extracts posterior predictive replications from the fitted IRTnorm model
#'   \item Computes the distribution of replicated raw scores for each MCMC draw
#'   \item Summarizes the replicated distributions using 5th, 50th, and 95th percentiles
#' }
#'
#' When \code{age_range} is specified, all steps are performed separately for each
#' age group, allowing assessment of IRTnorm model fit across different demographic subgroups.
#'
#' The \eqn{90\%} credible interval (between Q5 and Q95) provides a range of plausible
#' frequencies under the fitted IRTnorm model. Good IRTnorm model fit is indicated when observed
#' frequencies fall within or near these credible intervals.
#'
#' @note
#' \itemize{
#'   \item Age ranges use the convention (lower, upper], meaning lower is exclusive
#'     and upper is inclusive.
#'   \item Missing scores (scores with zero observations) are included in the output
#'     with frequency 0.
#'   \item The function assumes that \code{fit$y_rep} contains posterior
#'     predictive samples in the correct format.
#' }
#'
#' @examples
#' \dontrun{
#' # Load data
#' data(response_data)
#'
#' # Prepare data for IRTnorm model fit
#' info <- estimate_norming_model(data = response_data,
#' age_variable = "age",
#' int_variables = 1:50,
#' irt_model = "2PLnorm",
#' age__model = "polynom",
#' poly_mean = 3,
#' poly_sd = 2)
#'
#' # Single analysis (no age grouping)
#' result <- NC(fit_info = info)
#'
#' NC_visual(N_info = result)
#'
#' # Age-grouped analysis
#' result_age <- NC(fit_info = info, age_range = seq(-2,2,1))
#' 
#' NC_visual(result_age)
#'}
#'
#' @seealso
#' \code{\link{NC_visual}} for visualizing the results
#'
#' @references
#' Sinharay, S., Johnson, M. S., & Stern, H. S. (2006). Posterior predictive assessment
#' of item response theory models. Applied Psychological Measurement, 30(4), 298-321.
#'
#' @export
NC <- function(fit_info,
               age_range = NULL) {

  data <- fit_info$newdata
  model_info <- fit_info$model_info
  raw_data <- fit_info$raw_data
  int_variables <- fit_info$int_variables
  age_variable <- fit_info$age_variable
  
  # Validate inputs
  if (!is.null(age_range) && is.null(raw_data)) {
    stop("'raw_data' must be provided when 'age_range' is specified")
  }

  int_data <- raw_data[, int_variables, drop = FALSE]
  
  if(all(data$itemD == 1)){
    info_all <-  vector("list", 1)
  }
  if(any(data$itemD > 1)){
    info_all <-  vector("list", data$D +1)
  }
if(data$irt_mod !=  4){
  nitem <- length(int_variables)
}else{
  nitem <- sum(data$K_item)
}

  # Compute replicated raw score distributions (shared computation)
  score_rep <- apply(model_info$y_rep, 1, function(x) tapply(x, data$jj, sum))
  draws <- ncol(score_rep)

  # Compute observed scores
  observed_scores <- rowSums(int_data)

  if (is.null(age_range)) {
    info_all[[1]] <- compute_nc_overall(observed_scores, score_rep, nitem, draws)
  } else {
    info_all[[1]] <- compute_nc_by_age(observed_scores, raw_data, age_variable,
                                       age_range, score_rep, nitem, draws)
  }
  if(any(data$itemD > 1)){
    for(d in seq_len(data$D)){
      nitem <- sum(data$itemD == d)
      # Compute replicated raw score distributions (shared computation)
      score_rep <- apply(model_info$y_rep[,data$dd == d], 1, function(x) tapply(x, data$jj[data$dd == d], sum))
      draws <- ncol(score_rep)

      # Compute observed scores
      observed_scores <- rowSums(int_data[,data$itemD == d ])
      if (is.null(age_range)) {
        info_all[[d +1]] <- compute_nc_overall(observed_scores, score_rep, nitem, draws, dim = d)
      } else {
        info_all[[d + 1]] <- compute_nc_by_age(observed_scores, raw_data, age_variable,
                                               age_range, score_rep, nitem, draws, dim = d)
      }
    }
  }

  return(info_all)
}


# Helper: Create full observed distribution with zeros for missing scores
create_full_distribution <- function(score_obs, nitem) {
  if (length(score_obs) == 0) {
    return(data.frame(score = 0:nitem, obs = 0))
  }
  observed <- as.data.frame(table(score_obs))
  colnames(observed) <- c("score", "obs")
  observed$score <- as.numeric(as.character(observed$score))

  # Fill in missing scores with zero counts
  if (nrow(observed) != nitem + 1) {
    full_dist <- data.frame(score = 0:nitem, obs = 0)
    full_dist$obs[match(observed$score, full_dist$score)] <- observed$obs
    return(full_dist)
  }

  return(observed)
}

# Helper: Compute replicated score distributions
compute_replicated_distributions <- function(score_rep, age_filter = NULL,
                                             nitem, draws) {
  replicated <- matrix(nrow = draws, ncol = nitem + 1)

  for (i in 1:draws) {
    score_data <- if (is.null(age_filter)) {
      score_rep[, i]
    } else {
      score_rep[age_filter, i]
    }

    # Count occurrences of each score (0 to nitem)
    replicated[i, ] <- tabulate(score_data + 1, nbins = nitem + 1)
  }

  return(replicated)
}

# Helper: Compute quantiles from replicated distributions
compute_quantiles <- function(replicated) {
  quantiles <- t(apply(replicated, 2, stats::quantile, probs = c(0.05, 0.5, 0.95)))
  quantiles <- as.data.frame(quantiles)
  colnames(quantiles) <- c("Q5", "Q50", "Q95")
  return(quantiles)
}

# Helper: Process a single group (overall or age-specific)
process_group <- function(score_obs, score_rep, nitem, draws, age_filter = NULL) {
  # Observed distribution
  observed_full <- create_full_distribution(score_obs, nitem)

  # Replicated distributions
  replicated <- compute_replicated_distributions(score_rep, age_filter, nitem, draws)

  # Quantiles
  rep_quantiles <- compute_quantiles(replicated)

  # Combine results
  df <- cbind(observed_full, rep_quantiles)

  return(list(
    df = df,
    replicated = replicated,
    nitem = nitem
  ))
}

# Compute NC for overall sample
compute_nc_overall <- function(observed_scores, score_rep, nitem, draws, dim = 0) {
  result <- process_group(observed_scores, score_rep, nitem, draws)
  result$age <- "full age range"
  result$dim <- dim
  return(result)
}

# Compute NC by age groups
compute_nc_by_age <- function(observed_scores, raw_data, age_variable,
                              age_range, score_rep, nitem, draws, dim = 0) {
  age_groups <- length(age_range) - 1
  results <- vector("list", age_groups)

  ages <- raw_data[, age_variable]

  for (i in 1:age_groups) {
    # Define age filter
    age_filter <- ages > age_range[i] & ages <= age_range[i + 1]

    # Process this age group
    result <- process_group(
      score_obs = observed_scores[age_filter],
      score_rep = score_rep,
      nitem = nitem,
      draws = draws,
      age_filter = age_filter
    )

    # Add age label
    result$age <- sprintf("age: (%s,%s]", age_range[i], age_range[i + 1])
    result$dim <- dim
    results[[i]] <- result
  }

  return(results)
}

#' Visualize IRTnorm model-Implied vs Observed Raw Score Distributions
#'
#' @title Posterior Predictive IRTnorm model Check for Raw Score Distributions
#' @name NC_visual
#' @description
#' The function \code{NC_visual} creates a visualization comparing the IRTnorm model-implied
#' (posterior predictive) raw score distribution against the observed raw score
#' distribution. This serves as a posterior predictive IRTnorm model check (PPMC) to assess
#' how well the fitted IRTnorm model reproduces the observed data patterns.
#'
#' The function automatically detects whether the input comes from a single analysis
#' or from age-grouped analyses, and creates the appropriate visualization (single
#' plot or faceted plot by age groups).
#'
#' @param NC_info A list object returned by the \code{NC} function. This can be either:
#'   \itemize{
#'     \item A single result list (when \code{age_range = NULL} in \code{NC})
#'     \item A list of results by age groups (when \code{age_range} is specified in \code{NC})
#'   }
#'   Each result contains:
#'   \itemize{
#'     \item \code{df}: Data frame with observed frequencies and quantiles
#'     \item \code{replicated}: Matrix of replicated raw score distributions
#'     \item \code{nitem}: Number of items
#'     \item \code{age}: Age group label (only for age-grouped results)
#'   }
#'@param dim_plot \code{character} indicating which latent trait of the HO-2PL-norm
#'IRTnorm model should be visualized.
#'Can only be of the following characters: "general", "domain1"...
#'
#' @return Invisibly returns a ggplot2 object. The plot is printed to the current
#'   graphics device. The returned object can be further customized or saved.
#'
#' @details
#' The visualization includes:
#' \itemize{
#'   \item \strong{Orange jittered points}: Individual posterior predictive draws
#'     showing the distribution of replicated frequencies across MCMC iterations
#'   \item \strong{Dashed line}: Median (50th percentile) of the replicated distributions
#'   \item \strong{Dotted lines}: 5th and 95th percentiles of the replicated distributions,
#'     forming a \eqn{90\%}credible interval
#'   \item \strong{Black points}: Observed raw score frequencies from the actual data
#' }
#'
#' Good IRTnorm model fit is indicated when observed frequencies (black points) fall within
#' or near the credible intervals (between dotted lines) and follow the general pattern
#' of the median line.
#'
#' For age-grouped analyses, the function creates a faceted plot with separate panels
#' for each age group, allowing comparison of IRTnorm model fit across different age ranges.
#' @examples
#' \dontrun{
#' # Load data
#' data(response_data)
#'
#' # Prepare data for IRTnorm model fit
#' info <- estimate_norming_model(data = response_data,
#' age_variable = "age",
#' int_variables = 1:50,
#' irt_model = "2PLnorm",
#' age__model = "polynom",
#' poly_mean = 3,
#' poly_sd = 2)
#' 
#' # Single analysis (no age grouping)
#' result <- NC(fit_info = info_poly)
#'
#' NC_visual(N_info = result)
#'
#' # Age-grouped analysis
#' result_age <- NC(fit_info = info,
#'                  age_range = seq(-2,2,1))
#'
#' NC_visual(result_age)
#'}
#' @seealso \code{\link{NC}} for computing the raw score distributions
#'
#' @references
#' Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B. (2013).
#' Bayesian Data Analysis (3rd ed.). Chapman and Hall/CRC.
#'
#' @export
NC_visual <- function(NC_info, dim_plot = NULL) {
##TODO NEEDS TO BE CHECKED IF THIS IF CASE IS GOOD 
  if (length(NC_info) == 1 && is.list(NC_info[[1]]) && 
      !is.null(NC_info[[1]][[1]]$age)) {
    NC_info <- NC_info[[1]]
  }
  
  
process_single_result <- function(result, age_label = NULL, dim_label = NULL) {
  # Extract data
  df <- result$df
  replicated2 <- data.frame(result$replicated)
  colnames(replicated2) <- seq(0, result$nitem, 1)
  replicated2$id <- 1:nrow(replicated2)

  # Reshape to long format
  rep.long <- replicated2 |>
    tidyr::pivot_longer(cols = -id, names_to = "score", values_to = "frequency")
  rep.long$score <- as.numeric(rep.long$score)

  # Add age group label if provided
  if (!is.null(age_label)) {
    df$age_group <- age_label
    rep.long$age_group <- age_label
  }
  # Add domain label if provided
  if (!is.null(dim_label)) {
    df$domian <- dim_label
    rep.long$domian <- dim_label
  }
  return(list(df = df, rep.long = rep.long))
}


# Helper
remove_substr <- function(text, start, end) {
  paste0(substring(text, 1, start - 1),
         substring(text, end + 1))
}
# True if multiple dimensions are tested and if multiple ages are evaluated
if(length(NC_info) > 1 & !is.null(NC_info[[1]][[1]]$age)){
  if(is.null(dim_plot)){
    stop("'dim_plot' missing. Please provide which latent trait dimension  should be visualized. Due to
        the different age groups the latent trait dimensions have to be plotted individually.")
  }else{
    if(dim_plot == "general"){
      NC_info <- NC_info[[1]]
    }else{
      int_NC <- as.numeric(remove_substr(dim_plot, 1,6))
      if(!is.na(int_NC)){
        NC_info <- NC_info[[as.numeric(substr(dim_plot, 7,7))]] # +1
      }else{
        stop("'dim_plot' has the wrong format. ")
      }

    }
  }
}

# Detect if input is a list of lists or a single result
n_group <- length(NC_info)

if (n_group > 1) {
  # Process multiple  groups
  df_all <- NULL
  rep.long_all <- NULL

  for (g in 1:n_group) {
    processed <- process_single_result(result = NC_info[[g]],
                                       age_label =NC_info[[g]]$age,
                                       dim_label =NC_info[[g]]$dim)
    df_all <- rbind(df_all, processed$df)
    rep.long_all <- rbind(rep.long_all, processed$rep.long)

  }
}else {
  # Process single result (no age grouping)
  processed <- process_single_result(NC_info[[1]], age_label = "one", dim_label = "none")
  df_all <- processed$df
  rep.long_all <- processed$rep.long

}
n_domains <- length(unique(rep.long_all$domian)) > 1
if(n_domains){
  rep.long_all$domian <- factor(rep.long_all$domian,
                                labels = c("general", paste0("domain", seq_len(max(rep.long_all$domian)))))
  df_all$domian <- factor(df_all$domian,
                          labels = c("general", paste0("domain", seq_len(max(df_all$domian)))))
}
# Create plot
p <- ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::geom_jitter(
    data = rep.long_all,
    ggplot2::aes(x = score, y = frequency),
    color = "#fdb863",
    position = ggplot2::position_jitter(width = 0.12),
    alpha = 0.15
  ) +
  ggplot2::geom_line(
    data = df_all,
    ggplot2::aes(x = score, y = Q50), #+1 and other quantiles
    linewidth = 1,
    linetype = "dashed"
  ) +
  ggplot2::geom_line(
    data = df_all,
    ggplot2::aes(x = score , y = Q5),
    linewidth = 1,
    linetype = "dotted"
  ) +
  ggplot2::geom_line(
    data = df_all,
    ggplot2::aes(x = score , y = Q95),
    linewidth = 1,
    linetype = "dotted"
  ) +
  ggplot2::geom_point(
    data = df_all,
    ggplot2::aes(x = score , y = obs),
    size = 4
  ) +
  ggplot2::xlab("raw score") +
  ggplot2::ylab("frequency of individuals") +
  ggplot2::theme(
    text = ggplot2::element_text(size = 13, family = "serif"),
    axis.line = ggplot2::element_line(color = 'black'),
    plot.background = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    legend.position = "none"
  )

n_age <- length(unique(rep.long_all$age_group))
if(n_age > 1){
  p <- p + ggplot2::facet_wrap(~ age_group, scale = "free", ncol = round(sqrt(n_age)))
}
if(n_domains){
  n_domain <- length(unique(rep.long_all$domian))
  p <-  p + ggplot2::facet_wrap(~ domian, scale = "free", ncol = round(sqrt(n_domain)))
}
if(!is.null(dim_plot)){
  p <- p + ggplot2::ggtitle(dim_plot)
}
print(p)
return(invisible(p))
}

#' Visualize the distribution of norm scores across age ranges
#'
#' The \code{normscore_dist} function visualizes whether estimated norm scores
#' follow an approximately standard normal distribution across different age
#' ranges. For each specified age group, a histogram of estimated norm scores is
#' plotted together with the corresponding density of a standard normal
#' distribution.
#'
#' This diagnostic plot is useful for assessing whether the norming procedure
#' achieved the intended normalization of the latent trait within age groups.
#'
#' @title Check the distribution of norm scores across age ranges
#'
#' @name normscore_dist
#'
#' @param fit_info a named \code{list} as returned by \code{\link{estimate_norming_model}}
#'
#' @param age_range A \code{numeric} vector specifying the cut-off points for
#'   age groups. Consecutive values define the lower (exclusive) and upper
#'   (inclusive) bounds of each age group.
#'
#' @param binwidth A \code{numeric} value specifying the bin width of the
#'   histogram. Defaults to \code{0.5}.
#'
#' @param dim_plot only required for multidimensional IRTnorm models.
#' \code{character} indicating which latent trait of the higher-order model 
#' should be visualized.
#'Can only be of the following characters: "general", "domain1"...
#' @return A \pkg{ggplot2} object showing histograms of norm scores by age group,
#'   overlaid with the scaled density of a standard normal distribution.
#'
#' @details
#' Individuals are assigned to age groups defined by \code{age_range} using the
#' intervals \eqn{(a_i, a_{i+1}]}. For each age group, the histogram of estimated
#' norm scores (\code{normscore_est}) is plotted together with the expected
#' density of a standard normal distribution, scaled by the number of
#' observations and the histogram bin width.
#'
#' Faceting is used to display the distributions separately for each age group.
#'
#' @seealso \code{\link{extract_info}}, \code{\link{fit_IRTnorm}}
#' @examples
#' \dontrun{
#' # Load data
#' data(response_data)
#'
#' # Fit IRTnorm model and extract information
#' info <- estimate_norming_model(data = response_data,
#' age_variable = "age",
#' int_variables = 1:50,
#' irt_model = "2PLnorm",
#' age__model = "polynom",
#' poly_mean = 3,
#' poly_sd = 2)
#'
#' # Plot the density of the normed latent trait in age-groups
#' normscore_dist(fit_info = info, age_range = seq(-2,2,1))}
#' @export
normscore_dist <- function(fit_info, age_range, binwidth = 0.5, dim_plot = NULL){
  
  
  model_info <- fit_info$model_info 
  if(model_info$info_sample$irt_mod[1] <= 4){
    data <-  model_info$info_sample
  }
  if(model_info$info_sample$irt_mod[1] > 4){
    if(is.null(dim_plot)){
      stop("Your data is from a multidimensional IRT model. Choose which dimension
           you want to visualize by specifying 'dim_plot'.")
    }
    if(dim_plot == "general"){
      data <-  model_info$info_sample[model_info$info_sample$trait == "theta", ]
    }else{
      data <-  model_info$info_sample[model_info$info_sample$trait == paste0("normeddelta",substr(dim_plot, 7,7)), ]
    }
    data$normscore_est <- data$est_ability
  }
  if (min(data$normscore_est) == max(data$normscore_est)) {
    stop("All norm scores are identical. Cannot plot a meaningful distribution.")
  }
  
  data$age_group <- NA
  n_agegroup <- (length(age_range) - 1)
  n_group <-NULL
  for(age in 1:n_agegroup){
    data$age_group[data$age > age_range[age] & data$age <= age_range[age +1]] <- paste0("age: (", age_range[age], ",", age_range[age +1], "]")
    n_group <- c(n_group,sum(data$age > age_range[age] & data$age <= age_range[age +1]) )
  }
  if(any(is.na(data$age_group))){
    warning("Some individuals are not within the plotted age groups. These have been removed")
    data <- data[!is.na(data$age_group), ]
  }

  
  x_val <- seq(min(data$normscore_est), max(data$normscore_est), length.out = 300)  # standard x values

  curve_list <- lapply(seq_along(n_group), function(i) {
    data.frame(
      x = x_val,
      y = stats::dnorm(x_val, 0, 1) * binwidth * n_group[i],
      age_group = paste0("age: (", age_range[i], ",", age_range[i + 1], "]")
    )
  })

  curve_data <- do.call(rbind, curve_list)
  
  ggplot2::ggplot()  +
    ggplot2::geom_histogram(ggplot2::aes(x = normscore_est), binwidth = binwidth ,
                            data = data, color = "#fdb863", fill = "#fdb863") +
    ggplot2::geom_line(data = curve_data, ggplot2::aes(x, y), linewidth = 1) +
    ggplot2::facet_grid(~age_group) + ggplot2::xlab("norm score") +
    ggplot2::ylab("frequency of individuals")


}



# shared plotting theme, defined once
irt_plot_theme <- function(){
  ggplot2::theme(
    text             = ggplot2::element_text(size = 13, family = "serif"),
    legend.position  = "bottom",
    axis.line        = ggplot2::element_line(color = "black"),
    plot.background  = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    legend.title     = ggplot2::element_blank()
  )
}

# is this a multidimensional IRT model? (threshold defined once)
is_multidim <- function(irt_mod) irt_mod >= 5


#' Visualize observed data and IRT norming model–implied quantile curves
#'
#' \code{plot_quantile_curves} visualizes observed data across age together
#' with IRTnorm model-implied quantile curves derived from an age-dependent IRT
#' norming model. Depending on \code{perspective}, the function either compares
#' observed raw sum scores to model-implied raw score quantiles
#' (\code{perspective = "rawscore"}), or compares realized latent trait
#' estimates to model-implied latent trait quantiles
#' (\code{perspective = "ability"}).
#'
#' For the raw score perspective, quantile curves are computed from posterior
#' predictive replicated responses and smoothed as a function of age. It is
#' important to note that, based on the IRT-based norming model, the raw score
#' does not provide all information on the estimated normed latent trait.
#' Therefore the plot can only show the aggregated information of the model,
#' which does not completely coincide with the model-implied quantiles, as
#' identical raw scores can lead to different quantile positions if the
#' answer pattern differs.
#'
#' For the ability perspective, quantiles of the model-implied latent trait
#' distribution are computed for each individual and smoothed as a function
#' of age using local polynomial regression (LOESS), while realized latent
#' trait estimates (\code{theta_est}) are displayed as points.
#'
#' Both perspectives provide a diagnostic comparison between observed data
#' and the distribution implied by the fitted IRT norming model, allowing
#' assessment of whether the model reproduces observed patterns
#' (raw score or latent trait) across the entire age range.
#'
#' @title Visualize IRT norming model–implied quantiles (raw score or ability)
#' @name quantile_curves
#'
#' @param fit_info a named \code{list} as returned by \code{\link{estimate_norming_model}}
#'
#' @param perspective character string specifying which diagnostic to plot.
#'   One of \code{"rawscore"} or \code{"ability"}. With \code{"rawscore"}, the
#'   function visualizes observed raw sum scores against model-implied raw
#'   score quantiles across age. With \code{"ability"}, the function visualizes
#'   the internal consistency of the model by comparing realized latent trait
#'   estimates with model-implied quantiles across age. Default is "ability". 
#'@param probs vector indicating which quantiles should be plotted. 
#' Default is c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95). 
#' 
#' @return A \pkg{ggplot2} object. For \code{perspective = "rawscore"}, observed
#'   raw sum scores plotted against age together with smoothed model-implied
#'   quantile curves (\eqn{5\%}, \eqn{25\%}, \eqn{50\%}, \eqn{75\%}, \eqn{95\%}).
#'   For \code{perspective = "ability"}, realized latent trait estimates
#'   plotted against age together with smoothed model-implied latent trait
#'   quantile curves (\eqn{5\%}, \eqn{25\%}, \eqn{50\%}, \eqn{75\%}, \eqn{95\%}).
#'
#' @details
#' \strong{Raw score perspective:} Observed raw sum scores are computed as row
#' sums across the item response variables specified in \code{int_variables}.
#' Model-implied raw score distributions are obtained by summing posterior
#' predictive replicated item responses for each individual. Empirical
#' quantiles of the replicated raw scores are then smoothed as a function of
#' age using LOESS.
#'
#' \strong{Ability perspective:} For each individual, quantiles of the
#' model-implied latent trait distribution are computed from the
#' age-dependent normal distribution of the latent trait, and smoothed as a
#' function of age using LOESS. Realized latent trait estimates are overlaid
#' as points for comparison.
#'
#' In both cases, the resulting plot serves as a diagnostic tool for
#' evaluating whether the fitted IRT norming model reproduces observed
#' patterns across the entire age range.
#'
#' @examples
#' \dontrun{
#' # Load data
#' data(response_data)
#'
#' info <- estimate_norming_model(
#'   raw_data        = response_data,
#'   age_variable    = "age",
#'   int_variables   = 1:50,
#'   irt_model       = "2PLnorm",
#'   poly_mean       = 3,
#'   poly_sd         = 2,
#'   age_model       = "poly",
#'   iter_warmup     = 500,
#'   iter_sampling   = 1000,
#'   chains          = 2,
#'   parallel_chains = 2
#' )
#'
#' # Ability perspective
#' plot_quantile_curves(fit_info = info, perspective = "ability")
#'
#' # Raw score perspective
#' plot_quantile_curves(fit_info = info, perspective = "rawscore")
#' }
#' @export
quantile_curves <- function(fit_info, perspective = c("ability", "rawscore"), 
                            probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)){
  
  perspective <- match.arg(perspective)
  if(any(probs >= 1| probs <= 0)){
    stop("probs have to be between in (0,1)")
  }
  
  model_info  <- fit_info$model_info
  
  if (perspective == "rawscore"){
    data         <- fit_info$newdata
    raw_data     <- fit_info$raw_data
    age_variable <- fit_info$age_variable
    int_variable <- fit_info$int_variable 
    
    answer_pattern <- raw_data[, int_variable, drop = FALSE]
    
    rep_quantile_all <- lapply(unique(data$itemD), function(d){
      
      response_agg <- data.frame(y = Matrix::rowSums(answer_pattern[,data$itemD == d, drop = FALSE], na.rm = TRUE),
                                 age = raw_data[, age_variable])
      y_rep <- model_info$y_rep[, data$dd == d]
      person_jj <- data$jj[data$dd == d]
      score_rep <- apply(y_rep, 1, function(x) tapply(x, person_jj, sum))
      # apply drops dimensions when result is length 1 — ensure always a matrix
      if (!is.matrix(score_rep)) {
        score_rep <- matrix(score_rep, nrow = length(unique(person_jj)))
      }
      nperson <- nrow(score_rep)
      
      rep_quantile <- t(apply(score_rep, 1, stats::quantile, probs = probs))
      
      
      data.frame(rep       = as.vector(rep_quantile),
                 age       = rep(response_agg$age, length(probs)),
                 observed  = rep(response_agg$y, length(probs)),
                 quantiles = rep(as.character(probs), each = nperson),
                 domain    = paste0("domain ", d))
    })
    rep_quantile_all <- do.call(rbind, rep_quantile_all)
    
    p <- ggplot2::ggplot(rep_quantile_all) +
      ggplot2::geom_smooth(ggplot2::aes(age, rep, color = quantiles), formula = 'y ~ x', method = 'loess') +
      ggplot2::geom_point(ggplot2::aes(age, observed)) +
      ggplot2::ylab("sum score") + ggplot2::theme_bw() + irt_plot_theme()
    
    if (is_multidim(data$irt_mod)) p <- p + ggplot2::facet_grid(~domain)
    print(p)
    return(invisible(rep_quantile_all))
  }
  
  # perspective == "ability"
  irt_mod <- model_info$info_sample$irt_mod[1]
  
  if (!is_multidim(irt_mod)){
    person_data <- model_info$info_sample
    nperson     <- nrow(person_data)
    
    quan <- as.vector(mapply(
      function(mu, sigma) stats::qnorm(probs, mean = mu, sd = sigma),
      person_data$theta_mu, person_data$theta_sigma
    ))
    
    new_data <- data.frame(
      age   = rep(person_data$age, each = length(probs)),
      quan  = quan,
      qname = rep(paste0("p", probs), nperson)
    )
  } 
  else {
    person_data <- model_info$info_sample[substr(model_info$info_sample$trait, 1, 1) == "d", ]
    dims        <- unique(person_data$trait)
    nperson     <- fit_info$newdata$J
    
    new_data <- do.call(rbind, lapply(dims, function(d){
      int_data <- person_data[person_data$trait == d, ]
      quan <- as.vector(mapply(
        function(mu, sigma) stats::qnorm(probs, mean = mu, sd = sigma),
        int_data$theta_mu, int_data$theta_sigma
      ))
      data.frame(
        age    = rep(int_data$age, each = length(probs)),
        quan   = quan,
        qname  = rep(paste0("p", probs), nperson),
        domain = d
      )
    }))
    
    person_data$theta_est <- person_data$est_ability
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(age, theta_est), data = person_data) +
    ggplot2::geom_smooth(ggplot2::aes(age, quan, color = qname), data = new_data, formula = 'y ~ x', method = 'loess') +
    ggplot2::ylab("latent trait") + ggplot2::theme_bw() + irt_plot_theme()
  
  if (is_multidim(irt_mod)) p <- p + ggplot2::facet_grid(~domain)
  print(p)
  invisible(new_data)
}

#' Sim responses from a Graded Response Model
#'
#' Computes the category response probabilities for a single item under the
#' Graded Response (GRM).
#'
#' @param theta A \code{numeric} scalar or vector of latent trait values.
#' @param a A \code{numeric} scalar giving the item discrimination parameter.
#' @param b A \code{numeric} vector of item threshold parameters of length
#'   \eqn{m - 1} where \eqn{m} is the number of response categories.
#'
#' @return A \code{numeric} vector or \code{matrix} of category probabilities.
#'
#' @noRd
grm_probs <- function(theta, a,b){

  cum_prob <- stats::plogis(a * (theta- b))
if(length(theta) == 1){
  p <- c(1, cum_prob) - c(cum_prob, 0)
}else{
  p <- cbind(1, cum_prob) - cbind(cum_prob, 0)

}
  return(p)
}
#' Create a dummy coded response matrix
#'
#' Converts a matrix of polytomous item responses into a dummy coded matrix
#' where each response category of each item gets its own binary column.
#'
#' @param info A \code{matrix} or \code{data.frame} of item responses with
#'   persons in rows and items in columns. Each cell contains the response
#'   category for that person-item combination.
#' @param K_item A \code{numeric} vector of length equal to the number of items,
#'   giving the number of response categories for each item.
#'
#' @return A binary \code{matrix} of dimensions \eqn{N \times \sum K_i} where
#'   \eqn{N} is the number of persons and \eqn{\sum K_i} is the total number of
#'   response categories across all items. Column names follow the pattern
#'   \code{ItemI_CatK}. Each column contains 1 if the person responded in that
#'   category and 0 otherwise.
#'
#' @noRd
make_dummy_matrix <- function(info, K_item) {

  N <- nrow(info)          # persons
  I <- ncol(info)          # items
  total_K <- sum(K_item)

  # initialize matrix
  dummy <- matrix(0, nrow = N, ncol = total_K)

  # optional: column names
  col_names <- c()

  col_index <- 1

  for(i in 1:I){

    Ki <- K_item[i]
    responses <- info[, i]

    for(k in 1:Ki){
      dummy[, col_index] <- as.integer(responses == k)
      col_names <- c(col_names, paste0("Item", i, "_Cat", k))
      col_index <- col_index + 1
    }
  }

  colnames(dummy) <- col_names

  return(dummy)
}



#'@title Visualization of item and person fit
#'@name item_person_fit
#'@description Visualizes item-level or person-level fit by comparing observed
#'response patterns with IRTnorm model-implied response patterns.
#'
#' @param fit_info a named \code{list} as returned by \code{\link{estimate_norming_model}}
#'@param parameter a \code{character} indicating whether item fit or person fit
#' should be visualized. Must be either \code{"item"} or \code{"person"}.
#' @return A \pkg{ggplot2} object showing smoothed quantile curves of the
#'   IRTnorm model-implied latent trait distribution (\eqn{5\%},\eqn{ 25\%},\eqn{ 50\%},
#'   \eqn{75\%}, and \eqn{95\%}). together with realized latent trait estimates plotted against age.
#'@details For each person–item combination, the function computes the IRTnorm model-implied
#'probability of a correct response. Depending on \code{estimate}, responses are
#' either dichotomized at 0.5 (\code{"mostlikely"}) or retained as expected probabilities (\code{"expected"}).
#'
#' If \code{parameter = "item"}, the plot compares observed and the IRTnorm model-implied
#' the amount of correct answers per item.
#'
#' If \code{parameter = "person"}, the plot compares observed and the IRTnorm model-implied
#' raw score per person.
#'
#' A diagnonal reference line is added to facilitate visual assessment of fit.
#' The closer the observed points to the diagonal the better the fit.
#'@examples
#'\dontrun{
#' # Load data
#' data(response_data)
#'
#' # Fit IRTnorm model and extract information
#' info <- estimate_norming_model(data = response_data,
#' age_variable = "age",
#' int_variables = 1:50,
#' irt_model = "2PLnorm",
#' age__model = "polynom",
#' poly_mean = 3,
#' poly_sd = 2)
#' 
#'  # Example usage
#'  item_person_fit(
#'    fit_info = info,
#'    parameter = "item"
#'  )
#'
#'  item_person_fit(
#'    fit_info = info,
#'    parameter = "person"
#'  )
#'}
#'
#' @seealso
#' \code{\link{extract_info}}
#'@export
item_person_fit <- function(fit_info,  parameter){

model_info <- fit_info$model_info
answer_pattern <- fit_info$raw_data[, fit_info$int_variables]
data <- fit_info$newdata
  nperson <- data$J
  nitem <- data$I
  if(nrow(answer_pattern) != nperson){
    stop("'data$J' (", nperson, ") does not match the number of rows in 'answer_pattern' (",
         nrow(answer_pattern), "). Check that data$J reflects the actual number of persons.")
  }
  if(ncol(answer_pattern) != nitem){
    stop("'data$I' (", nitem, ") does not match the number of columns in 'answer_pattern' (",
         ncol(answer_pattern), "). Check that data$I reflects the actual number of items.")
  }
  
  
  if(data$irt_mod <= 4){
    info_person <- matrix(model_info$info_sample$theta_est, ncol = 1)
  }
  if(data$irt_mod >= 5 ){
    if(any(!data$itemD %in% seq_len(data$D))){
      stop("'data$itemD' contains dimension indices that do not match the number of dimensions in 'data$D'. Check that itemD values run from 1 to D.")
    }
    temp <- model_info$info_sample[substr(model_info$info_sample$trait,1,1) == "d",]
    temp$id <- rep(seq_len(nperson), data$D)
    info_person <- matrix(temp$est_ability, ncol = data$D)
  }
  info_item   <- model_info$info_item

 
  
  # Take the most likely answer pattern
  if(data$irt_mod != 4){
    expected <- matrix(
      0,
      nrow =nperson,
      ncol = nitem
    )
    for(i in seq_len(nitem)){
      expected[, i] <- info_item$gamma[i] + (1 - info_item$gamma[i])*(exp(info_item$alpha[i] * (info_person[,data$itemD[i]] - info_item$beta[i])) /
                                                                        (1 + exp(info_item$alpha[i] * (info_person[,data$itemD[i]] - info_item$beta[i]))))
    }
    answer_pattern_new <- answer_pattern
  }else{
    expected <- matrix(
      0,
      nrow =nperson,
      ncol = nitem*(max(data$K_item))
    )
    int_threshold <- seq(1,nitem*data$K_max, nitem)[1:(data$K_item[1] -1)]
    expected[, 1:(sum(data$K_item[1]))] <- grm_probs(theta = info_person, a = info_item$alpha[1], b = info_item$threshold[int_threshold])
    for(i in seq_len(nitem)[-1]){
      int_threshold <- seq(i,nitem*data$K_max, nitem)[1:(data$K_item[i] -1)]
      expected[, (1+ sum(data$K_item[1:(i-1)])):(sum(data$K_item[1:i]))] <- grm_probs(theta = info_person, a = info_item$alpha[i], b = info_item$threshold[int_threshold])
    }
    answer_pattern_new <- make_dummy_matrix(info = answer_pattern, K_item = data$K_item)
  }
  if(any(is.na(answer_pattern_new))){
    expected[is.na(answer_pattern_new)] <- NA
  }
  if(parameter == "item"){
    
    if(data$irt_mod <= 3 ){
      df <- data.frame(
        observed = colSums(answer_pattern_new, na.rm = TRUE) /colSums(!is.na(expected)),
        irt_mod    = colSums(expected, na.rm = TRUE) / colSums(!is.na(expected))
      )
    }
    if(data$irt_mod == 4){
      df <- data.frame(value = c(colSums(answer_pattern_new, na.rm = TRUE) /colSums(!is.na(expected)),
                                 colSums(expected, na.rm = TRUE) / colSums(!is.na(expected))),
                       type = rep(c("observed","irt_mod"), each = sum(data$K_item)),
                       item = rep(rep(1:nitem, data$K_item), 2),
                       level = rep(sequence(data$K_item),2))
    }
    if(data$irt_mod > 4){
      df <- data.frame(
        observed = colSums(answer_pattern, na.rm = TRUE) /colSums(!is.na(expected)),
        irt_mod    = colSums(expected, na.rm = TRUE) / colSums(!is.na(expected)),
        domain = paste0("domain ",data$itemD))
    }
    if(any(is.nan(df$observed)) || any(is.nan(df$irt_mod))){
      warning("One or more items have no valid responses. These items are excluded from the plot.")
      df <- df[!is.nan(df$observed) & !is.nan(df$irt_mod), ]
    }

    xname <- "Proportion correct per item observed"
    yname <- "Proportion correct per item based on model"
    titlename <- "Item fit"
  }

  if(parameter == "person"){
    if(data$irt_mod <= 4){
if(data$irt_mod <=3){
  data$K_item <- rep(1, nitem)
}
      df <- data.frame(
        observed = rowSums(answer_pattern, na.rm = TRUE),
        irt_mod    = rowSums(expected *matrix(rep(sequence(data$K_item), nperson), nrow = nperson, byrow = TRUE
        ), na.rm = TRUE)
      )
    }else{
      df <- NULL
      for(d in seq_len(data$D)){
        df <- rbind(df, data.frame(
          observed = rowSums(answer_pattern[data$itemD == d, ], na.rm = TRUE),
          irt_mod    = rowSums(expected[data$itemD == d, ], na.rm = TRUE),
          domain = paste0("domain ", d)))
      }
    }
    xname <- "Realized raw score"
    yname <- "model-implied raw score"
    titlename <- "Person fit"
  }
  if(data$irt_mod == 4 & parameter == "item"){
    p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(type), y = value, fill = factor(level))) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::facet_wrap(~ paste0("item ", item)) +
      #scale_fill_viridis_d(option = "plasma") +
      ggplot2::labs(x = NULL, y = "Proportion", fill = "Level") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  }else{
    p <- ggplot2::ggplot(df, ggplot2::aes(x = observed, y = irt_mod)) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(
        intercept = 0, slope = 1,
        linewidth = 1,
        color = "#CC4678FF"
      ) +   ggplot2::labs(
        x = xname,
        y = yname,
        title =titlename
      ) +
      ggplot2::theme_bw()
  }


  if(data$irt_mod > 4){
    p <- p +  ggplot2::facet_grid(~domain)

  }
  return(p)
  invisible(df)
}




#' @title Plot the estimated age effect of IRT-based continuous norming model
#'@name effect_age
#'@description
#' Visualizes the IRTnorm model-implied relationship between age and the
#' latent trait distribution mean and standard deviation.
#'
#'@param fit_info a list containing IRTnorm model output from the \code{extract_info()} function.
#'It must include the elements \code{info_sample} (person-level information,
#' including latent trait estimates) and \code{info_item} (item parameters).
#'@return A \pkg{ggplot2} object showing the estimated mean and standard deviation of the
#' latent trait as functions of age.
#' @seealso
#' \code{\link{extract_info}}
#'@examples
#'\dontrun{
#' # Load data
#' data(response_data)
#'
#' # Fit IRTnorm model and extract information
#' info <- estimate_norming_model(data = response_data,
#' age_variable = "age",
#' int_variables = 1:50,
#' irt_model = "2PLnorm",
#' age__model = "polynom",
#' poly_mean = 3,
#' poly_sd = 2)
#'
#' # Plot estimated age effect
#' effect_age(fit_info = info)
#'}
#'@export
effect_age <- function(fit_info){
  info_person <- fit_info$model_info$info_sample
  
  if(info_person$irt_mod[1] <= 4){
    df_long <- data.frame(
      age = rep(info_person$age, 2),
      value = c(info_person$theta_mu, info_person$theta_sigma),
      moment = factor(
        rep(c("mean", "standard deviation"), each = nrow(info_person)),
        levels = c("mean", "standard deviation")
      )
    )
  }
  if(info_person$irt_mod[1] > 4){
    if(!any(substr(info_person$trait, 1, 1) == "d")){
      stop("No traits starting with 'd' found in 'fit_info$info_sample$trait'. ",
           "Expected trait names like 'delta1', 'delta2' for multidimensional irt_mods. ",
           "Check that 'fit_info$info_sample' contains the correct trait labels.")
    }
    data_temp <- info_person[substr(info_person$trait,1,1) == "d", ]
    D <- length(unique(data_temp$trait))
    nperson <- nrow(data_temp) /D
    delta_mean <- data_temp$theta_mu
    if(info_person$irt_mod[1] %in% c(5,6,7,8)){
      delta_sd <-  data_temp$theta_sigma
      correlation <- rep(fit_info$info_relation[["nu"]], each = nperson )
    }
    if(info_person$irt_mod[1] %in% c(9,10)){
      delta_sd <-  data_temp$theta_sigma
      correlation <- rep(fit_info$info_relation[["nu"]], each = nperson ) * rep(fit_info$info_relation[["lambda"]], D)
    }
    if(info_person$irt_mod[1] %in% c(11,12)){
      delta_sd <-  sqrt(data_temp$theta_sigma^2 + data_temp$info_relation[["unique_sigma"]]^2)
      correlation <- data_temp$theta_sigma /delta_sd
    }
    df_long <- data.frame(
      age = rep(data_temp$age, 3),
      value = c(delta_mean,delta_sd, correlation),
      moment = factor(
        rep(c("mean", "standard deviation", "correlation"), each = D*nperson),
        levels = c("mean",  "standard deviation", "correlation")),
      trait = rep(data_temp$trait,3))
  }
  p <- ggplot2::ggplot(df_long) +
    ggplot2::geom_line(ggplot2::aes(age, value), linewidth = 1.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size=13,family = "serif"),
                   axis.line = ggplot2::element_line(color='black'),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank()) +
    ggplot2::ylab("latent trait")
  if(info_person$irt_mod[1] <= 4){
    p <- p +ggplot2::facet_wrap(~moment, scales = "free", nrow = 1)
  }
  if(info_person$irt_mod[1] > 4){
    p <- p + ggplot2::facet_wrap(moment~trait, scales = "free_y", 2)
  }

  print(p)
  invisible(list(plot = p, data = df_long))
  
}
