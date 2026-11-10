#' Observed and model-implied raw score distributions
#'
#' Compares the observed distribution of raw scores with the distribution
#' implied by a fitted IRT norming model (posterior predictive check). This
#' shows how well the model reproduces the raw scores in the data, either for
#' the whole sample or separately for age groups. Use
#' \code{\link{visual_raw_score_ppc}} to plot the result.
#'
#' @param fit_info An object of class \code{"IRTnorm"}, as returned by
#'   \code{\link{estimate_norming_model}}.
#' @param age_range Optional numeric vector with at least two values giving the
#'   boundaries of the age groups, on the same scale as the age variable in the
#'   data (e.g., \code{c(5, 7, 9, 11)} creates the groups \[5,7\], \[7,9\], and
#'   \[9,11\]). If \code{NULL} (default), the whole sample is analyzed. Persons
#'   outside the range of \code{age_range} are excluded from the age-group
#'   results, and a warning reports how many.
#'
#' @details
#' The observed raw score is the sum of a person's item responses. For each
#' raw score, the observed frequency is compared with the frequencies in data
#' replicated from the fitted model. The 90\% credible interval (Q5 to Q95) of
#' the replicated frequencies gives the range of frequencies that is plausible
#' under the model. The model fits well if the observed frequencies lie within
#' or close to this interval across the score range. For background on
#' posterior predictive checks of IRT models, see Sinharay, Johnson, and
#' Stern (2006).
#'
#' For multidimensional models, the analysis is run for the total score over
#' all items and, in addition, for the raw score of each domain.
#'
#' Raw scores without any observations are kept in the output with frequency 0.
#'
#' @return An object of class \code{"IRTnorm_ppc"}, a list with one element
#'   for the total score and, for multidimensional models, one additional
#'   element per domain (in domain order).
#'
#'   If \code{age_range = NULL}, each element is a list with:
#'   \describe{
#'     \item{\code{df}}{Data frame with one row per possible raw score and the
#'       columns \code{score}, \code{obs} (observed frequency), \code{Q5},
#'       \code{Q50}, and \code{Q95} (5th, 50th, and 95th percentile of the
#'       model-implied frequencies).}
#'     \item{\code{replicated}}{Matrix of model-implied frequencies, with one
#'       row per posterior draw and one column per raw score.}
#'     \item{\code{nitem}}{Number of items the score is based on.}
#'   }
#'   If \code{age_range} is specified, each element is itself a list with one
#'   entry per age group. Each entry has the same components as above, plus
#'   \code{age}, a label describing the age group (e.g., \code{"age: [30,60]"}).
#'
#' @examples
#' \dontrun{
#' data(response_data)
#' info <- estimate_norming_model(
#'   raw_data = response_data,
#'   age_variable = "age",
#'   item_variables = 1:50,
#'   irt_model = "2PLnorm",
#'   age_model = "poly",
#'   poly_mean = 3,
#'   poly_sd = 2
#' )
#'
#' # Whole sample
#' result <- raw_score_ppc(fit_info = info)
#' visual_raw_score_ppc(result)
#'
#' # Separately by age group
#' result_age <- raw_score_ppc(fit_info = info, age_range = c(5, 7, 9, 11))
#' visual_raw_score_ppc(result_age)
#' }
#'
#' @seealso \code{\link{visual_raw_score_ppc}} for plotting the results,
#'   \code{\link{estimate_norming_model}}
#'
#' @references
#' Sinharay, S., Johnson, M. S., & Stern, H. S. (2006). Posterior predictive
#' assessment of item response theory models. \emph{Applied Psychological
#' Measurement, 30}(4), 298-321.
#'
#' @export
raw_score_ppc <- function(fit_info,
               age_range = NULL) {
  if(!inherits(fit_info, "IRTnorm")){
    stop("fit_info has to be generated from estimate_norming_model()")
  }
  data <- fit_info$newdata
  model_info <- fit_info$model_info
  raw_data <- fit_info$raw_data
  item_variables <- fit_info$item_variables
  age_variable <- fit_info$age_variable

  if (!is.null(age_range)) {
    if(!is.numeric(age_range) || length(age_range) <2){
      stop("age_range must be a numeric vector with at least two values.")
    }
    ages <- raw_data[[age_variable]]
    lower <- min(age_range)
    upper <- max(age_range)
    n_out <- sum(!is.na(ages) & (ages < lower | ages > upper))
    if (n_out > 0) {
      warning(n_out, " individual(s) fall outside 'age_range' [", lower, ",", upper,
              "] and are excluded from the age-group results.", call. = FALSE)
    }
  }
  
  
  int_data <- raw_data[, item_variables, drop = FALSE]
  
  if(all(data$itemD == 1)){
    info_all <-  vector("list", 1)
  }
  if(any(data$itemD > 1)){
    info_all <-  vector("list", data$D +1)
  }
if(data$irt_mod !=  4){
  nitem <- length(item_variables)
}else{
  nitem <- sum(data$K_item)
}

  # Compute replicated raw score distributions (shared computation)
  score_rep <- apply(model_info$y_rep, 1, function(x) tapply(x, data$jj, sum))
  draws <- ncol(score_rep)

  # Compute observed scores
  observed_scores <- rowSums(int_data, na.rm = T)

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
      score_rep <- apply(model_info$y_rep[,data$dd == d, drop = FALSE], 1, function(x) tapply(x, data$jj[data$dd == d], sum))
      draws <- ncol(score_rep)

      # Compute observed scores
      observed_scores <- rowSums(int_data[,data$itemD == d , drop = FALSE], na.rm = T)
      if (is.null(age_range)) {
        info_all[[d +1]] <- compute_nc_overall(observed_scores, score_rep, nitem, draws, dim = d)
      } else {
        info_all[[d + 1]] <- compute_nc_by_age(observed_scores, raw_data, age_variable,
                                               age_range, score_rep, nitem, draws, dim = d)
      }
    }
  }
  class(info_all) <- "IRTnorm_ppc"
  return(info_all)
}


# Helper: Create full observed distribution with zeros for missing scores
#create_full_distribution <- function(score_obs, nitem) {
#  if (length(score_obs) == 0) {
#    return(data.frame(score = 0:nitem, obs = 0))
#  }
#  observed <- as.data.frame(table(score_obs))
#  colnames(observed) <- c("score", "obs")
#  observed$score <- as.numeric(as.character(observed$score))
#
#  # Fill in missing scores with zero counts
#  if (nrow(observed) != nitem + 1) {
#    full_dist <- data.frame(score = 0:nitem, obs = 0)
#    full_dist$obs[match(observed$score, full_dist$score)] <- observed$obs
#    return(full_dist)
#  }
#
#  return(observed)
#}


create_full_distribution <- function(score_obs, nitem) {
  data.frame(
    score = 0:nitem,
    obs   = as.vector(table(factor(score_obs, levels = 0:nitem)))
  )
}

# Helper: Compute replicated score distributions
compute_replicated_distributions <- function(score_rep, age_filter = NULL,
                                             nitem, draws) {
  replicated <- matrix(nrow = draws, ncol = nitem + 1)

  for (i in seq_len(draws)) {
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

  return(list(
    df = cbind(observed_full, rep_quantiles),
    replicated = replicated,
    nitem = nitem
  ))
}

# Compute raw_score_ppc for overall sample
compute_nc_overall <- function(observed_scores, score_rep, nitem, draws, dim = 0) {
  result <- process_group(observed_scores, score_rep, nitem, draws)
  result$age <- "full age range"
  result$dim <- dim
  return(result)
}

# Compute raw_score_ppc by age groups
compute_nc_by_age <- function(observed_scores, raw_data, age_variable,
                              age_range, score_rep, nitem, draws, dim = 0) {
  age_groups <- length(age_range) - 1
  results <- vector("list", age_groups)

  ages <- raw_data[[age_variable]]

  for (i in seq_len(age_groups)) {
    # Define age filter
    age_filter <- ages >= age_range[i] & ages <= age_range[i + 1]

    # Process this age group
    result <- process_group(
      score_obs = observed_scores[age_filter],
      score_rep = score_rep,
      nitem = nitem,
      draws = draws,
      age_filter = age_filter
    )

    # Add age label
    result$age <- sprintf("age: [%s,%s]", age_range[i], age_range[i + 1])
    result$dim <- dim
    results[[i]] <- result
  }

  return(results)
}

#' Plot observed and model-implied raw score distributions
#'
#' Plots the observed raw score distribution together with the distribution
#' implied by the fitted IRT norming model (posterior predictive check), as
#' computed by \code{\link{raw_score_ppc}}. This shows how well the model
#' reproduces the raw scores in the data. The plot is printed and also
#' returned.
#'
#' @param raw_score_ppc_info An object of class \code{"IRTnorm_ppc"}, as
#'   returned by \code{\link{raw_score_ppc}}, either for the whole sample or
#'   for age groups.
#' @param dim_plot Character string indicating which latent trait to plot:
#'   \code{"general"} or \code{"domainK"}, where \code{K} is the domain number
#'   (e.g., \code{"domain1"}). Only required for multidimensional models in
#'   combination with age groups, because then only one trait can be shown at
#'   a time. It is ignored for all other models.
#'
#' @details
#' The plot shows, for each raw score:
#' \itemize{
#'   \item \strong{Black points}: the observed frequencies in the data.
#'   \item \strong{Orange points}: the model-implied frequencies, one point
#'     per posterior draw.
#'   \item \strong{Dashed line}: the median of the model-implied frequencies.
#'   \item \strong{Dotted lines}: the 5th and 95th percentiles of the
#'     model-implied frequencies, which together form a 90\% credible
#'     interval.
#' }
#' The model fits well if the observed frequencies lie within or close to the
#' credible interval and follow the pattern of the median line.
#'
#' The layout depends on the input:
#' \itemize{
#'   \item Whole sample, unidimensional model: a single panel.
#'   \item Age groups: one panel per age group.
#'   \item Whole sample, multidimensional model: one panel per trait (the
#'     general trait and each domain).
#'   \item Age groups, multidimensional model: one panel per age group for the
#'     trait selected with \code{dim_plot}.
#' }
#' Each panel has its own axis scales.
#'
#' @return A \code{ggplot} object, returned invisibly. It can be customized
#'   with additional \pkg{ggplot2} layers or saved with
#'   \code{\link[ggplot2]{ggsave}}.
#'
#' @examples
#' \dontrun{
#' data(response_data)
#' info <- estimate_norming_model(
#'   raw_data = response_data,
#'   age_variable = "age",
#'   item_variables = 1:50,
#'   irt_model = "2PLnorm",
#'   age_model = "poly",
#'   poly_mean = 3,
#'   poly_sd = 2
#' )
#'
#' # Whole sample
#' result <- raw_score_ppc(fit_info = info)
#' visual_raw_score_ppc(result)
#'
#' # Separately by age group
#' result_age <- raw_score_ppc(fit_info = info, age_range = c(-2,-1,0,1,2))
#' visual_raw_score_ppc(result_age)
#'
#' }
#'
#' @seealso \code{\link{raw_score_ppc}} for computing the raw score
#'   distributions
#'
#' @references
#' Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B. (2013).
#' \emph{Bayesian Data Analysis} (3rd ed.). Chapman and Hall/CRC.
#'
#' Sinharay, S., Johnson, M. S., & Stern, H. S. (2006). Posterior predictive
#' assessment of item response theory models. \emph{Applied Psychological
#' Measurement, 30}(4), 298-321.
#'
#' @export
visual_raw_score_ppc <- function(raw_score_ppc_info, dim_plot = NULL) {

    if(!inherits(raw_score_ppc_info, "IRTnorm_ppc")){
    stop("raw_score_ppc_info has to be generated from raw_score_ppc()")
  }

  # existing unwrap block
  if (length(raw_score_ppc_info) == 1 && is.list(raw_score_ppc_info[[1]]) &&
      !is.null(raw_score_ppc_info[[1]][[1]]$age)) {
    raw_score_ppc_info <- raw_score_ppc_info[[1]]
  }
 
process_single_result <- function(result, age_label = NULL, dim_label = NULL) {
  # Extract data
  df <- result$df
  replicated2 <- data.frame(result$replicated)
  colnames(replicated2) <- seq(0, result$nitem, 1)
  replicated2$id <- 1:nrow(replicated2)

  # Reshape to long format
  rep.long <- replicated2 |>
    tidyr::pivot_longer(cols = -dplyr::all_of("id"), names_to = "score", values_to = "frequency")
  rep.long$score <- as.numeric(rep.long$score)


  # Add age group label if provided
  if (!is.null(age_label)) {
    df$age_group <- age_label
    rep.long$age_group <- age_label
  }
  # Add domain label if provided
  if (!is.null(dim_label)) {
    df$domain <- dim_label
    rep.long$domain <- dim_label
  }
  return(list(df = df, rep.long = rep.long))
}

if (!(length(raw_score_ppc_info) > 1 && !is.null(raw_score_ppc_info[[1]][[1]]$age)) &&
    !is.null(dim_plot)) {
  warning("`dim_plot` is only used for multidimensional models with age groups ",
          "and is ignored here.")
  dim_plot <- NULL
}
# True if multiple dimensions are tested and if multiple ages are evaluated
if(length(raw_score_ppc_info) > 1 & !is.null(raw_score_ppc_info[[1]][[1]]$age)){
  if(is.null(dim_plot)){
    stop("'dim_plot' missing. Please provide which latent trait dimension  should be visualized. Due to
        the different age groups the latent trait dimensions have to be plotted individually.")
  }else{
    possible_dim <- c("general", paste0("domain", seq_len(length(raw_score_ppc_info) -1)))
    if(!dim_plot %in% possible_dim){
      dim_ok <- paste0(possible_dim, collapse = ", ")
      stop(paste0("dim_plot has to be one of the following values ",dim_ok ))
    }
    if(dim_plot == "general"){
      raw_score_ppc_info <- raw_score_ppc_info[[1]]
    }else{
        raw_score_ppc_info <- raw_score_ppc_info[[as.numeric(as.numeric(sub("domain", "", dim_plot))) +1]] 

    }
  }
}

# Detect if input is a list of lists or a single result
n_group <- length(raw_score_ppc_info)

if (n_group > 1) {
  # Process multiple  groups
  df_all <- NULL
  rep.long_all <- NULL

  for (g in seq_len(n_group)) {
    processed <- process_single_result(result = raw_score_ppc_info[[g]],
                                       age_label =raw_score_ppc_info[[g]]$age,
                                       dim_label =raw_score_ppc_info[[g]]$dim)
    df_all <- rbind(df_all, processed$df)
    rep.long_all <- rbind(rep.long_all, processed$rep.long)

  }
}else {
  # Process single result (no age grouping)
  processed <- process_single_result(raw_score_ppc_info[[1]], 
                                     age_label = raw_score_ppc_info[[1]]$age,
                                     dim_label = "none")
  df_all <- processed$df
  rep.long_all <- processed$rep.long

}
n_domains <- length(unique(rep.long_all$domain)) > 1
if(n_domains){
  rep.long_all$domain <- factor(rep.long_all$domain,
                                labels = c("general", paste0("domain", seq_len(max(rep.long_all$domain)))))
  df_all$domain <- factor(df_all$domain,
                          labels = c("general", paste0("domain", seq_len(max(df_all$domain)))))
}
# Create plot
p <- ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::geom_jitter(
    data = rep.long_all,
    ggplot2::aes(x = .data$score, y = .data$frequency),
    color = "#fdb863",
    position = ggplot2::position_jitter(width = 0.12),
    alpha = 0.15
  ) +
  ggplot2::geom_line(
    data = df_all,
    ggplot2::aes(x = .data$score, y = .data$Q50), #+1 and other quantiles
    linewidth = 1,
    linetype = "dashed"
  ) +
  ggplot2::geom_line(
    data = df_all,
    ggplot2::aes(x = .data$score , y = .data$Q5),
    linewidth = 1,
    linetype = "dotted"
  ) +
  ggplot2::geom_line(
    data = df_all,
    ggplot2::aes(x = .data$score , y = .data$Q95),
    linewidth = 1,
    linetype = "dotted"
  ) +
  ggplot2::geom_point(
    data = df_all,
    ggplot2::aes(x = .data$score , y = .data$obs),
    size = 4
  ) +
  ggplot2::xlab("raw score") +
  ggplot2::ylab("frequency of individuals") +
  ggplot2::theme(
    text = ggplot2::element_text(size = 13),
    axis.line = ggplot2::element_line(color = 'black'),
    plot.background = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    legend.position = "none"
  )

n_age <- length(unique(rep.long_all$age_group))
if(n_age > 1){
  p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$age_group), scale = "free", ncol = round(sqrt(n_age)))
}
if(n_domains){
  n_domain <- length(unique(rep.long_all$domain))
  p <-  p + ggplot2::facet_wrap(ggplot2::vars(.data$domain), scale = "free", ncol = round(sqrt(n_domain)))
}
if(!is.null(dim_plot)){
  p <- p + ggplot2::ggtitle(dim_plot)
}
print(p)
return(invisible(p))
}

#' Visualize the distribution of norm scores across age groups
#'
#' Plots the estimated norm scores of the individuals in the data as a
#' histogram for each age group, together with the density of a standard
#' normal distribution. This shows whether the norming achieved the intended
#' normalization of the latent trait within age groups: if it did, the
#' histogram should follow the curve in every age group.
#'
#' @param fit_info An object of class \code{"IRTnorm"}, as returned by
#'   \code{\link{estimate_norming_model}}.
#' @param age_range Numeric vector with at least two values giving the
#'   boundaries of the age groups, on the same scale as the age variable in the
#'   data (e.g., \code{c(5, 7, 9, 11)} creates the groups \[5,7\], \[7,9\], and
#'   \[9,11\]). The lower boundary of a group is excluded and the upper boundary
#'   is included.
#' @param binwidth Numeric. Width of the histogram bins, in norm score units.
#'   Default is 0.5.
#' @param dim_plot Only required for multidimensional models. Character string
#'   indicating which trait to plot: \code{"general"} or \code{"domainK"},
#'   where \code{K} is the domain number (e.g., \code{"domain1"}). For
#'   multidimensional models, the domain scores are shown on the standardized
#'   scale.
#'
#' @details
#' The norm score of each individual is the posterior mean of the
#' standardized latent trait. The normal curve is scaled to the number of
#' individuals in the age group and to \code{binwidth}, so that it is directly
#' comparable to the histogram frequencies. Each age group is shown in its own
#' panel.
#'
#' Individuals outside the range of \code{age_range} are removed from the
#' plot, with a warning.
#'
#' @return A \code{ggplot} object, which can be printed directly or further
#'   customized with additional \pkg{ggplot2} layers.
#'
#' @seealso \code{\link{estimate_norming_model}},
#'   \code{\link{raw_score_ppc}}
#'
#' @examples
#' \dontrun{
#' load(response_data)
#' info <- estimate_norming_model(
#'   raw_data = response_data,
#'   age_variable = "age",
#'   item_variables = 1:50,
#'   irt_model = "2PLnorm",
#'   age_model = "poly",
#'   poly_mean = 3,
#'   poly_sd = 2
#' )
#'
#' # Norm score distribution in age groups
#' normscore_dist(fit_info = info, age_range = c(-2,1,0,1,2))
#'
#' }
#'
#' @export
normscore_dist <- function(fit_info, age_range, binwidth = 0.5, dim_plot = NULL){
  if(!inherits(fit_info, "IRTnorm")){
    stop("fit_info has to be generated from estimate_norming_model()")
  }
    if(!is.numeric(age_range) || length(age_range) <2){
      stop("age_range must be a numeric vector with at least two values.")
    }
  
  age_range <- sort(age_range)
  model_info <- fit_info$model_info 
  if(model_info$info_sample$irt_mod[1] <= 4){
    if(!is.null(dim_plot)){
      warning("Your model is unidimensional, the argument dim_plot is ignored.")
      dim_plot <- NULL
    }
    data <-  model_info$info_sample
  }
  if(model_info$info_sample$irt_mod[1] > 4){
    if(is.null(dim_plot)){
      stop("Your data is from a multidimensional IRT model. Choose which dimension
           you want to visualize by specifying 'dim_plot'.")
    }
    D <- fit_info$newdata$D
    if (!dim_plot %in% c("general", paste0("domain", seq_len(D)))) {
      stop("`dim_plot` must be one of: ",
           paste(c("general", paste0("domain", seq_len(D))), collapse = ", "), ".")
    }
    
    if(dim_plot == "general"){
      data <-  model_info$info_sample[model_info$info_sample$trait == "theta", ]
    }else{
      data <-  model_info$info_sample[model_info$info_sample$trait == paste0("normeddelta",as.numeric(sub("domain", "", dim_plot))), ]
    }
    data$normscore_est <- data$est_ability
  }
  if (min(data$normscore_est) == max(data$normscore_est)) {
    stop("All norm scores are identical. Cannot plot a meaningful distribution.")
  }
  
  
  data$age_group <- NA
  n_agegroup <- (length(age_range) - 1)
  n_group <-NULL
  for(age in seq_len(n_agegroup)){
    data$age_group[data$age >= age_range[age] & data$age <= age_range[age +1]] <- paste0("age: [", age_range[age], ",", age_range[age +1], "]")
    n_group <- c(n_group,sum(data$age >= age_range[age] & data$age <= age_range[age +1]) )
  }
  if(any(is.na(data$age_group))){
    warning("Some individuals are not within the plotted age groups. These have been removed")
    data <- data[!is.na(data$age_group), ]
  }
  extrem <- max(abs(data$normscore_est))
  x_val <- seq(-extrem,extrem, length.out = 300)  # standard x values
  
  curve_list <- lapply(seq_along(n_group), function(i) {
    data.frame(
      x = x_val,
      y = stats::dnorm(x_val, 0, 1) * binwidth * n_group[i],
      age_group = paste0("age: [", age_range[i], ",", age_range[i + 1], "]")
    )
  })
  
  curve_data <- do.call(rbind, curve_list)
  
  ggplot2::ggplot()  +
    ggplot2::geom_histogram(ggplot2::aes(x = .data$normscore_est), binwidth = binwidth ,
                            data = data, color = "#fdb863", fill = "#fdb863") +
    ggplot2::geom_line(data = curve_data, ggplot2::aes(.data$x, .data$y), linewidth = 1) +
    ggplot2::facet_grid(~factor(age_group,levels = paste0("age: [", age_range[1:(n_agegroup)], ",", age_range[-1], "]"))) + 
    ggplot2::xlab("norm score") +
    ggplot2::ylab("frequency of individuals") + ggplot2::theme_bw()
  
  
}



# shared plotting theme, defined once
irt_plot_theme <- function(){
  ggplot2::theme(
    text             = ggplot2::element_text(size = 13),
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


#' Plot latent trait estimates and model-implied quantile curves
#'
#' Plots the estimated latent trait of each individual against age, together
#' with the quantile curves implied by the fitted IRT norming model. This shows
#' whether the model-implied distribution of the latent trait matches the
#' spread of the estimates across the age range: if it does, the expected
#' share of points should fall below each curve (e.g., about 25\% below the
#' 25\% curve).
#'
#' @param fit_info An object of class \code{"IRTnorm"}, as returned by
#'   \code{\link{estimate_norming_model}}.
#' @param probs Numeric vector of probabilities strictly between 0 and 1,
#'   giving the quantiles to plot. Default is
#'   \code{c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)}.
#'
#' @details
#' The points are the estimated latent trait values of the individuals in the
#' data. For multidimensional models, the points are the domain-specific
#' traits, and each domain is shown in its own panel.
#'
#' The curves show the quantiles of the age-dependent normal distribution
#' of the latent trait implied by the fitted model. They are drawn as smoothed
#' curves across age. If there are too few distinct ages for smoothing, the
#' unsmoothed quantiles are connected by lines instead.
#'
#' @return Invisibly, a data frame with the model-implied quantiles for each
#'   individual, with the columns \code{age}, \code{quan} (quantile value),
#'   \code{qname} (the probability, e.g., \code{"p0.05"}) and, for
#'   multidimensional models, \code{domain}. It can be used to create custom
#'   plots.
#'
#' @seealso \code{\link{estimate_norming_model}},
#'   \code{\link{normscore_dist}}
#'
#' @examples
#' \dontrun{
#' data(response_data)
#' info <- estimate_norming_model(
#'   raw_data = response_data,
#'   age_variable = "age",
#'   item_variables = 1:50,
#'   irt_model = "2PLnorm",
#'   age_model = "poly",
#'   poly_mean = 3,
#'   poly_sd = 2
#' )
#'
#' # Default quantiles
#' quantile_curves(fit_info = info)
#'
#' # Selected quantiles; keep the underlying values
#' q <- quantile_curves(fit_info = info, probs = c(0.025, 0.5, 0.975))
#' }
#'
#' @export
quantile_curves <- function(fit_info, 
                            probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)){
  if(!inherits(fit_info, "IRTnorm")){
    stop("fit_info has to be generated from estimate_norming_model()")
  }
  if(any(probs >= 1| probs <= 0)){
    stop("probs have to be in (0,1)")
  }
  if(!is.numeric(probs) || any(probs < 0 | probs > 1)){
    stop("probs have to be numerical values in (0,1)")
  }
  
  model_info  <- fit_info$model_info
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
      qname = factor(rep(paste0("p", probs), nperson), levels = paste0("p", probs))
    )
  } else {
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
        qname  = factor(rep(paste0("p", probs), nperson), levels = paste0("p", probs)),
        domain = d
      )
    }))
    person_data$theta_est <- person_data$est_ability
      }
  
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(.data$age, .data$theta_est), data = person_data) +
    ggplot2::ylab("latent trait") + ggplot2::theme_bw() + irt_plot_theme()
  
  # check if enough data to use lm with loess
  n_per_group <- tapply(new_data$age, new_data[c("qname", if (is_multidim(irt_mod)) "domain")],
                        function(a) length(unique(a)))
  span <- 0.75; degree <- 2
  min_n <- ceiling((degree + 1) / span) + 1 
  if (min(n_per_group) >= min_n){
    p <- p + ggplot2::geom_smooth(ggplot2::aes(.data$age, .data$quan, color = .data$qname), data = new_data, formula = 'y ~ x', method = 'loess')
  }else{
   p <- p + ggplot2::geom_line(ggplot2::aes(.data$age,.data$quan, color = .data$qname), data = new_data) 
  }

  if (is_multidim(irt_mod)) p <- p + ggplot2::facet_grid(~.data$domain)
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
#grm_probs <- function(theta, a,b){

#  cum_prob <- stats::plogis(a * (theta- b))
#if(length(theta) == 1){
#  p <- c(1, cum_prob) - c(cum_prob, 0)
#}else{
#  p <- cbind(1, cum_prob) - cbind(cum_prob, 0)
#}
#  return(p)
#}

grm_probs <- function(theta, a, b) {
  theta <- as.vector(theta)                         # also accepts a J x 1 matrix
  cum   <- stats::plogis(a * outer(theta, b, "-"))  # J x (K-1): P(Y >= k)
  cbind(1, cum) - cbind(cum, 0)                     # J x K:     P(Y = k)
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

  for(i in seq_len(I)){

    Ki <- K_item[i]
    responses <- info[, i]

    for(k in seq_len(Ki)){
      dummy[, col_index] <- as.integer(responses == k)
      col_names <- c(col_names, paste0("Item", i, "_Cat", k))
      col_index <- col_index + 1
    }
  }

  colnames(dummy) <- col_names

  return(dummy)
}



#' Plot item fit and person fit
#'
#' Compares the observed responses with the responses implied by the fitted
#' IRT norming model, either per item (item fit) or per person (person fit).
#' The closer the points lie to the diagonal reference line, the better the
#' model reproduces the data.
#'
#' @param fit_info An object of class \code{"IRTnorm"}, as returned by
#'   \code{\link{estimate_norming_model}}.
#' @param parameter Character string, either \code{"item"} or
#'   \code{"person"}.
#'   \itemize{
#'     \item \code{"item"}: compares the observed and model-implied
#'       proportion of correct responses for each item. See Details.
#'     \item \code{"person"}: compares the observed and model-implied raw
#'       score of each person.
#'   }
#'
#' @details
#' The model-implied values are computed from the item parameter estimates
#' and the estimated latent trait of each person (posterior means). Missing
#' responses are ignored in both the observed and the model-implied values.
#'
#' \strong{Item fit.} For the 1PL, 2PL, and 3PL models, each point is one
#' item, with the observed proportion of correct responses on the x-axis and
#' the model-implied proportion on the y-axis. For multidimensional models,
#' the items are shown in one panel per domain. Items without any valid
#' response are excluded, with a warning. For the graded response model, the
#' plot shows, for each item, the observed and model-implied proportion in
#' each response category as stacked bars.
#'
#' \strong{Person fit.} Each point is one person, with the observed raw score
#' on the x-axis and the model-implied expected raw score on the y-axis. For
#' multidimensional models, the raw score is computed separately for each
#' domain and shown in one panel per domain.
#'
#' @return Invisibly, a list with two components: \code{plot}, the
#'   \code{ggplot} object, and \code{data}, the data frame underlying the
#'   plot. The plot is also printed when the function is called.
#'
#' @seealso \code{\link{estimate_norming_model}},
#'   \code{\link{raw_score_ppc}}
#'
#' @examples
#' \dontrun{
#' data(response_data)
#' info <- estimate_norming_model(
#'   raw_data = response_data,
#'   age_variable = "age",
#'   item_variables = 1:50,
#'   irt_model = "2PLnorm",
#'   age_model = "poly",
#'   poly_mean = 3,
#'   poly_sd = 2
#' )
#'
#' # Item fit
#' item_person_fit(fit_info = info, parameter = "item")
#'
#' # Person fit, keeping the plot and the underlying data
#' res <- item_person_fit(fit_info = info, parameter = "person")
#' res$plot
#' head(res$data)
#' }
#'
#' @export
item_person_fit <- function(fit_info,  parameter){
  if(!inherits(fit_info, "IRTnorm")){
    stop("fit_info has to be generated from estimate_norming_model()")
  }
  parameter <- match.arg(parameter, c("item", "person"))
  model_info <- fit_info$model_info
  answer_pattern <- fit_info$raw_data[, fit_info$item_variables,drop = FALSE]
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
      ncol = sum(data$K_item) 
    )
    if(is.list(info_item$threshold)){
      expected[, 1:(sum(data$K_item[1]))] <- grm_probs(theta = info_person, a = info_item$alpha[1], b = info_item$threshold[[1]])
      for(i in seq_len(nitem)[-1]){
        expected[, (1+ sum(data$K_item[1:(i-1)])):(sum(data$K_item[1:i]))] <- grm_probs(theta = info_person, a = info_item$alpha[i], b = info_item$threshold[[i]])
      }
    }else{
      int_threshold <- seq(1,nitem*data$K_max, nitem)[1:(data$K_item[1] -1)]
      expected[, 1:(sum(data$K_item[1]))] <- grm_probs(theta = info_person, a = info_item$alpha[1], b = info_item$threshold[int_threshold])
      for(i in seq_len(nitem)[-1]){
        int_threshold <- seq(i,nitem*data$K_max, nitem)[1:(data$K_item[i] -1)]
        expected[, (1+ sum(data$K_item[1:(i-1)])):(sum(data$K_item[1:i]))] <- grm_probs(theta = info_person, a = info_item$alpha[i], b = info_item$threshold[int_threshold])
      }
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
      df$item_label <- factor(paste0("item ", df$item),
                              levels = paste0("item ", seq_len(nitem)))
    }
    if(data$irt_mod > 4){
      df <- data.frame(
        observed = colSums(answer_pattern, na.rm = TRUE) /colSums(!is.na(expected)),
        irt_mod    = colSums(expected, na.rm = TRUE) / colSums(!is.na(expected)),
        domain = paste0("domain ",data$itemD))
    }
    value_cols <- if (data$irt_mod == 4) "value" else c("observed", "irt_mod")
    bad <- Reduce(`|`, lapply(df[value_cols], is.nan))
    if(any(bad)){
      warning("One or more items have no valid responses. These items are excluded from the plot.")
      df <- df[!bad, ]
    }
    
    xname <- "observed"
    yname <- "model-implied"
    titlename <- "item fit: proportion correct per item"
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
          observed = rowSums(answer_pattern[,data$itemD == d ], na.rm = TRUE),
          irt_mod    = rowSums(expected[,data$itemD == d ], na.rm = TRUE),
          domain = paste0("domain ", d)))
      }
    }
    xname <- "observed"
    yname <- "model-implied"
    titlename <- "person fit: raw score"
  }
  if(data$irt_mod == 4 & parameter == "item"){
    p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(.data$type),
                                          y = .data$value, fill = factor(.data$level))) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::facet_wrap(ggplot2::vars(.data$item_label)) +
      ggplot2::labs(x = NULL, y = "proportion", fill = "Level") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
  }else{
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$observed, y = .data$irt_mod)) +
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
  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
  
  if(data$irt_mod > 4){
    p <- p + ggplot2::facet_grid(cols = ggplot2::vars(.data$domain))
  }
  print(p)
  invisible(list(plot = p, data = df))
}




#' Plot the estimated age effect of an IRT norming model
#'
#' Plots how the fitted IRT norming model relates age to the distribution of
#' the latent trait. For unidimensional models, the plot shows the
#' age-dependent mean and standard deviation of the latent trait. For
#' multidimensional models, it shows the age-dependent mean and standard
#' deviation of each domain and the correlation between each domain and the
#' general trait.
#'
#' @param fit_info An object of class \code{"IRTnorm"}, as returned by
#'   \code{\link{estimate_norming_model}}.
#'
#' @details
#' The curves show the model-implied values, based on the posterior mean
#' estimates, evaluated at the ages of the individuals in the data. They are
#' therefore only drawn over the age range that is observed in the data.
#'
#' For unidimensional models, the mean and the standard deviation are shown
#' in two panels.
#'
#'  For multidimensional models, the plot is a grid with one row per quantity
#' (mean, standard deviation, correlation) and one column per domain. The
#' y-axis scale is shared across the domains within a row, so the domains can
#' be compared directly. The correlation is the correlation between the
#' domain-specific trait and the general trait.
#'
#' @return Invisibly, a list with two components: \code{plot}, the
#'   \code{ggplot} object, and \code{data}, the data frame underlying the
#'   plot (columns \code{age}, \code{value}, \code{moment} and, for
#'   multidimensional models, \code{trait}).
#'   The plot is printed when the function is called.
#'
#' @seealso \code{\link{estimate_norming_model}},
#'   \code{\link{quantile_curves}}
#'
#' @examples
#' \dontrun{
#' info <- estimate_norming_model(
#'   raw_data = response_data,
#'   age_variable = "age",
#'   item_variables = 1:50,
#'   irt_model = "2PLnorm",
#'   age_model = "poly",
#'   poly_mean = 3,
#'   poly_sd = 2
#' )
#'
#' # Plot the estimated age effect
#' effect_age(fit_info = info)
#'
#' # Keep the plot and the underlying data
#' res <- effect_age(fit_info = info)
#' res$plot
#' head(res$data)
#' }
#'
#' @export
effect_age <- function(fit_info){
  if(!inherits(fit_info, "IRTnorm")){
    stop("fit_info has to be generated from estimate_norming_model()")
  }
  info_person <- fit_info$model_info$info_sample
  info_relation<- fit_info$model_info$info_relation
  
  if(info_person$irt_mod[1] <= 4){
    df_long <- data.frame(
      age = rep(info_person$age, 2),
      value = c(info_person$theta_mu, info_person$theta_sigma),
      moment = factor(
        rep(c("mean", "standard deviation"), each = nrow(info_person)),
        levels = c("mean", "standard deviation")
      )
    )
  }else{
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
      correlation <- rep(info_relation[["nu"]], each = nperson )
    }
    if(info_person$irt_mod[1] %in% c(9,10)){
      delta_sd <-  data_temp$theta_sigma
      correlation <- rep(info_relation[["nu"]], each = nperson ) * rep(info_relation[["lambda"]], D)
    }
    if(info_person$irt_mod[1] %in% c(11,12)){
      delta_sd <-  sqrt(matrix(data_temp$theta_sigma^2, ncol =D) + matrix(info_relation[["unique_sigma"]]^2, ncol =D))
      correlation <- matrix(data_temp$theta_sigma,ncol =D) /delta_sd
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
    ggplot2::geom_line(ggplot2::aes(.data$age, .data$value), linewidth = 1.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size=13),
                   axis.line = ggplot2::element_line(color='black'),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank()) +
    ggplot2::ylab("latent trait")
  if(info_person$irt_mod[1] <= 4){
    p <- p +ggplot2::facet_wrap(ggplot2::vars(.data$moment), scales = "free", nrow = 1)
  }
  if(info_person$irt_mod[1] > 4){
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data$moment),
                                 cols = ggplot2::vars(.data$trait),
                                 scales = "free_y")
  }

  print(p)
  invisible(list(plot = p, data = df_long))
  
}



#' Calculates the average correlation across and within domains. 
#'
#' @param Q3  A \code{matrix} of the item residual correlation
#' @param itemD A \code{numeric} vector of length equal to the number of items,
#'   giving the domain each item measures.
#'
#' @return A data.frame with the average item residual correlation within 
#' domains and across domains. 
#'
#' @noRd
compute_domain_averages <- function(Q3, itemD) {
  
  #J <- nrow(Q3)
  #within_vals <- c()
  #between_vals <- c()
  
  #for (i in seq_len(J - 1)) {
  #  for (j in (i + 1):J) {
   #   val <- Q3[i, j]
   #   if (is.na(val)) next
   #   if (itemD[i] == itemD[j]) {
   ##     within_vals <- c(within_vals, val)
   #   } else {
   #     between_vals <- c(between_vals, val)
    #  }
   # }
 # }
  
  ut <- upper.tri(Q3) & !is.na(Q3)
  same <- outer(itemD, itemD, "==")
  within_vals  <- Q3[ut & same]
  between_vals <- Q3[ut & !same]
  data.frame(mean_correlation = c(mean(within_vals, na.rm = TRUE),
                                  mean(between_vals, na.rm = TRUE)), 
             domains = c("within", "between"), 
             n_pairs = c(length(within_vals), length(between_vals)))
}

#' Visualizes Q3 for hierarchical models, highlighting which items belong to 
#' which domain by including horizontal and vertical lines. 
#'
#' @param Q3  A \code{matrix} of the item residual correlation
#' @param itemD A \code{numeric} vector of length equal to the number of items,
#'   giving the domain each item measures.
#'
#' @return A a heat map with lines indicating the items measuring the same domain 
#'
#' @noRd
visual_multi <- function(Q3 = Q3, itemD = itemD){
  breaks <- utils::head(cumsum(rle(itemD)$lengths), -1) + 1
  n <- ncol(Q3)
  
  # Horizontal lines: from the diagonal to the right edge
  h_seg <- data.frame(
    x    = breaks - 0.5,
    xend = n - 0.5,
    y    = breaks + 0.5,
    yend = breaks + 0.5
  )
  
  # Vertical lines: from the bottom edge up to the diagonal
  v_seg <- data.frame(
    x    = breaks - 0.5,
    xend = breaks - 0.5,
    y    = 0.5,
    yend = breaks + 0.5
  )
  
  p <- ggcorrplot::ggcorrplot(
    Q3,
    hc.order = FALSE,
    type = "lower",
    outline.color = "white",
    lab = FALSE,
    colors = c("#0D0887FF", "white", "#D98F3A")
  ) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::geom_segment(
      data = h_seg,
      ggplot2::aes(x = .data$x, xend = .data$xend, y =.data$ y, yend = .data$yend),
      inherit.aes = FALSE, color = "black", linewidth = 0.8
    ) +
    ggplot2::geom_segment(
      data = v_seg,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      inherit.aes = FALSE, color = "black", linewidth = 0.8
    )
  
  p
}



#' Compute item residual correlations (Q3)
#'
#' Checks the local independence assumption of the fitted IRT norming model.
#' For each person and item, the residual is the observed response minus the
#' model-implied expected response. The residual correlations between all item
#' pairs (Yen's Q3) should be close to zero if the model captures the
#' dependencies in the data. The function prints a summary of the correlations
#' between items and, optionally, plots them as a heat map.
#'
#' @param fit_info An object of class \code{"IRTnorm"}, as returned by
#'   \code{\link{estimate_norming_model}}.
#' @param visual Logical. If \code{TRUE} (default), the residual correlation
#'   matrix is plotted as a heat map.
#' @param domain Optional integer indicating which domain to analyze. Only
#'   usable for multidimensional models; the analysis is then restricted to the
#'   items of this domain. If \code{NULL} (default), all items are analyzed.
#'   For unidimensional models, the argument is ignored with a warning.
#'
#' @details
#' The expected responses are computed from the posterior mean item parameters
#' and the estimated latent trait of each person (posterior means). Missing
#' responses are ignored when computing the correlations (pairwise complete
#' observations).
#'
#' A summary (minimum, quartiles, mean, maximum) of the correlations between
#' all distinct item pairs is printed. Correlations far from zero, in
#' particular strongly positive ones, indicate pairs of items that are related
#' beyond what the model explains.
#'
#' For multidimensional models without a specified \code{domain}, the
#' function additionally prints the average residual correlation within and
#' between domains, and the heat map is organized by domain.
#'
#' The heat map and the printed output are produced when the function is
#' called. The numerical results are returned invisibly.
#'
#' @return Invisibly, a list with the components:
#' \describe{
#'   \item{\code{Q3}}{Matrix of residual correlations between all items
#'     analyzed.}
#'   \item{\code{summary_info}}{Summary of the correlations between distinct
#'     items (the off-diagonal elements of \code{Q3}).}
#'   \item{\code{itemD}}{Domain of each item. Only for multidimensional
#'     models without a specified \code{domain}.}
#'   \item{\code{compare_domains}}{Average residual correlations within and
#'     between domains. Only for multidimensional models without a specified
#'     \code{domain}.}
#' }
#'
#' @seealso \code{\link{estimate_norming_model}},
#'   \code{\link{item_person_fit}}
#'
#' @references
#' Yen, W. M. (1984). Effects of local item dependence on the fit and
#' equating performance of the three-parameter logistic model.
#' \emph{Applied Psychological Measurement, 8}(2), 125-145.
#'
#' @examples
#' \dontrun{
#' data(response_data)
#' info <- estimate_norming_model(
#'   raw_data = response_data,
#'   age_variable = "age",
#'   item_variables = 1:50,
#'   irt_model = "2PLnorm",
#'   age_model = "poly",
#'   poly_mean = 3,
#'   poly_sd = 2
#' )
#'
#' # Summary and heat map
#' item_residual_cor(info)
#'
#' # Keep the residual correlation matrix, without the plot
#' res <- item_residual_cor(info, visual = FALSE)
#' res$Q3[1:5, 1:5]
#'
#' # Multidimensional model: restrict to one domain
#' # item_residual_cor(info_multi, domain = 1)
#' }
#'
#' @export
item_residual_cor <- function(fit_info, visual = TRUE, domain = NULL) {
  if(!inherits(fit_info, "IRTnorm")){
    stop("fit_info has to be generated from estimate_norming_model()")
  }
  data <- fit_info$newdata
  J <- data$J
  I <- data$I
  raw_data <- fit_info$raw_data
  item_variables <- fit_info$item_variables
  raw_data  <- as.matrix(raw_data[,item_variables])
  irt_mod <- data$irt_mod
  
  if(irt_mod != 4){
    a_vec <- colMeans(fit_info$norm_info$alpha)
    b_vec <- colMeans(fit_info$norm_info$beta)
    c_vec <- colMeans(fit_info$norm_info$gamma)
  }
  if(irt_mod < 5){
    theta <- fit_info$model_info$info_sample$theta_est
  }
  if(!is.null(domain) && irt_mod <= 4){
    warning("The domain argument is only for hierarchical norming models. 
            You indicated a unidimensional model, therefore domain is ignored.")
  }
  if(!is.null(domain) && irt_mod > 4){
    itemD <- data$itemD
    if(!(domain %in% itemD) ||length(domain) != 1 ){
      stop("domain has to be a single numerical value indicating which of the 
           domains of the test are supposed to be inspected.")
    }
    irt_mod <- 2
    a_vec <- a_vec[itemD == domain]
    b_vec <- b_vec[itemD == domain]
    c_vec <- c_vec[itemD == domain]
    samp_info <-fit_info$model_info$info_sample
    theta <- samp_info$est_ability[samp_info$trait == paste0("delta", domain)]
    I <- sum(itemD == domain)
    raw_data <- raw_data[,itemD == domain]
  }
  
  E <- matrix(NA_real_, nrow = J, ncol = I)
  if(irt_mod < 4){
    for (i in 1:I) {
      E[, i]  <- c_vec[i] + (1-c_vec[i])*stats::plogis(a_vec[i] * (theta - b_vec[i] ))
    }
  }
  if(irt_mod == 4){
    a_vec <- fit_info$model_info$info_item$alpha
    b_vec <- fit_info$model_info$info_item$threshold
    #p <- NULL
    if(is.list(b_vec)){
      for (i in seq_len(I)) {
        k_i <- seq_len(data$K_item[i] - 1)
        bi <- b_vec[[i]]
        cum <- stats::plogis(a_vec[i] * outer(theta, bi, "-"))   # N x (K-1): P(Y >= k)
        p   <- cbind(1, cum) - cbind(cum, 0)                 # N x K:     P(Y = k)
        E[, i] <- drop(p %*% (seq_len(ncol(p)) - 1))
      }
    }else{
      for (i in seq_len(I)) {
        #k_i <- seq_len(data$K_item[i] - 1)
        #bi <- b_vec[paste0("b[", i, ",", k_i, "]")]
        bi <- (b_vec[seq(i,length(b_vec),I)])[data$K_item[i] -1]
        cum <- stats::plogis(a_vec[i] * outer(theta, bi, "-"))   # N x (K-1): P(Y >= k)
        p   <- cbind(1, cum) - cbind(cum, 0)                 # N x K:     P(Y = k)
        E[, i] <- drop(p %*% (seq_len(ncol(p)) - 1))
      }
    }

  }
  if(irt_mod > 4){
    D <- fit_info$newdata$D
    samp_info <-fit_info$model_info$info_sample
    delta <- matrix(samp_info$est_ability[substr(samp_info$trait,1,5) == "delta"], ncol = D)
    itemD <- fit_info$newdata$itemD
    for(i in 1:I){
      E[, i] <- c_vec[i] + (1-c_vec[i])*stats::plogis(a_vec[i] * (delta[,itemD[i]] - b_vec[i]))
    }
  }
  
  DIFF <-  raw_data  - E   
  
  Q3 <- stats::cor(DIFF, use = "pairwise.complete.obs")
  colnames(Q3) <- rownames(Q3) <- colnames( raw_data )
  off_diag <- Q3 [row(Q3 ) != col(Q3 )]
  
  summary_info <- summary(off_diag)
  print(summary_info)
  
  if(irt_mod > 4){
    cat("\n\n")
    compare_domains <- compute_domain_averages(Q3 = Q3, itemD = itemD)
    print(compare_domains)
  }
  if(visual){
    if(irt_mod > 4){
      p1 <- visual_multi(Q3 = Q3, itemD = itemD)
    }else{
      p1 <- ggcorrplot::ggcorrplot(Q3, hc.order = FALSE, type = "lower",
                                   lab = F, colors = c( "#0D0887FF", "white", "#D98F3A")) +
        ggplot2::theme(text = ggplot2::element_text(size=12),
                       axis.text.x = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank()) 
    }
    print(p1)
  }
  if(irt_mod > 4){
    invisible(list(Q3 = Q3, itemD = itemD,
                   compare_domains =compare_domains, summary_info = summary_info))
    
  }else{
    invisible(list(Q3 = Q3, summary_info = summary_info))
  }
}


#' Explore the model-implied correlations between domains
#'
#' Explores the relationship between the domains of a multidimensional
#' IRT norming model by returning the
#' model-implied correlations between the domain-specific traits.
#'
#' How the correlations are reported depends on the model:
#' \itemize{
#'   \item \strong{Constant correlations} (\code{"Testlet1PLnorm"},
#'     \code{"Testlet2PLnorm"}, \code{"HO1PLnorm-C"}, \code{"HO2PLnorm-C"}):
#'     the model-implied correlation of each pair of domains is returned
#'     together with the correlation of the observed domain sum scores, so
#'     that the two can be compared.
#'   \item \strong{Age-dependent correlations} (\code{"HO1PLnorm-P"},
#'     \code{"HO2PLnorm-P"}, \code{"HO1PLnorm-F"}, \code{"HO2PLnorm-F"}):
#'     the model-implied correlation of each pair of domains is returned for
#'     every individual, and can be plotted against age.
#' }
#'
#' @param fit_info An object of class \code{"IRTnorm"}, as returned by
#'   \code{\link{estimate_norming_model}}, fitted with a multidimensional
#'   model.
#' @param visual Logical. Only used for models with age-dependent
#'   correlations. If \code{TRUE} (default), the model-implied correlations
#'   are plotted against age, with one color per domain pair.
#'
#' @details
#' The model-implied correlations are based on the posterior mean estimates.
#' Domain pairs are labeled by their domain numbers (e.g., \code{"1&2"} for
#' domains 1 and 2).
#'
#' @return A data frame.
#'
#'   For models with constant correlations, it has one row per domain pair and
#'   the columns \code{sumscore_cor} (correlation of the observed domain sum
#'   scores), \code{model_cor} (model-implied correlation), and
#'   \code{domains} (the domain pair).
#'
#'   For models with age-dependent correlations, it has one row per
#'   individual, one column per domain pair with the model-implied correlation
#'   (named, e.g., \code{"cor(d1,d2)"}), and the column \code{age}.
#'
#' @seealso \code{\link{estimate_norming_model}},
#'   \code{\link{effect_age}}
#'
#' @examples
#' \dontrun{
#' data(response_data)
#' info <- estimate_norming_model(
#'   raw_data = response_data,
#'   age_variable = "age",
#'   item_variables = 1:50,
#'   irt_model = "HO2PLnorm-C",
#'   itemD = rep(1:3, each = 17)[-1],
#'   age_model = "poly",
#'   poly_mean = 3,
#'   poly_sd = 2
#' )
#'
#' # Model-implied vs. sum score correlations between the domains
#' cor_structure(info)
#' }
#'
#' @export
cor_structure <- function(fit_info, visual = TRUE){
  if(!inherits(fit_info, "IRTnorm")){
    stop("fit_info has to be generated from estimate_norming_model()")
  }
  
  normdata <- fit_info$newdata
  response_data <- fit_info$raw_data
  if(normdata$irt_mod <= 4){
    stop("this function is only for hierarchical IRTnorm models, you fit object 
         is from a unidimensional model")
  }
  itemD <- normdata$itemD
  D <- normdata$D
  i <- rep.int(seq_len(D - 1), (D - 1):1)
  j <- sequence((D - 1):1, from = 2:D)
  cor_mod_temp <- fit_info$model_info$info_relation[["nu"]]
  if(normdata$irt_mod < 9){
    temp <- as.matrix(response_data[,fit_info$item_variables, drop = FALSE])
    sumscore <- as.data.frame(t(rowsum(t(temp), group = itemD, na.rm= TRUE)))
    cor_sum <- stats::cor(sumscore)
    cor_mod <- cor_mod_temp[i] * cor_mod_temp[j]
    info <- data.frame(sumscore_cor = cor_sum[cbind(i, j)], 
                       model_cor = cor_mod, 
                       domains = paste0(i, "&", j))
    return(info)
  }
  
  if(normdata$irt_mod %in% c(9, 10)){
    cor_mod_temp <- outer(fit_info$model_info$info_relation[["lambda"]], cor_mod_temp)
  }
  if(normdata$irt_mod %in% c(11,12)){
    unique_sigma <- fit_info$model_info$info_relation[["unique_sigma"]]
    delta_sigma <- matrix(ncol = D, nrow = normdata$J)
    for(d in seq_len(D)){
      delta_sigma[,d] <- fit_info$model_info$info_sample$theta_sigma[fit_info$model_info$info_sample$trait == paste0("delta",d)]
    }
    cor_mod_temp <- delta_sigma / sqrt(delta_sigma^2 + unique_sigma^2)
  }
  cor_mod_temp2 <- matrix(NA, nrow = nrow(cor_mod_temp), ncol = length(i))
  for(ti in seq_along(i)){
    cor_mod_temp2[,ti] <- cor_mod_temp[,i[ti]] *cor_mod_temp[,j[ti]]
  }
  
  if(normdata$irt_mod > 8){
    cor_mod <- data.frame(cor_mod_temp2, age = fit_info$raw_data[fit_info$age_variable])
    colnames(cor_mod) <-  c(paste0("cor(d",i, ",d", j, ")"), "age")
    if(visual){
      pair_labels <- paste0("cor(d", i, ",d", j, ")")
      
      cor_long <- data.frame(
        combi = factor(rep(pair_labels, each = nrow(cor_mod)),
                       levels = pair_labels),
        cor_value = as.vector(as.matrix(cor_mod[, pair_labels])),
        age = rep(cor_mod$age, length(i))
      )
      p <- ggplot2::ggplot(cor_long) + 
        ggplot2::geom_point(ggplot2::aes(.data$age, .data$cor_value, color = .data$combi)) + 
        ggplot2::theme_bw()+ 
        ggplot2::ylab("correlation")
      print(p)
    }
    return(cor_mod)
  }
  
}

