#' Response data example dataset
#'
#' A dataset containing item responses and age information.
#' @docType data
#' @usage data(response_data)
#' @format A \code{data.frame} with \eqn{n} = 500 rows and \eqn{p}  = 56 columns:
#' \describe{
#'   \item{item1}{Numeric. Response to item 1 (0 = incorrect, 1 = correct).}
#'   \item{item2}{Numeric. Response to item 2 (0 = incorrect, 1 = correct).}
#'   ...
#'   \item{item50}{Numeric. Response to item 56 (0 = incorrect, 1 = correct).}
#'   \item{age}{Numeric. Age of the respondent in years.}
#'   
#' }
#' @source Simulated data.
"response_data"

#' Presaved conv_plot in order for the vignette to run fast.
#'
#' @docType data
#' @usage data(conv_plot)
#' @format A \code{ggplot2} object.
#' @source Simulated data.
"conv_plot"


#' Presaved info for the vignette to run fast.
#'
#' @docType data
#' @usage data(info)
#' @format A \code{list} with information of the fitted model.
#' @source Simulated data.
"info"

#' Presaved information of a newly normed individual for the vignette to run fast.
#'
#' @docType data
#' @usage data(int_info_newperson)
#' @format Output form loo compare function.
#' @source Simulated data.
"int_info_newperson"

#' Presaved NC_info for the vignette to run fast.
#'
#' @docType data
#' @usage data(NC_info)
#' @format A \code{list} with information for norming new individuals.
#' @source Simulated data.
"NC_info"


#' Presaved NC_info_age for the vignette to run fast.
#'
#' @docType data
#' @usage data(NC_info_age)
#' @format A \code{list} with information for norming new individuals.
#' @source Simulated data.
"NC_info_age"

#' Presaved different_poly_2PL for the vignette to run fast.
#'
#' @docType data
#' @usage data(different_poly_2PL)
#' @format A \code{list} with information for norming new individuals.
#' @source Simulated data.
"different_poly_2PL"

#' Presaved different_IRTmodels for the vignette to run fast.
#'
#' @docType data
#' @usage data(different_IRTmodels)
#' @format A \code{list} with information for norming new individuals.
#' @source Simulated data.
"different_IRTmodels"
