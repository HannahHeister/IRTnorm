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

#' A data set containing item responses and age information based on the GRM.
#' @docType data
#' @usage data(data_GRM)
#' @source Simulated data.
"data_GRM"

#' A data set containing item responses and age information based on the 1PLnorm model.
#' @docType data
#' @usage data(my_data)
#' @source Simulated data.
"my_data"

#' Presaved info for the vignette to run fast.
#'
#' @docType data
#' @usage data(info_small)
#' @format A \code{list} with information of the fitted model.
#' @source output from the estimate_norming_model function.
"info"

#' Presaved information of a newly normed individual for the vignette to run fast.
#'
#' @docType data
#' @usage data(int_info_newperson)
#' @format list containing the norm score estimates of a newly normed individual.
#' @source output from the norm_new_individual function.
"int_info_newperson"

#' Presaved different_poly_2PL for the vignette to run fast.
#'
#' @docType data
#' @usage data(agemodel_2PL)
#' @format A \code{list} with information for model comparison. 
#' @source output from the compare_age_models function.
"agemodel_2PL"

#' Presaved different_poly_2PL for the vignette to run fast.
#'
#' @docType data
#' @usage data(agemodel_2PL_excl)
#' @format A \code{list} with information for model comparison. 
#' @source output from the compare_age_models function based on data after 
#' excluding the two highly influentional item responses.
"agemodel_2PL_excl"


#' Presaved different_IRTmodels for the vignette to run fast.
#'
#' @docType data
#' @usage data(IRTmods)
#' @format A \code{list} with information for model comparison. 
#' @source output from the compare_IRTnorm_models function.
"IRTmods"

#' Presaved output of item_residual_cor for the 2PL model. 
#'
#' @docType data
#' @usage data(info_residualcor)
#' @source output of item_residual_cor for the 2PL model.
"info_residualcor"

#' Presaved testdata .
#'
#' @docType data
#' @usage data(test)
#' @format A \code{list} with information for model comparison. 
#' @source output from the compare_IRTnorm_models function.
"test"

#' Presaved testdata .
#'
#' @docType data
#' @usage data(HOex_cor_info)
#' @format A \code{list} with information for model comparison. 
#' @source output from the compare_IRTnorm_models function.
"cor_info"


#' A dataset containing item responses and age information based on the HO1PLnorm-C model.
#' @docType data
#' @usage data(HOdata)
#' @source Simulated data.
"HOdata"

#' Presaved HOex_res_info_all in order for the vignette to run fast.
#'
#' @docType data
#' @usage data(HOex_res_info_all)
#' @format A \code{list} with information returned from item_residual_cor.
#' @source based on data from the HO example 
"HOex_res_info_all"

#' Presaved HOex_res_info1 in order for the vignette to run fast.
#'
#' @docType data
#' @usage data(HOex_res_info1)
#' @format A \code{list} with information returned from item_residual_cor.
#' @source based on data from the HO example  for the first domain
"HOex_res_info1"


#' Presaved GRM model fit in order for the vignette to run fast.
#'
#' @docType data
#' @usage data(infoGRM)
#' @format A \code{list} with seleted elements of the information returned from estimate_norming_model.
#' @source based on data from the GRM example
"infoGRM"