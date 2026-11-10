build_spline <- function(age_df, k, fixed_knots) {
  if (!is.null(fixed_knots)) {
    if (identical(fixed_knots, "equidistant")) {
      fixed_knots <- list(age = seq(min(age_df$age), max(age_df$age), length.out = k))
    } else {
      fixed_knots <- list(age = fixed_knots)
    }
    k <- length(fixed_knots$age)
    age <- NULL
    sm <- mgcv::smoothCon(
      mgcv::s(age, bs = "cr", k = k),
      knots       = fixed_knots,
      data        = age_df,
      absorb.cons = FALSE
    )[[1]]
  } else {
    sm <- mgcv::smoothCon(
      mgcv::s(age, bs = "cr", k = k),
      data        = age_df,
      absorb.cons = FALSE
    )[[1]]
  }
  
  re <- mgcv::smooth2random(sm, names(age_df), type = 2)
  X  <- re$Xf        # (J x 2) -- intercept + linear trend, UNPENALIZED
  Z  <- re$rand$Xr   # (J x k-2) -- wiggly deviations, PENALIZED
  
  Xraw <- mgcv::PredictMat(sm, age_df)
  R_Z  <- MASS::ginv(Xraw) %*% Z   # rotation for penalized part
  R_X  <- MASS::ginv(Xraw) %*% X   # rotation for unpenalized part
  
  list(
    sm = sm, X = X, Z = Z,
    S_fixed = ncol(X), K = ncol(Z),
    R_Z = R_Z, R_X = R_X
  )
}



#' Prepare data for fitting an IRT-based continuous norming model
#'
#' @description
#' Converts person-by-item response data and an age variable into the named
#' list required to fit an age-dependent Bayesian IRT norming model. Missing
#' item responses are dropped, age is standardized, and the age design
#' matrices (polynomial or spline) are constructed.
#'
#' This function is mainly called internally by the model-fitting function.
#' It is exported so that you can inspect the data list that is passed to the
#' sampler, which can be useful when interpreting the fitted stan object, 
#' such as in the example of the vingette: How_to_use_IRTnorm.
#'
#' @param raw_data A \code{data.frame} or \code{matrix} with one row per
#'   person. It must contain the age variable and the item response variables.
#'   Item responses may contain \code{NA}.
#' @param age_variable A column name (\code{character}) or column index
#'   (\code{numeric}) in \code{raw_data} giving the person's age. Age is
#'   treated as continuous. Persons with missing age are removed with a warning.
#' @param item_variables A \code{character} or \code{numeric} vector with the
#'   columns of \code{raw_data} that contain the item responses.
#' @param irt_model A \code{character} scalar giving the norming model. The
#'   first option is the default. All models have an age-dependent mean and
#'   standard deviation of the latent trait.
#'   \describe{
#'     \item{\code{"1PLnorm"}}{One-parameter logistic (Rasch) model.}
#'     \item{\code{"2PLnorm"}}{Two-parameter logistic model.}
#'     \item{\code{"3PLnorm"}}{Three-parameter logistic model.}
#'     \item{\code{"GRMnorm"}}{Graded response model for polytomous items.
#'       Requires \code{K_item}.}
#'     \item{\code{"Testlet1PLnorm"}, \code{"Testlet2PLnorm"}}{Higher-order
#'       1PL / 2PL model with age-dependent mean and standard deviation that
#'       are constant across domains. Requires \code{itemD}.}
#'     \item{\code{"HO1PLnorm-C"}, \code{"HO2PLnorm-C"}}{Higher-order 1PL / 2PL
#'       model with a \emph{constant} correlation matrix across age.
#'       Requires \code{itemD}.}
#'     \item{\code{"HO1PLnorm-P"}, \code{"HO2PLnorm-P"}}{Higher-order 1PL / 2PL
#'       model with a correlation matrix that is \emph{proportional} across
#'       age. Requires \code{itemD}.}
#'     \item{\code{"HO1PLnorm-F"}, \code{"HO2PLnorm-F"}}{Higher-order 1PL / 2PL
#'       model with a \emph{free} correlation matrix across age.
#'       Requires \code{itemD}.}
#'   }
#' @param age_model A \code{character} scalar specifying how age effects are
#'   modeled. The first option is the default.
#'   \describe{
#'     \item{\code{"splines"}}{Penalized B-spline basis functions of
#'       standardized age (see \code{k_mean}, \code{k_sd}, \code{k_sd2} and the
#'       \code{fixed_knots_*} arguments).}
#'     \item{\code{"poly"}}{Raw polynomials of standardized age (see
#'       \code{poly_mean}, \code{poly_sd}, \code{poly_sd_unique}).}
#'   }
#' @param itemD Only required for the Testlet and higher-order (HO) models.
#'   An integer vector of length equal to the number of items indicating the
#'   domain each item loads on. Domains must be numbered consecutively
#'   starting at 1 (e.g., \code{1, 2, 3, 4} for four domains).
#' @param K_item Only required for \code{irt_model = "GRMnorm"}. An integer
#'   vector of length equal to the number of items giving the number of
#'   response categories per item (minimum 2 for binary items). Responses must
#'   be coded starting at 1; if the smallest response is below 1, all
#'   responses are shifted and a warning is issued.
#' @param poly_mean,poly_sd Only required if \code{age_model = "poly"}.
#'   Single non-negative integers giving the polynomial degree for the
#'   age-dependent mean and standard deviation of the latent trait. A degree
#'   of 0 fits an intercept only (no age effect).
#' @param poly_sd_unique Only used if \code{age_model = "poly"} and
#'   \code{irt_model} is an HO-P or HO-F model. Polynomial degree for the
#'   age-dependent unique (domain-specific) standard deviation component.
#'   If \code{NULL}, the design matrix for \code{poly_sd} is reused.
#' @param k_mean,k_sd Only used if \code{age_model = "splines"}. Number of
#'   basis functions/knots for the age-dependent mean and standard deviation of
#'   the latent trait. Defaults are 10 and 8.
#' @param k_sd2 Only used if \code{age_model = "splines"} and \code{irt_model}
#'   is an HO-P or HO-F model. Number of basis functions/knots for the
#'   additional age-dependent component (the lambda term for HO-P, the loading
#'   of the unique age effect for HO-F). Default is 8.
#' @param fixed_knots_mean,fixed_knots_sd,fixed_knots_sd2 Optional. By default
#'   knots are placed based on the data, as in \code{mgcv::s()}. Use
#'   \code{"equidistant"} for equidistant knots, or a numeric vector of knot
#'   locations (on the standardized age scale) for fully user-specified
#'   knots. They apply to the mean, standard deviation and (for HO-P/HO-F
#'   models) the additional component, respectively.
#' @param prior_knowledge Optional. A named list used to supply informative
#'   priors. Each element is named after a parameter group (e.g.,
#'   \code{"alpha"}, \code{"beta"}, \code{"gamma"}, \code{"nu"}, \code{"b"})
#'   and is itself a list of two vectors: the prior means first, then the prior
#'   standard deviations. The names of these two inner elements (e.g.,
#'   \code{alpha_prior_mu}, \code{alpha_prior_sd}) are used as the names in
#'   the returned list, so they must match what the model expects.
#'   \code{alpha}, \code{beta} and \code{gamma} need one value per item;
#'   \code{nu} needs one value per domain in the Testlet/HO models. For
#'   \code{b} in the GRM, each inner element is a list of length \code{I}
#'   containing ordered threshold values (means, then standard deviations).
#'   See the examples.
#' @param parameter_fixed Optional. A named list of parameters to fix. Names
#'   are parameter groups (as for \code{prior_knowledge}) and elements are the
#'   values the parameters are fixed to. A parameter cannot appear in both
#'   \code{prior_knowledge} and \code{parameter_fixed}.
#'
#' @return A named list with the data structures for model fitting.
#' \describe{
#'   \item{Response data}{
#'     \code{I} (number of items), \code{J} (number of persons),
#'     \code{N} (number of observed responses),
#'     \code{ii} and \code{jj} (item and person index of each observed
#'     response), \code{y} (the observed responses, in the same order as
#'     \code{ii} and \code{jj}).}
#'   \item{Age information}{
#'     \code{age_variable} (original age), \code{m_age} and \code{s_age} (mean
#'     and standard deviation used for standardization, rounded to two
#'     decimals), and \code{age_mod} (1 = polynomial, 2 = splines).}
#'   \item{Polynomial age model}{
#'     \code{S} and \code{W} (number of mean and sd coefficients),
#'     \code{age} and \code{age2} (design matrices for the mean and sd), and
#'     for HO-P/HO-F models \code{M} and \code{age3}.}
#'   \item{Spline age model}{
#'     For the mean (index 1), sd (index 2) and, for HO-P/HO-F models, the
#'     additional component (index 3): the number of fixed effects
#'     (\code{S_fixed}, \code{W_fixed}, \code{M_fixed}), the number of
#'     penalized basis functions (\code{K1}, \code{K2}, \code{K3}), the
#'     fixed-effect and penalized design matrices (\code{X1}-\code{X3},
#'     \code{Z1}-\code{Z3}), the associated transformation matrices
#'     (\code{R1_X}-\code{R3_X}, \code{R1_Z}-\code{R3_Z}) and the smooth
#'     objects (\code{sm1}-\code{sm3}) produced when building the basis.}
#'   \item{Model indicator}{\code{irt_mod}: 1 = 1PL, 2 = 2PL, 3 = 3PL,
#'     4 = GRM, 5 = Testlet1PL, 6 = Testlet2PL, 7 = HO1PL-C, 8 = HO2PL-C,
#'     9 = HO1PL-P, 10 = HO2PL-P, 11 = HO1PL-F, 12 = HO2PL-F.}
#'   \item{Multidimensional models}{
#'     \code{itemD}, \code{dd} (domain of each observed response), \code{D}
#'     (number of domains), \code{U} (number of person-domain combinations)
#'     and \code{JD_loadings} (index of the person-domain combination of each
#'     response).}
#'   \item{GRM}{\code{K_max} and \code{K_item}.}
#'   \item{Priors / fixed values}{Any elements supplied via
#'     \code{prior_knowledge} and \code{parameter_fixed}, under their given
#'     names.}
#' }
#'
#' @details
#' Each row of \code{raw_data} is one person and each column in
#' \code{item_variables} is one item (dichotomous, or polytomous for
#' \code{"GRMnorm"}). Missing responses are removed from the data. The
#' remaining responses are stored in long format, with \code{ii}, \code{jj}
#' mapping each response to its item and person.
#'
#' Age is standardized as \code{(age - m_age) / s_age}, where \code{m_age}
#' and \code{s_age} are the mean and standard deviation rounded to two
#' decimals. Polynomial or spline bases are built on the standardized age.
#' Keep \code{m_age} and \code{s_age} if you want to evaluate the fitted
#' model at new ages, since they must be applied in the same way.
#'
#' @seealso \code{\link[mgcv]{s}} for the default knot placement of the spline
#'   bases.
#'
#' @examples
#' ## Simulate example data
#' set.seed(1234)
#' n_person <- 200
#' n_item <- 5
#' example_responses <- matrix(rbinom(n_person * n_item, size = 1, prob = 0.5),
#'                             ncol = n_item)
#' example_data <- data.frame(example_responses,
#'                            age_years = runif(n_person, 6, 18))
#'
#' ## Polynomial age model
#' normdata_poly <- data_prep(
#'   raw_data = example_data,
#'   age_variable = "age_years",
#'   item_variables = paste0("X", 1:5),
#'   irt_model = "2PLnorm",
#'   age_model = "poly",
#'   poly_mean = 2,
#'   poly_sd = 1
#' )
#'
#' ## Inspect the prepared data
#' normdata_poly$N
#' dim(normdata_poly$age)
#' normdata_poly$irt_mod
#'
#' ## Spline age model with prior knowledge on item difficulty and
#' ## discrimination
#' normdata_spline <- data_prep(
#'   raw_data = example_data,
#'   age_variable = "age_years",
#'   item_variables = paste0("X", 1:5),
#'   irt_model = "2PLnorm",
#'   age_model = "splines",
#'   prior_knowledge = list(
#'     alpha = list(alpha_prior_mu = c(1, 1, 0.9, 1.1, 0.8),
#'                  alpha_prior_sd = c(0.5, 0.4, 0.6, 0.5, 0.3)),
#'     beta  = list(beta_prior_mu = c(2, -1, 0, 1.5, -0.5),
#'                  beta_prior_sd = c(0.3, 0.9, 1, 0.2, 0.9))
#'   )
#' )
#' @export
data_prep <- function(raw_data, age_variable, item_variables, 
                      irt_model = c("1PLnorm","2PLnorm","3PLnorm","GRMnorm",
                                    "Testlet1PLnorm","Testlet2PLnorm",
                                    "HO1PLnorm-C", "HO2PLnorm-C",
                                    "HO1PLnorm-P", "HO2PLnorm-P",
                                    "HO1PLnorm-F", "HO2PLnorm-F"),
                      age_model  = c("splines","poly"),
                      itemD = NULL, K_item = NULL,
                      poly_mean = NULL, poly_sd = NULL, poly_sd_unique = NULL,
                      k_mean = 10, k_sd = 8,k_sd2 = 8, fixed_knots_mean = NULL, fixed_knots_sd = NULL,fixed_knots_sd2 = NULL,
                      prior_knowledge = NULL, parameter_fixed = NULL){

  age_model <- match.arg(age_model)
  irt_model <- match.arg(irt_model)

  age <- as.vector(raw_data[, age_variable])
  
  if(any(is.na(age))){
    warning("In your dataset are observations with missing age_variable, 
             these individuals been removed from the dataset.")
    raw_data <- raw_data[!is.na(age),]
    age <- as.vector(raw_data[, age_variable])
  }
  answerpattern <- raw_data[, item_variables]
  nperson <- nrow(raw_data)
  nitem <- length(item_variables)
  jj_list <- lapply(1:nitem, function(i) which(!is.na(answerpattern[,i])))

  R <- as.vector(as.matrix(answerpattern))
  m_age <- round(mean(age),2)
  s_age <- round(stats::sd(age),2)
  age_scale <- (age - m_age)/s_age
  age_df <- data.frame(age = age_scale)
  if(age_model == "poly"){
    if(is.null(poly_sd) | is.null(poly_mean)){
      stop("If age_model = poly is selected poly_mean and poly_sd have to be provided.")
    }
    if(poly_sd == 0){
      age_sd <- matrix(rep(1, length(age_scale)), ncol = 1)
    }else{
      age_sd <- as.matrix(stats::model.matrix(~ poly(age,max(poly_sd), raw = TRUE), data = data.frame(age = age_scale)))
    }
    if(poly_mean == 0){
      age_mean <- matrix(rep(1, length(age_scale)), ncol = 1)
    }else{
      age_mean <- as.matrix(stats::model.matrix(~ poly(age,max(poly_mean), raw = TRUE), data = data.frame(age = age_scale)))
    }
    age_add  <- list( S = ncol(age_mean),
                      W = ncol(age_sd),
                      age = age_mean,
                      age2 = age_sd,
                      age_mod = 1)
    if(irt_model %in% c("HO1PLnorm-P", "HO2PLnorm-P","HO1PLnorm-F", "HO2PLnorm-F")){
      if(is.null(poly_sd_unique)){
        age_sd_u <-  age_sd
      }
      else{
        age_sd_u <- as.matrix(stats::model.matrix(~ poly(age,max(poly_sd_unique), raw = TRUE), data = data.frame(age = age_scale)))
      }
      age_add <- list( S = ncol(age_mean),
                       W = ncol(age_sd),
                       age = age_mean,
                       age2 = age_sd,
                       age_mod = 1, 
                       M = ncol(age_sd_u), 
                       age3 = age_sd_u)
                      
    }
  }
  if(age_model == "splines"){

    mean_sp <- build_spline(age_df, k_mean, fixed_knots_mean)
    sd_sp   <- build_spline(age_df, k_sd,   fixed_knots_sd)
    
    age_add <- list(
      S_fixed = mean_sp$S_fixed,
      K1      = mean_sp$K,
      W_fixed = sd_sp$S_fixed,
      K2      = sd_sp$K,
      X1      = mean_sp$X,
      Z1      = mean_sp$Z,
      X2      = sd_sp$X,
      Z2      = sd_sp$Z,
      R1_Z    = mean_sp$R_Z,
      R1_X    = mean_sp$R_X,
      R2_Z    = sd_sp$R_Z,
      R2_X    = sd_sp$R_X,
      age_mod = 2,
      sm1     = mean_sp$sm,
      sm2     = sd_sp$sm)

  if(irt_model %in% c("HO1PLnorm-P", "HO2PLnorm-P","HO1PLnorm-F", "HO2PLnorm-F")){
    sd_sp2   <- build_spline(age_df, k_sd2,  fixed_knots_sd2)
    age_add[["M_fixed"]] <- sd_sp2$S_fixed
    age_add[["K3"]] <- sd_sp2$K
    age_add[["X3"]] <- sd_sp2$X
    age_add[["Z3"]] <- sd_sp2$Z
    age_add[["R3_Z"]] <- sd_sp2$R_Z
    age_add[["R3_X"]] <- sd_sp2$R_X
    age_add[["sm3"]] <- sd_sp2$sm

  }
  }
 # if(length(fixed_item) > 1 & model != "HO2PLnorm" & model != "HO1PLnorm"){
  #  warning("There are more items indicated to fix then necessary. For the model only the first
 #               item will be fixed")
#    fixed_item <- fixed_item[1]
#  }
  normingdata_IRT <- c(list(I = nitem,
                          J = nperson,
                          N = sum(colSums(!is.na(answerpattern))), #nperson*nitem,
                          ii = rep(1:nitem, times = colSums(is.na(answerpattern) == F)),
                          jj = unlist(jj_list),
                          y = R[!is.na(R)],
                          m_age = m_age,
                          s_age = s_age,
                          itemD = rep(1, nitem),
                          dd =  rep(1,sum(colSums(!is.na(answerpattern)))),
                          age_variable = age), age_add)

 # if(!is.null(fixed_item)){
  #  normingdata_IRT <-  c(normingdata_IRT, list(fixed_beta = fixed_item,
 #                      free_beta = (1:nitem)[-fixed_item],
  #                     fixed_alpha = fixed_item,
  #                     free_alpha = (1:nitem)[-fixed_item]))
  #}

  if(irt_model == "1PLnorm"){
    normingdata_IRT[["irt_mod"]] <- 1
  }
  if(irt_model == "2PLnorm"){
    normingdata_IRT[["irt_mod"]] <- 2
  }
  if(irt_model == "3PLnorm"){
    normingdata_IRT[["irt_mod"]] <- 3
  }

  if(any(c("Testlet1PLnorm","Testlet2PLnorm", "HO1PLnorm-C", "HO2PLnorm-C",
           "HO1PLnorm-P", "HO2PLnorm-P", "HO1PLnorm-F", "HO2PLnorm-F") == irt_model)){
    if(is.null(itemD)){
      stop("itemD missing. Please specify the dimensions each item is loading on.")
    }
    if(length(itemD) != length(item_variables)){
      stop("length of itemD and item_variables differ. itemD should indicate the domain each
           item in item_variables measures.")
    }
    D <- length(unique(itemD))
    #if(model == "HO2PLnorm" | model == "HO1PLnorm"){
    #  if(length(unique(itemD[fixed_item])) != D){
     #   stop("fixed_item is not sufficent. Please fix at least one item per domain for identification.")
     # }
   # }
    idx_list <- lapply(seq_len(nitem), function(i) {
      which(!is.na(answerpattern[, i])) + (itemD[i] - 1) * nperson
    })

    JD <- unlist(idx_list)
    dd_data <- rep(itemD, colSums(!is.na(answerpattern)))
    normingdata_IRT[["dd"]] <- dd_data
    normingdata_IRT[["D"]] <- D
    normingdata_IRT[["U"]] <- D*nperson
    normingdata_IRT[["JD_loadings"]] <- JD
    normingdata_IRT[["irt_mod"]] <- 5
    normingdata_IRT[["itemD"]] <- itemD

    if(irt_model == "Testlet2PLnorm"){
      normingdata_IRT[["irt_mod"]] <-  6
    }
    if(irt_model == "HO1PLnorm-C"){
      normingdata_IRT[["irt_mod"]] <-  7
    }
    if(irt_model == "HO2PLnorm-C"){
      normingdata_IRT[["irt_mod"]] <-  8
    }
    if(irt_model == "HO1PLnorm-P"){
      normingdata_IRT[["irt_mod"]] <-  9
    }
    if(irt_model == "HO2PLnorm-P"){
      normingdata_IRT[["irt_mod"]] <-  10
    }
    if(irt_model == "HO1PLnorm-F"){
      normingdata_IRT[["irt_mod"]] <-  11
    }
    if(irt_model == "HO2PLnorm-F"){
      normingdata_IRT[["irt_mod"]] <-  12
    }
  }
  
  if(irt_model == "GRMnorm"){
    if(is.null(K_item)){
      stop("K_item missing. Please specify the number of response categories per item.")
    }
    if(length(K_item) != length(item_variables)){
      stop("length of K_item and item_variables differ. K_item should indicate the number of response categories per
           item in item_variables.")
    }
    if(min(normingdata_IRT$y) < 1){
      normingdata_IRT$y <- normingdata_IRT$y + abs(min(normingdata_IRT$y)) +1
      warning("The smallest item response value was below 1. For the GRMnorm model item responses need to be coded positivly. Therefore responses been recoded by adding the smallest value to all responses.")
    }
    normingdata_IRT[["K_max"]] <- max(K_item)
    normingdata_IRT[["K_item"]] <- K_item
    normingdata_IRT[["irt_mod"]] <- 4
  }
  
  max_score <- if(irt_model == "GRMnorm") K_item else rep(1, nitem)
  if(any(normingdata_IRT$y > max_score[normingdata_IRT$ii])){
    stop("There are observations with higher item scores than allowed ",
         "(1 for binary models, K_item for GRMnorm). ",
         "Resolve this issue and run the function again.")
  }
  prior_names <- names(prior_knowledge)
  para_names <- names(parameter_fixed)
  info_names <- c(para_names, prior_names)
  if(any(duplicated(info_names))){
    stop(paste0("The model parameters",info_names[duplicated(info_names)],
                "are indicated with prior information and as fixed, you can
                  only indicate one."))
  }
  
  if(!is.null(prior_knowledge)){
    for(j in seq_len(length(prior_knowledge))){
      if(any(prior_names[j] == c("alpha","beta","gamma")) &
         #length(prior_knowledge[[j]]) != 2 &
         any(lengths(prior_knowledge[[j]]) != nitem)){
        stop(paste0("When you choose to incoperate prior knowledge for ", prior_names[j],
                    " you have to indicate the distribution mean and standard deviation for all items.
                        Your information is in the wrong format."))
      }
      if(prior_names[j] == "nu" & any(c("Testlet1PLnorm","Testlet2PLnorm", "HO1PLnorm-C", "HO2PLnorm-C",
                                           "HO1PLnorm-P", "HO2PLnorm-P") == irt_model)){
        if(any(lengths(prior_knowledge[[j]]) != D)){
        stop(paste0("When you choose to incoperate prior knowledge for ", prior_names[j],
                    " you have to indicate the distribution mean and standard deviation for all D domains.
                        Your information is in the wrong format."))
        }
      }
      if(irt_model == "GRMnorm" & prior_names[j] == "b"){
        if(length(prior_knowledge[[j]][[1]])!= nitem |
           length(prior_knowledge[[j]][[2]])!= nitem |
           any(lengths(prior_knowledge[[j]][[1]]) != max(K_item) -1) |
           any(lengths(prior_knowledge[[j]][[2]]) != max(K_item) -1) |
           any(sapply(prior_knowledge[[j]][[1]], function(x) any(diff(x) <0)))){
          stop(paste0("When you choose to incoperate prior knowledge for ", prior_names[j],
                      " you have to indicate those for all thresholds. Therefore you need
                        to indicate a list of length I with each list element having ordered values of length K_max.
                        Your information is in the wrong format."))
        }
      }
      
      
      normingdata_IRT[[names(prior_knowledge[[j]])[1]]] <- prior_knowledge[[j]][[1]]
      normingdata_IRT[[names(prior_knowledge[[j]])[2]]] <- prior_knowledge[[j]][[2]]
      
    }
  }
  
  if(!is.null(parameter_fixed)){
    for(j in seq_len(length(parameter_fixed))){
      if(para_names[j] %in% c("alpha","beta","gamma") & length(parameter_fixed[[j]]) != nitem){
        stop(paste0("When you choose to fix ", para_names[j],
                    " you have to indicate those for all items."))
      }
      #  if(para_names[j] == "nu" & length(parameter_fixed[[j]]) != D)){
      #      stop(paste0("When you choose to fix ", prior_names[j],
      #                  " you have to indicate those for all D domains.")
      #    }
        if(irt_model == "GRMnorm" & para_names[j] == "b"){
          if(length(parameter_fixed[[j]]) != nitem |
             !is.list(parameter_fixed[[j]]) |
             all(lengths(parameter_fixed[[j]]) != max(K_item) -1 ) |
             any(sapply(parameter_fixed[[j]], function(x) any(diff(x) <0)))){
               stop(paste0("When you choose to fix ", para_names[j],
                           " you have to indicate those for all thresholds. Therefore you need
                        to indicate a list of length I with each list element having ordered values of length K_max."))
             }
        }
           
      
      normingdata_IRT[[para_names[j]]] <- parameter_fixed[[j]]
      
    }
  }
  

  return(normingdata_IRT)
}


check_rawinput <- function(raw_data, item_variables, age_variable){
  if(is.numeric(item_variables)){
    if(max(item_variables) > ncol(raw_data) | min(item_variables) < 1){
      stop("`item_variables` contains an index outside the range of `raw_data`.")
    }
    item_variables <- colnames(raw_data)[item_variables]
  }
  if(is.numeric(age_variable)){
    if(age_variable > ncol(raw_data) | age_variable < 1){
      stop("`age_variable` is outside the range of `raw_data`.")
    }
    age_variable <- colnames(raw_data)[age_variable]
  }
  # Input checks
  missing_cols <- setdiff(c(age_variable, item_variables), names(raw_data))
  if (length(missing_cols) > 0) {
    stop("Columns not found in raw_data: ", paste(missing_cols, collapse = ", "))
  }
  
}

