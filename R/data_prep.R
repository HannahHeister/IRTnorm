build_spline <- function(age_df, k, fixed_knots) {
  if (!is.null(fixed_knots)) {
    if (identical(fixed_knots, "equidistant")) {
      fixed_knots <- list(age = seq(min(age_df$age), max(age_df$age), length.out = k))
    } else {
      fixed_knots <- list(age = fixed_knots)
    }
    k <- length(fixed_knots$age)
    
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


#' @title Prepare data for fitting a IRT-based continuous norming model
#'@name data_prep
#'
#'@description The data_prep function prepares the data such that it can be
#'used to norm the data with a Bayesian IRT-based norming model.
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
#'     \item{\code{"splines"}}{Penalized B-spline basis functions of standardized age are
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
#' @param itemD only required if \code{model} is a multidimensional model.
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
#' latent trait mean
#'
#' @param k_sd only required if \code{age_model = "splines"}.
#' A \code{num} indicating the number of knots used for modelling the
#' latent trait standard deviation.
#' 
#' @param k_sd2 only required if \code{age_model = "splines"}  and \code{irt_model} a HO-P or HO-C model.
#' A \code{num} indicating the number of knots used for modelling the lambda for HO-P 
#' and the loading of the unique age effect for HO-F. 
#'
#'@param fixed_knots_mean only required if knots should not be chosen based on
#' data such as the default in \code{s} from \code{mgcv}. It can be choose between
#' equidistant knots by setting fixed_knots_sd =  "equidistant" or fully user
#' specified which a numerical vector of knot values.
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
#' @return A named list containing all data structures required for fitting
#'   the specified IRT norming model.
#'
#' @details
#' Prepare data for fitting an age-dependent IRT-based norming model.
#'
#' The \code{data_prep} function restructures person–item response data into
#' a list format suitable for fitting age-dependent IRT norming models.
#' The function extracts item responses, handles missing data, scales age,
#' and constructs age design matrices according to the specified age model.
#'
#' Each row in \code{raw_data} corresponds to one individual and each column in
#' \code{item_variables} corresponds to one dichotomous test item.

#' Age is internally standardized using its mean and standard deviation prior
#' to constructing polynomial or spline basis functions. Missing item responses
#' are removed and handled via index vectors (\code{ii}, \code{jj}) that map
#' responses to persons and items.
#'
#' @examples
#' ## Simulate example data
#' set.seed(1234)
#'
#' n_person <- 200
#' n_item <- 5
#'
#' example_responses <- matrix(rbinom(n_person*n_item, size = 1, prob = 0.5),
#'  ncol =n_item)
#' example_data <- data.frame(example_responses,
#'                            age_years = runif(n_person, 6, 18))
#'
#' ## Prepare data using a polynomial age model
#' normdata_poly <- data_prep(
#'   raw_data = example_data,
#'   age_variable = "age_years",
#'   item_variables = paste0("X", 1:5),
#'   irt_model = "2PLnorm",
#'   age_model = "poly",
#'   poly_mean = 2,
#'   poly_sd = 1,
#' )
#'
#' ## Prepare data using a spline-based age model with prior knowledge about the
#' ## item difficulty and item discrimination
#' normdata_spline <- data_prep(
#'   raw_data = example_data,
#'   age_variable = "age_years",
#'   item_variables = paste0("X", 1:5),
#'   irt_model = "2PLnorm",
#'   age_model = "splines",
#'   prior_knowledge = list(alpha = list(alpha_prior_mu = c(1,1,0.9,1.1,0.8),
#'                                       alpha_prior_sd = c(0.5, 0.4, 0.6, 0.5, 0.3)),
#'                          beta = list(beta_prior_mu = c(2,-1,0,1.5,-0.5),
#'                                      beta_prior_sd = c(0.3,0.9,1,0.2,0.9)))
#' )
#'@export
#'@keywords internal

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
         length(prior_knowledge[[j]]) != 2 &
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
        stop(paste0("When you choose to fix ", prior_names[j],
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
  
  if(irt_model == "GRMnorm"){
    if(is.null(K_item)){
      stop("K_item missing. Please specify the number of response categories per item.")
    }
    if(length(K_item) != length(item_variables)){
      stop("length of K_item and item_variables differ. K_item should indicate the number of response categories per
           item in item_variables.")
    }
    if(min(normingdata_IRT$y) < 1){
      normingdata_IRT$y <- normingdata_IRT$y + abs(min(normingdata_IRT$y))
      print("The smallest item response value was below 1. For the GRMnorm model item responses need to be coded positivly. Therefore responses been recoded by adding the smallest value to all responses.")
    }
    normingdata_IRT[["K_max"]] <- max(K_item)
    normingdata_IRT[["K_item"]] <- K_item
    normingdata_IRT[["irt_mod"]] <- 4

  }
  return(normingdata_IRT)
}



