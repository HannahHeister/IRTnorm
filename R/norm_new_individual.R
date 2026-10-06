# compute age-effects from posterior draws
compute_age_effects <- function(post_draw, data, age_scale, nperson, irt_mod, multidim) {
  
  # calculating age effects for age_new
  if(data$age_mod == 1){
    p1_theta_post <- post_draw$p1_theta
    p2_theta_post <- post_draw$p2_theta
    if(data$S == 1){
      age_mean <- matrix(rep(1,length(age_scale)), ncol = 1)
    }else{
      age_mean <- as.matrix(stats::model.matrix(~ poly(age,data$S - 1, raw = TRUE), data = data.frame(age = age_scale)))
    }
    if(data$W == 1){
      age_sd <- matrix(rep(1,length(age_scale)), ncol = 1)
    }else{
      age_sd <- as.matrix(stats::model.matrix(~ poly(age,data$W - 1, raw = TRUE), data = data.frame(age = age_scale)))
    }
        if(multidim == FALSE){
      D <- 1
      mu_post <- p1_theta_post %*% t(age_mean)
      sigma_post <- exp(p2_theta_post %*% t(age_sd))
    }else{
      D <- data$D
      mu_post <- array(dim = c(nrow(p1_theta_post), nperson, D))
      sigma_post <-array(dim = c(nrow(p1_theta_post), nperson, D))
      if(irt_mod == 5 |irt_mod == 6){
        for(d in seq_len(D)){
          mu_post[,,d] <- p1_theta_post %*% t(age_mean)
          sigma_post[,,d] <- exp(p2_theta_post %*% t(age_sd))
        }
      }
      if(irt_mod >= 7){
        for(d in seq_len(D)){
          mu_post[,,d] <- p1_theta_post[,seq(d,D*data$S,D)] %*% t(age_mean)
          sigma_post[,,d] <- exp(p2_theta_post[,seq(d,D*data$W,D)] %*% t(age_sd))
        }
        if(irt_mod >= 9){
          if(data$M == 1){
            age3 <- matrix(rep(1,length(age_scale)), ncol = 1)
            
          }else{
            age3 <- as.matrix(stats::model.matrix(~ poly(age,data$M - 1, raw = TRUE), data = data.frame(age = age_scale)))
          }
          p3_theta_post <- post_draw$p3_theta
        }
        if(irt_mod == 9 | irt_mod == 10){
          lambda_post <- p3_theta_post %*% t(age3)
        }
        if(irt_mod == 11 | irt_mod == 12){
          unique_post <- array(dim = c(nrow(p3_theta_post), nperson, D))
          for(d in seq_len(D)){
            unique_post[,,d] <- log1p(exp(p3_theta_post[,seq(d,D*data$M,D)] %*% t(age3)))
          }
        }
      }
      
    }
  }
  if(data$age_mod == 2){
    l1_post <- post_draw$l1
    l2_post <- post_draw$l2
    u1_post <- post_draw$u1
    u2_post <- post_draw$u2
    sm1 <- data$sm1
    sm2 <- data$sm2
    R1_X <- data$R1_X
    R2_X <- data$R2_X
    R1_Z <- data$R1_Z
    R2_Z <- data$R2_Z
    data_temp <- data.frame(age = c(age_scale))
    Xraw1_new <- mgcv::PredictMat(sm1, data_temp)
    X1_new    <- Xraw1_new %*% R1_X          # correct unpenalized
    Z1_new    <- Xraw1_new %*% R1_Z          # correct penalized
    
    Xraw2_new <- mgcv::PredictMat(sm2, data_temp)
    X2_new    <- Xraw2_new %*% R2_X          # correct unpenalized
    Z2_new    <- Xraw2_new %*% R2_Z          # correct penalized
    if(multidim == FALSE){
      D <- 1
      mu_post <- l1_post %*% t(X1_new) + u1_post %*% t(Z1_new)
      sigma_post <- exp(l2_post %*% t(X2_new) + u2_post %*% t(Z2_new))
    }else{
      D <- data$D
      mu_post <- array(dim = c(nrow(l1_post), nperson, D))
      sigma_post <-array(dim = c(nrow(l2_post), nperson, D))
      if(irt_mod == 5 |irt_mod == 6){
        for(d in seq_len(D)){
          mu_post[,,d] <- l1_post %*% t(X1_new) + u1_post %*% t(Z1_new)
          sigma_post[,,d] <- exp(l2_post %*% t(X2_new) + u2_post %*% t(Z2_new))
        }
      }
      if(irt_mod >=7){
        for(d in seq_len(D)){
          mu_post[,,d] <- l1_post[,seq(d,D*data$S,D)] %*% t(X1_new) + u1_post[,seq(d,D*data$S,D)]%*% t(Z1_new)
          sigma_post[,,d] <- exp(l2_post[,seq(d,D*data$W,D)] %*% t(X2_new)+ u2_post[,seq(d,D*data$W,D)]) %*% t(Z2_new)
        }
        if(irt_mod >= 9){
          l3_post <- post_draw$l3
          u3_post <- post_draw$u3
          sm3 <- data$sm3
          R3_X <- data$R3_X
          R3_Z <- data$R3_Z
          Xraw3_new <- mgcv::PredictMat(sm3, data_temp)
          X3_new    <- Xraw3_new %*% R3_X          # correct unpenalized
          Z3_new    <- Xraw3_new %*% R3_Z          # correct penalized
        }
        if(irt_mod == 9 | irt_mod == 10){
          lambda_post <- l3_post %*% t(X3_new) + u3_post %*% t(Z3_new)
        }
        if(irt_mod == 11 | irt_mod == 12){
          unique_post <- array(dim = c(nrow(l3_post), nperson, D))
          for(d in seq_len(D)){
            unique_post[,,d] <- log1p(exp( l3_post[,seq(d,D*data$M,D)] %*% t(X3_new) + u3_post[,seq(d,D*data$M,D)] %*% t(Z3_new)))
          }
        }
      }
      
    }
  }
  if(!irt_mod %in% c(9,10)){
    lambda_post <- NULL
  }
  if(!irt_mod %in% c(11,12)){
    unique_post <- NULL
  }
  
  return(list(mu_post = mu_post, sigma_post = sigma_post, 
              lambda_post = lambda_post, unique_post = unique_post))
  
}

# log posterior for one MH evaluation
# Extract the existing log_post_draw_theta inner function to the top level
# Helper: log post_draw for theta given one draw of params (vectorized theta)
log_post_draw_theta <- function(data, theta, a_vec, b_vec, c_vec,  y_vec, mu, sigma,
                                unique =NULL, nu_vec = NULL, lambda_vec = NULL,
                                unique_vec = NULL, irt_mod, itemD = NULL) {
  nitem <- length(y_vec)
  if(irt_mod > 4){
    D <- data$D
  }else{
    D <- 1
  } 
  if(is.null(c_vec) | all(c_vec == 0)){
    c_vec <- rep(0, nitem)
  }
  if(is.null(a_vec)| all(a_vec == 1)){
    a_vec <- rep(1,nitem)
  }
  # compute log prior: normal(mu, sigma)
  #unidim 
  if(irt_mod <= 4){
    lp_prior <- stats::dnorm(theta, mean = mu, sd = sigma, log = TRUE)
    # compute log-likelihood row by row
    if(irt_mod != 4){
      p <- c_vec + (1-c_vec)*stats::plogis(a_vec * (theta - b_vec))
      ll <- sum( stats::dbinom(as.numeric(y_vec), size = 1, prob = p, log = TRUE))
    }
    if(irt_mod == 4){
      y_mat <- make_dummy_matrix(info = y_vec, K_item = data$K_item)
      all_p <- NULL
      for (i in seq_len(nitem)) {
        #int_threshold <-((1 + (data$K_max -1)*(i-1)):((data$K_max -1)*(i)))[seq_len(data$K_item[i] -1)]
        p <- grm_probs(theta = theta, a = a_vec[i], b = b_vec[[i]][seq_len(data$K_item[i] -1)]) # b_vec[int_threshold])
        all_p <- c(all_p, p)
      }
      ll<- sum(y_mat * log(all_p))
      
    }
    
  }
  #multidim
  if(irt_mod > 4){
    lp_prior <-  stats::dnorm(theta, mean = 0, sd = 1, log = TRUE) + 
      sum( stats::dnorm(unique, mean = 0, sd = 1, log = TRUE))
    if(irt_mod >4 & irt_mod < 9){
      delta <- mu + nu_vec * sigma * theta + sqrt(1 - nu_vec^2)*sigma *unique
    }
    if(irt_mod == 9 | irt_mod == 10){
      delta <- mu + nu_vec * lambda_vec * sigma * theta + sqrt(1 - (nu_vec *lambda_vec)^2 )*sigma *unique
    }
    if(irt_mod == 11 | irt_mod == 12){
      delta <- mu + sigma * theta + unique_vec *unique
    }
    
    p <- c_vec + (1-c_vec)*stats::plogis(a_vec * (delta[itemD] - b_vec))
    ll <- sum(stats::dbinom(y_vec, size = 1, prob = p, log = TRUE))
    
  }
  
  lp <- lp_prior + ll
  return(lp)
}


# Fully Bayesian estimation of latent trait for a new observation
# given post_draw draws of 2PL item parameters and latent distribution params.


#'@title Estimate norm score for new observations
#'@name norm_new_individual
#'@description Performs Bayesian estimates of the norm score for a new individual using
#' Metropolis-Hastings sampling. This function integrates over
#' posterior uncertainty in item parameters and norming parameters from a
#' previously fitted IRT model.
#' @param new_y Numeric vector of binary item responses (0/1) for a single individual.
#'   Can contain NAs for missing responses.
#' @param new_age Numeric scalar. The age of the individual for whom theta is being estimated.
#' @param fit_info A named \code{list} as returned by \code{\link{estimate_norming_model}}.
#' @param mh_args List of Metropolis-Hastings algorithm tuning parameters:
#'   \itemize{
#'     \item \code{n_iter}: Total number of MCMC iterations (default: 250).
#'     \item \code{burnin}: Number of initial iterations to discard (default: 125).
#'     \item \code{thin}: Thinning interval - keep every nth sample (default: 1).
#'     \item \code{proposal_sd}: Standard deviation of proposal distribution (default: 0.6).
#'   }
#'
#' @return A list with the following components:
#'   \item{theta_samples}{Numeric vector of all posterior samples of theta (post-burnin, post-thinning)}.
#'   \item{theta_zsample}{Numeric vector of standardized theta samples}.
#'   \item{by_draw_samples}{List of length M, each element containing theta samples for one posterior draw}.
#'   \item{acceptance_rates}{Numeric vector of length M with acceptance rates for each posterior draw}.
#'   \item{summary}{List containing: mean, sd, median, and 95% credible interval of theta}.
#'   \item{mh_args}{List of MH arguments used}.
#'   \item{obs_items}{Integer vector of indices for observed (non-NA) items}.
#'   \item{mean_mu}{Mean of the posterior distribution for mu (theta population mean)}.
#'   \item{mean_sigma}{Mean of the posterior distribution for sigma (theta population SD)}. 
#'
#' @details
#' This function runs a Metropolis-Hastings algorithm to sample theta conditional on the
#'   item parameters (a, b) and population parameters (mu, sigma) draws of the fitted model.
#'
#' Only observed items (non-NA responses) are used in the likelihood calculation.
#'
#' @section Tuning the MH Algorithm:
#' The \code{proposal_sd} parameter controls the random-walk step size. Target acceptance
#' rates are typically 20-50%. If acceptance rates are too high (>70%), increase
#' \code{proposal_sd}. If too low (<15%), decrease it.
#'
#'@seealso \code{\link{info_new_norming}}
#' @examples
#' \dontrun{
#' # Load data
#' data(response_data)
#'
#' # Fit IRTnorm model and extract important information 
#' info <- estimate_norming_model(
#'   raw_data        = response_data,
#'   age_variable    = "age",
#'   item_variables   = 1:50,
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
#' # new individual has same answer pattern as individual 20 in response_data
#' y <- unlist(response_data[20,1:50])
#' age <- response_data$age[20]
#'
#' # Estimate theta for the new individual
#' result <- estimate_theta_new(new_y = y,
#'                              new_age = age,
#'                              fit_info = info)
#'
#' EAP estimate of the individuals normscore
#' mean(result$theta_zsample)
#'
#' # Check acceptance rates
#' summary(result$acceptance_rates)
#' }
#' @export
norm_new_individual <- function(new_y, new_age,fit_info, mh_args = list(n_iter=250, burnin=125, thin=1, proposal_sd=0.6)) {
  
  data <- fit_info$newdata
  post_draw <- fit_info$norm_info
  
  if(data$irt_mod > 4){
    D <- data$D
  }else{
    D <- 1
  } 
  
  ## extract post_draw draws
  if(data$irt_mod == 4){
    a_post <- post_draw$alpha  # post_draw$a 
    b_post <- post_draw$threshold
    c_post <- NULL
    nu_post <- NULL
  }else{
    a_post <- post_draw$alpha  # post_draw$a 
    b_post <- post_draw$beta 
    c_post <- post_draw$gamma
    nu_post <- post_draw$nu
  }

  
  age_scale <- (new_age - data$m_age)/data$s_age
  nperson <-  length(new_age)
  multidim <- length(unique(data$itemD)) != 1
  irt_mod <- data$irt_mod

  ageeffects <- compute_age_effects(post_draw = post_draw, data = data, 
                                    age_scale = age_scale, nperson = nperson,
                                    irt_mod = irt_mod, multidim = multidim)
  mu_post <- ageeffects$mu_post
  sigma_post <- ageeffects$sigma_post
  lambda_post <- ageeffects$lambda_post
  unique_post <- ageeffects$unique_post
  
  M <- nrow(b_post)
  # Handle NAs in y: only use observed items
  obs_idx <- which(!is.na(new_y))
  y_obs <- new_y[obs_idx]
  
  if(is.data.frame(y_obs)){
    y_obs <- as.matrix(y_obs)
  }
  if(irt_mod == 4){
    obs_idx_GRM <- NULL
    for(i in  which(!is.na(new_y))){
      obs_idx_GRM <- c(obs_idx_GRM, which(rep(1:(data$I), data$K_max -1) == i))
      
    }
  }else{
    obs_idx_GRM <- obs_idx
  }
  
  # MH arguments
  n_iter <- mh_args$n_iter
  burnin  <- mh_args$burnin
  thin    <- mh_args$thin
  proposal_sd <- mh_args$proposal_sd
  
  
  # Pre-allocate
  all_ztheta_samples <- numeric(0)
  all_theta_samples <- numeric(0)
  all_unique_samples <- NULL
  all_zdelta_samples <- NULL
  all_delta_samples<- NULL
  by_draw_samples <- vector("list", M)
  acceptance_rates <- numeric(M)
  
  itemD_obs <- data$itemD[obs_idx]
  
  # Loop over post_draw draws
  for (m in seq_len(M)) {
    a_m <- a_post[m, obs_idx]
    b_m <- b_post[m, obs_idx_GRM]
    c_m <- c_post[m, obs_idx]
    
    # init
    accepts <- 0L
    if(!multidim){
      mu_m <- mu_post[m,]
      sigma_m <- sigma_post[m,]
      
      chain <- numeric(n_iter)
      
      theta_curr <- mu_m  # start at population mean
      lp_curr <- log_post_draw_theta(data = data, theta = theta_curr, a_vec = a_m, b_vec = b_m,
                                     c_vec = c_m,
                                     y_vec = y_obs, mu = mu_m, sigma = sigma_m, 
                                     irt_mod = irt_mod)
      
      for (it in seq_len(n_iter)) {
        theta_prop <- stats::rnorm(1, mean = theta_curr, sd = proposal_sd)
        lp_prop <- log_post_draw_theta(data = data, theta_prop, a_m, b_m, c_m,
                                       y_obs, mu_m, sigma_m,irt_mod = irt_mod)
        log_accept_ratio <- lp_prop - lp_curr
        if (is.finite(log_accept_ratio) && log(stats::runif(1)) < log_accept_ratio) {
          theta_curr <- theta_prop
          lp_curr <- lp_prop
          accepts <- accepts + 1L
        }
        chain[it] <- theta_curr
      }
      acceptance_rates[m] <- accepts / n_iter
      # apply burnin and thin
      keep_idx <- seq(from = burnin + 1, to = n_iter, by = thin)
      samples_m <- chain[keep_idx]
      by_draw_samples[[m]] <- samples_m
      all_theta_samples <- c(all_theta_samples, samples_m)
      all_ztheta_samples <- c(all_ztheta_samples, (samples_m - mean(mu_m))/mean(sigma_m))
      
    }
    if(multidim){
      mu_m <- mu_post[m,,]
      sigma_m <- sigma_post[m,,]
      unique_m <- unique_post[m,,]
      nu_m <- as.numeric(nu_post[m, ])
      lambda_m <- lambda_post[m,]
      
      chain <- matrix(nrow = n_iter, ncol = D*2 +1)
      
      theta_curr <- rep(0, 1)  # general latent trait then unique
      unique_curr <- rep(0, D) 
      lp_curr <- log_post_draw_theta(data = data, theta = theta_curr, unique = unique_curr, 
                                     a_vec = a_m, b_vec = b_m, c_vec = c_m,
                                     y_vec = y_obs, mu = mu_m, sigma = sigma_m,
                                     nu_vec = nu_m, lambda_vec = lambda_m,
                                     unique_vec = unique_m, irt_mod = irt_mod,
                                     itemD =  itemD_obs)
      
      for (it in seq_len(n_iter)) {
        theta_prop <- stats::rnorm(1, mean = theta_curr, sd = proposal_sd)
        unique_prop <- stats::rnorm(D, mean = unique_curr, sd = rep(proposal_sd,D))
        for(d in seq_len(D + 1)){
          unique_temp <- unique_curr
          if(d == 1){
            theta_temp <- theta_prop
          }else{
            theta_temp <- theta_curr
            unique_temp[d -1] <- unique_prop[d-1]
          }
          lp_prop <- log_post_draw_theta(data = data, theta = theta_temp,unique = unique_temp, 
                                         a_m, b_m, c_m, y_obs, mu_m, sigma_m,
                                         nu_vec = nu_m, lambda_vec = lambda_m,
                                         unique_vec = unique_m, irt_mod = irt_mod,
                                         itemD =  itemD_obs)
          log_accept_ratio <- lp_prop - lp_curr
          if (is.finite(log_accept_ratio) && log(stats::runif(1)) < log_accept_ratio) {
            theta_curr <- theta_temp
            unique_curr <- unique_temp
            lp_curr <- lp_prop
          }
        }
        if(irt_mod < 9){
          delta_curr <- mu_m + nu_m * sigma_m * theta_curr + sqrt(1 - nu_m^2)*sigma_m *unique_curr
          z_delta_curr <- (delta_curr - mu_m)/sigma_m
        }
        if(irt_mod == 9 | irt_mod == 10){
          delta_curr <- mu_m + nu_m * lambda_m * sigma_m * theta_curr+ sqrt(1 - (nu_m *lambda_m)^2 )*sigma_m *unique_curr
          z_delta_curr <- (delta_curr - mu_m)/sigma_m
        }
        if(irt_mod == 11 | irt_mod == 12){
          delta_curr <- mu_m + sigma_m * theta_curr + unique_m *unique_curr
          z_delta_curr <- (delta_curr - mu_m)/sqrt(sigma_m^2 + unique_curr^2)
        }
        chain[it,] <- c(theta_curr, z_delta_curr, delta_curr)
      }
      # apply burnin and thin
      keep_idx <- seq(from = burnin + 1, to = n_iter, by = thin)
      samples_m <- chain[keep_idx, ]
      by_draw_samples[[m]] <- samples_m
      all_theta_samples <- c(all_theta_samples, samples_m[,1])
      all_zdelta_samples  <- rbind(all_zdelta_samples, samples_m[,2:(D + 1)])
      all_delta_samples <- rbind(all_delta_samples, samples_m[,(D +2):(D*2 + 1)])
    }
  }
  
  if(!multidim){
    # Summaries
    ztheta_mean <- mean(all_ztheta_samples)
    ztheta_sd <- stats::sd(all_ztheta_samples)
    ztheta_med <- stats::median(all_ztheta_samples)
    ci <- stats::quantile(all_ztheta_samples, probs = c(0.025, 0.975))
    
    res <- list(
      theta_samples = all_theta_samples,
      theta_zsample = all_ztheta_samples,
      by_draw_samples = by_draw_samples,
      acceptance_rates = acceptance_rates,
      summary = list(mean = ztheta_mean, sd = ztheta_sd, median = ztheta_med, ci = ci),
      mh_args = mh_args,
      obs_items = obs_idx
    )
  }
  if(multidim){
    # Summaries
    delta_mean <- colMeans(all_zdelta_samples)
    delta_sd <- apply(all_zdelta_samples, 2, stats::sd)
    delta_med <- apply(all_zdelta_samples, 2, stats::median)
    ci <- apply(all_zdelta_samples, 2, function(x){return(stats::quantile(x,  probs = c(0.025, 0.975)))})
    
    theta_mean <- mean(all_theta_samples)
    theta_sd <- stats::sd(all_theta_samples)
    theta_med <- stats::median(all_theta_samples)
    theta_ci <- stats::quantile(all_theta_samples,  probs = c(0.025, 0.975))
    
    res <- list(
      theta_samples = all_theta_samples,
      zdelta_sample = all_zdelta_samples,
      delta_sample = all_delta_samples,
      by_draw_samples = by_draw_samples,
      summary = data.frame(mean = c(theta_mean, delta_mean),
                           sd = c(theta_sd, delta_sd),
                           median = c(theta_med, delta_med),
                           lower2_5 =  c(theta_ci[1], ci[1,]),
                           upper97_5 =  c(theta_ci[2], ci[2,]),
                           domain = c("general", paste0("domain", seq_len(D)))),
      mh_args = mh_args,
      obs_items = obs_idx)
  }
  return(res)
}





