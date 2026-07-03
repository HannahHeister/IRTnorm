# Helper functions for testthat 

#' Minimal 1-D data list
make_data_1d <- function(n_persons = 50, n_items = 6) {
  list(
    irt_mod  = 1,
    D      = 1,
    itemD  = rep(1L, n_items),
    dd     = rep(1L, n_items * n_persons),
    jj     = rep(seq_len(n_persons), each = n_items),
    K_item = NULL
  )
}

#' 2-D data list
make_data_multidim <- function(n_persons = 40, n_items_per_dim = 5) {
  n_items <- n_items_per_dim * 2
  list(
    irt_mod  = 5,
    D      = 2,
    itemD  = rep(1:2, each = n_items_per_dim),
    dd     = rep(rep(1:2, each = n_items_per_dim), n_persons),
    jj     = rep(seq_len(n_persons), each = n_items),
    K_item = NULL
  )
}

#' fit_info whose y_rep is all a single constant value
#' fit_info with constant y_rep (all draws = value)
make_fit_info_constant <- function(data, value = 0L, n_draws = 100) {
  list(y_rep = matrix(value, nrow = n_draws, ncol = length(data$jj)))
}
#' fit_info with random binary responses
make_fit_info_random <- function(data, n_draws = 200, seed = 42) {
  n_obs <- length(data$jj)
  set.seed(seed)
  y_rep <- matrix(
    sample(0:1, size = n_draws * n_obs, replace = TRUE),
    nrow = n_draws, ncol = n_obs
  )
  list(y_rep = y_rep)
}

#' Raw data where every person scores the same fixed value on every item
make_raw_data_constant <- function(data, int_variables,
                                   item_value   = 1L,
                                   age_variable = "age",
                                   age_min = 6, age_max = 18) {
  n_persons <- max(data$jj)
  n_items   <- length(int_variables)
  df <- as.data.frame(
    matrix(item_value, nrow = n_persons, ncol = n_items,
           dimnames = list(NULL, int_variables))
  )
  df[[age_variable]] <- seq(age_min, age_max, length.out = n_persons)
  df
}

#' Raw data with random binary responses
make_raw_data_random <- function(data, int_variables,
                                 age_variable = "age",
                                 age_min = 6, age_max = 18,
                                 seed = 7) {
  n_persons <- max(data$jj)
  n_items   <- length(int_variables)
  set.seed(seed)
  df <- as.data.frame(
    matrix(
      sample(0:1, size = n_persons * n_items, replace = TRUE),
      nrow = n_persons, ncol = n_items,
      dimnames = list(NULL, int_variables)
    )
  )
  df[[age_variable]] <- runif(n_persons, age_min, age_max)
  df
}


# settings to evaluate for pipeline tests
irt_models <- c("1PLnorm","2PLnorm","3PLnorm","GRMnorm",
                "Testlet1PLnorm","Testlet2PLnorm",
                "HO1PLnorm-C", "HO2PLnorm-C",
                "HO1PLnorm-P", "HO2PLnorm-P",
                "HO1PLnorm-F", "HO2PLnorm-F")

age_models <- c("splines","poly")

all_combi <- list()
for(mod in irt_models){
  if(mod %in% c("2PLnorm", "Testlet2PLnorm","HO2PLnorm-C", "HO2PLnorm-P", "HO2PLnorm-F", "GRMnorm")){
    if(mod != "GRMnorm"){
      if(mod == "2PLnorm"| mod == "HO2PLnorm-F"){
        x <- c("P_beta","P_alpha","F_beta","F_alpha")
      }else{
        x <- c("P_beta","P_alpha","P_nu", "F_beta","F_alpha")
      }
    }else{
      x <- c("P_b", "P_alpha","F_b", "F_alpha")
    }
  }
  if(mod %in% c("1PLnorm", "Testlet1PLnorm","HO1PLnorm-C", "HO1PLnorm-P", "HO1PLnorm-F")){
    x <- c("P_beta", "F_beta")
    
    if(mod != "1PLnorm" & mod != "HO1PLnorm-F"){
      x <- c(x, "P_nu")
    }
  }
  if(mod == "3PLnorm"){
    x <- c("P_beta", "P_alpha", "P_gamma","F_beta", "F_alpha", "F_gamma")
  }
  combos <- lapply(seq_along(x), function(k) combn(x, k, simplify = FALSE))
  prior_fixed_list <- unlist(combos, recursive = FALSE)
  prior_fixed_list_clean <- list()
  prior_fixed_list_clean[[1]] <- " "
  u <- 2
  for(i in seq_len(length(prior_fixed_list))){
    if(!any(duplicated(substr(prior_fixed_list[[i]], 2,10)))){
      prior_fixed_list_clean[[u]] <- prior_fixed_list[[i]]
      u <- u + 1
    }
  }
  all_combi[[mod]] <- prior_fixed_list_clean
}


# Simulate data for pipeline tests 
# Simulate data Testletnorm ---
sim_response <- function(delta, nperson = 1000){
  D <- 3
  itemD <- rep(1:3, each = 30)
  alpha <- NULL; beta <- NULL
  for(d in 1:D){
    alpha <- c(alpha, rlnorm(30, 0,0.25)) #/mean(delta_age[,d]))
    beta <- c(beta, runif(30, quantile(delta[,d], pnorm(-3)),quantile(delta[,d], pnorm(3))))
  }
  nitem <- 90
  psolve <- matrix(0, nperson, nitem)
  for(j in 1:nitem){
    for (i in 1:nperson){
      psolve[i, j] <- exp(alpha[j] * (delta[i,itemD[j]] - beta[j]))/
        (1 + exp(alpha[j] * (delta[i,itemD[j]]  - beta[j])))
    }
  }
  R <- (matrix(runif(nitem * nperson), nperson, nitem) < psolve) * 1
  
  sim_score <- data.frame(score = c(rowSums(R[,1:30]),rowSums(R[,31:60]),rowSums(R[,61:90])),
                          domain = rep(1:3, each = nperson))
  return(list(sim_score,R))
}
sim_setting <- function(example, seed = 1234, nperson){
  set.seed(seed)
  age <- runif(nperson, -2,2)
  # with intercept
  age_poly <-  model.matrix( ~ poly(age, 3, raw = T), data = data.frame(age = age))
  D <- 3
  itemD <- rep(1:3, each = 30)
  nperson <- nperson
  normed_theta <- rnorm(nperson)
  normed_unique <- matrix(rnorm(nperson*D), ncol = 3)
  if(example != "Testletnorm"){
    mean_age <- cbind(age_poly[,1:4] %*%  c(0,0.25, -0.1, 0.01),
                      age_poly[,1:4] %*%  c(-0.4,0.3,0,0),
                      age_poly[,1:4] %*%  c(-0.8,0.25,0.12, 0.02))
    sd_age <- cbind(exp(age_poly[,1:3] %*% c(-0.3, -0.08, 0.0)),
                    exp(age_poly[,1:3] %*% c(-0.3 ,0,0)),
                    exp(age_poly[,1:3] %*% c(-0.2, 0.07,  -0.02)))
    xi <- c(0.9,0.5,0.75)
  }
  if(example == "HOnorm_C"){
    correlation <- rep(xi, each = length(age))
    
    loading_theta <- correlation * sd_age
    loading_epsilon <- sqrt(1 - correlation^2) * sd_age
  }
  if(example == "HOnorm_F"){
    weight <- cbind((exp(age_poly[,1:3] %*% c(-0.2, -0.06, -0.015))),
                    (exp(age_poly[,1:3] %*% c(-0.6, 0.1, -0.015))),
                    exp(age_poly[,1:3] %*% c(-0.5, 0.05, 0.03)))
    sd_age <- sd_age
    theta_age <- sd_age * weight
    epsilon_age <- sd_age *sqrt(1-weight^2)
    
    correlation <- theta_age/sqrt(theta_age^2 + epsilon_age^2)
    
    loading_theta <- theta_age
    loading_epsilon <- epsilon_age
  }
  if(example == "HOnorm_P"){
    lambda_age <- exp(age_poly[,2:3] %*% c(0.06, -0.05))
    correlation <- cbind(xi[1]* lambda_age,
                         xi[2]*lambda_age,
                         xi[3]*lambda_age)
    
    loading_theta <- correlation * sd_age
    loading_epsilon <- sqrt(1 - correlation^2) * sd_age
    
  }
  if(example == "Testletnorm"){
    mean_age <- matrix(rep(age_poly[,1:4] %*%  c(0,0.25, -0.1, 0.01),3), ncol = 3)
    sd_age <- matrix(rep(exp(age_poly[,1:3] %*% c(-0.3, -0.08, 0.0)),3), ncol = 3)
    correlation  <-  rep(c(0.9,0.5,0.75), each = nperson)
    
    loading_theta <- correlation * sd_age
    loading_epsilon <- sqrt(1 - correlation^2) * sd_age
  }
  delta <- normed_theta*loading_theta + normed_unique*loading_epsilon + mean_age
  normed_delta <- (delta - mean_age)/sd_age
  
  info <- data.frame(age = rep(age, 15),
                     value = c(as.vector(sd_age), as.vector(mean_age),
                               as.vector(correlation), as.vector(loading_theta),
                               as.vector(loading_epsilon)),
                     moment = rep(c("sd", "mu","cor", "load_theta", "load_epsilon"), each = nperson *3),
                     domain = rep(rep(c("A","B","C"), each =nperson),5),
                     example = example)
  
  info2 <- data.frame(values = c(normed_theta,
                                 as.vector(normed_delta),
                                 as.vector(delta)),
                      age = rep(age, 7),
                      parameter = c(rep("theta_n", length(age)), rep(c("delta_n", "delta"), each = length(age)*3)),
                      domain = c(rep("general", nperson), rep(rep(LETTERS[1:3], each = nperson),2)))
  
  sim_score <- sim_response(delta, nperson = nperson)
  
  R <- sim_score[[2]]
  sim_score <- sim_score[[1]]
  sim_score$age <- rep(age, 3)
  infos <- list(response = data.frame(R, age = age),
                itemD = itemD,
                true = data.frame(age = rep(age,7),
                                  true_value = c(delta, normed_delta, normed_theta),
                                  trait = c(rep("delta", nperson*3),
                                            rep("n_delta", nperson*3),
                                            rep("n_theta", nperson)),
                                  domain = c(rep(LETTERS[1:3], each = nperson),
                                             rep(LETTERS[1:4], each = nperson)),
                                  n_trait = c(rep(paste0("delta",1:3), each = nperson),
                                              rep(paste0("n_delta",1:3), each = nperson),
                                              rep("n_theta", nperson))))
  return(infos)
}

set.seed(1234)
data_T <- sim_setting(example = "Testletnorm", seed = 1234, nperson = 1000)

# Simulate data 1PLnorm ----
Sim_data <- function (nperson = 1000, nitem = 50, seed = 1234, unidim = T,
                      Option, alpha, beta,
                      formular_mean =  ~ poly(age,3, raw = TRUE),
                      formular_sd = ~ poly(age,3, raw = TRUE),
                      int_mean = F, int_sd = F){
  
  set.seed(seed)
  age <- runif(nperson, -2, 2)
  normedability <- rnorm(nperson,0,1)
  age_poly <-  model.matrix( ~ poly(age, 3, raw = T), data = data.frame(age = age))
  if(Option == "A"){
    sd_age <- 0.9
    mean_age <-  age_poly[,1:2] %*% (c(0,0.4))
  }
  if(Option == "B"){
    mean_age <- age_poly[,1:4] %*% (c(0.2076,0.5906,-0.15, 0.01))
    sd_age <- age_poly[,1:3] %*% (c(0.35,0, 0.1)/0.65)
  }
  if(Option == "C"){
    sd_age <- age_poly[,1:2] %*%(c(0.42,0.05)/0.65)
    mean_age <- age_poly[,1:3] %*% (c(0.7,0.88,-0.55)*0.65)
  }
  if(int_mean == F){
    age_poly <-  model.matrix(formular_mean, data = data.frame(age = age))
    int_mean <- 1:ncol(age_poly)
  }
  if(int_sd == F){
    age_poly <-  model.matrix(formular_sd, data = data.frame(age = age))
    int_sd <- 1:ncol(age_poly)
  }
  ability <- normedability*sd_age  + mean_age
  psolve <- matrix(0, nperson, nitem)
  for (i in 1:nperson){
    for(j in 1:nitem){
      psolve[i, j] <- exp(alpha[j] * (ability[i] - beta[j]))/
        (1 + exp(alpha[j] * (ability[i]  - beta[j])))
    }
  }
  R <- (matrix(runif(nitem * nperson), nperson, nitem) < psolve) * 1
  response <- data.frame(R, age = age)
  return(response)
}
set.seed(1234)
data_1PL <- Sim_data(Option = "B", alpha = rep(1, 50), beta = runif(50,-3,3))


# Simulate data GRMnorm ----
sim_response_GRM <- function(seed, nperson = 1000, nitem =40){
  logistic <- function(x){1/(1 + exp(-x))}
  set.seed(seed)
  age <- runif(nperson, -2,2)
  age_poly <-  model.matrix( ~ poly(age, 3, raw = T), data = data.frame(age = age))
  normed_theta <- rnorm(nperson)
  mean_age <- age_poly[,1:4] %*%  c(-0.8,0.25,0.12, 0.02)
  sd_age <- exp(age_poly[,1:3] %*% c(-0.2, 0.07,  -0.02))
  
  theta <- (normed_theta * sd_age) + mean_age
  
  alpha <- rlnorm(nitem, 0,0.25)
  n_categories_per_item <- sample(3:5, size = nitem, replace = TRUE)
  
  b <- list()
  for(n in seq_len(nitem)){
    b[[n]] <- sort(runif(n_categories_per_item[n], -3.3,3))
  }
  R <- matrix(NA, nrow = nperson, ncol = nitem)
  for (i in seq_len(nitem)) {
    
    thresholds <- b[[i]]
    m <- length(thresholds) + 1  # number of categories for this item
    
    # P(X >= k) for k = 1..m-1, plus boundary conditions P(X >= 0) = 1, P(X >= m) = 0
    p_cumulative <- matrix(NA, nrow = nperson, ncol = m + 1)
    p_cumulative[, 1] <- 1                     # P(X >= category 1) = 1
    p_cumulative[, m + 1] <- 0                 # P(X >= category m+1) = 0
    
    for (k in seq_len(m - 1)) {
      p_cumulative[, k + 1] <- logistic(alpha[i] * (theta - thresholds[k]))
    }
    
    # category probabilities are successive differences
    p_category <- p_cumulative[, 1:m] - p_cumulative[, 2:(m + 1)]
    
    # sample a category for each person based on these probabilities
    R[, i] <- apply(p_category, 1, function(probs) {
      sample(seq_len(m), size = 1, prob = probs)
    })
  }
  response <- cbind(R, age = age)
  return(list(response = response,
              m_categories = n_categories_per_item  + 1))
}

data_GRM <- sim_response_GRM(seed = 1234)

rm(Sim_data, sim_response, sim_response_GRM, sim_setting)

