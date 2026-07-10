
set.seed(1234)
random_subset <- list()
for(i in irt_models){
  temp <- all_combi[[i]]  
  int <- sample( seq_len(length(temp)), 1)
  random_subset[[i]] <- temp[[int]]
}

test_that("IRT models can be fitted and information can be extracted", {
  skip_if_not_installed("cmdstanr")
  
  skip_on_cran()
  
for(irtmod in irt_models){
  temp <- random_subset[[irtmod]]
  if(irtmod %in% c("1PLnorm","2PLnorm","3PLnorm")){
    nitem <- 50
  }
  if(irtmod %in% c("GRMnorm")){
    nitem <- 40
  }
  if(irtmod %in% c("Testlet1PLnorm","Testlet2PLnorm",
                   "HO1PLnorm-C", "HO2PLnorm-C",
                   "HO1PLnorm-P", "HO2PLnorm-P",
                   "HO1PLnorm-F", "HO2PLnorm-F")){
    nitem <- 90
  }
  for(agemod in age_models){
      if(any(substr(temp,1, 1) == "P")){
        prior_names <- substr(temp[which(substr(temp,1, 1) == "P")],3,10) 
        prior <- list()
        for(n in prior_names){
          if(n != "nu" & n != "b"){
            temp2 <- list("A" = runif(nitem, 0.5, 2),
                          "B" = runif(nitem, 0.2, 0.7))
          }
          if(n == "nu"){
            temp2 <- list("A" = runif(3, 0.2,0.9),
                          "B" = runif(3, 0.1, 0.1))
          }
          if(n != "b"){
            temp2 <-  setNames(temp2, c(paste0(n, "_prior_mu"),
                                        paste0(n, "_prior_sd")))
          }
          if (n == "b"){
            b_prior_mu <- list()
            b_prior_sd <- list()
            for(nit in seq_len(nitem)){
              b_prior_mu[[nit]] <- sort(runif(max(data_GRM$m_categories) -1, 0.5, 2))
              b_prior_sd[[nit]] <- sort(runif(max(data_GRM$m_categories) -1, 0.2, 0.7))
              
            }
            temp2 <- list("b_prior_mu" =  b_prior_mu,
                          "b_prior_sd" = b_prior_sd)
          }
          prior[[n]] <- temp2
          
        }
      }else{prior <- NULL; prior_names <- NULL}
      if(any(substr(temp,1, 1) == "F")){
        fixed_names <- substr(temp[which(substr(temp,1, 1) == "F")],3,10) 
        fixed<- list()
        
        for(n in fixed_names){
          if(n != "nu" & n != "b"){
            fixed[[n]] <- runif(nitem, 0.5, 2)
          }
          if(n == "nu"){
            fixed[[n]] <-  runif(3, 0.2,0.9)
          }
          if(n == "b"){
            b <- list()
            for(nit in seq_len(nitem)){
              b[[nit]] <- sort(runif(max(data_GRM$m_categories) -1, 0.5, 2))
            }
            fixed[["b"]] <- b
            
          }
        }}else{fixed <- NULL; fixed_names <- NULL}
    
    print(paste(irtmod, agemod, prior_names, fixed_names))
      if(irtmod %in% c("1PLnorm","2PLnorm","3PLnorm")){
        data_temp <- data_prep(raw_data = data_1PL, age_variable = "age", item_variables = 1:50, irt_model = irtmod, age_model = agemod, 
                               poly_mean = 4, poly_sd = 3,
                               prior_knowlegde = prior, parameter_fixed = fixed)
      }
      if(irtmod %in% c("GRMnorm")){
        data_temp <- data_prep(raw_data = data_GRM$response, age_variable = "age", item_variables = 1:40,
                               irt_model = "GRMnorm", age_model = agemod, 
                               poly_mean = 4, poly_sd = 3, 
                               K_item =  data_GRM$m_categories,
                               prior_knowlegde = prior, parameter_fixed = fixed)
      }
      if(irtmod %in% c("Testlet1PLnorm","Testlet2PLnorm",
                       "HO1PLnorm-C", "HO2PLnorm-C",
                       "HO1PLnorm-P", "HO2PLnorm-P",
                       "HO1PLnorm-F", "HO2PLnorm-F")){
        data_temp <- data_prep(raw_data = data_T$response, age_variable = "age", item_variables = 1:90,
                               irt_model = irtmod, age_model = agemod, itemD = data_T$itemD,
                               poly_mean = 4, poly_sd = 3, 
                               prior_knowlegde = prior, parameter_fixed = fixed)
      }
      
     
      expect_no_error(fit <- fit_IRTnorm(data = data_temp, seed = NULL,
                         iter_warmup = 5, iter_sampling = 5,chains = 2))
      expect_no_error(extract_info (fit = fit, data = data_temp))
      expect_no_error(info_new_norming(fit = fit, data = data_temp))
      
    }
  }
})


