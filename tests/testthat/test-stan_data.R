
test_that("stan code compiles", {
  skip_on_cran()
  irt_models <- c("1PLnorm","2PLnorm","3PLnorm","GRMnorm",
                  "Testlet1PLnorm","Testlet2PLnorm",
                  "HO1PLnorm-C", "HO2PLnorm-C",
                  "HO1PLnorm-P", "HO2PLnorm-P",
                  "HO1PLnorm-F", "HO2PLnorm-F")
  
  age_models <- c("splines","poly")
    for (irtmod in irt_models) {
      temp <- all_combi[[irtmod]]
      for(agemod in age_models){
        for(set in seq_len(length(temp))){
          if(any(substr(temp[[set]],1, 1) == "P")){
            prior <- substr(temp[[set]][which(substr(temp[[set]],1, 1) == "P")],3,10) 
          }else{prior <- NULL}
          if(any(substr(temp[[set]],1, 1) == "F")){
            fixed <- substr(temp[[set]][which(substr(temp[[set]],1, 1) == "F")],3,10) 
          }else{fixed <- NULL}
          print(paste("Failed for",agemod, "_", irtmod, "_P_", prior, "_F_", fixed))
        expect_no_error(get_compiled_model(age_model = agemod,
                                    irt_model = irtmod,
                                    prior_knowlegde = prior,
                                    parameter_fixed  = fixed))
      }
      }
    }
  })
  