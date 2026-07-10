

#check if the data resulting form data_prep aligns with the stan model input 

split_stan_statements <- function(data_block_text) {
  # Remove comments
  text <- gsub("//[^\n]*", "", data_block_text)
  text <- gsub("/\\*.*?\\*/", "", text, perl = TRUE)
  
  # Strip outer data { ... }
  text <- sub("^\\s*data\\s*\\{", "", text)
  text <- sub("\\}\\s*$", "", text)
  
  # Split on semicolons that are NOT inside <...> or [...]
  # Walk character by character tracking bracket depth
  chars <- strsplit(text, "")[[1]]
  depth <- 0
  stmt <- character(0)
  buf <- c()
  
  for (ch in chars) {
    if (ch %in% c("<", "[")) depth <- depth + 1
    if (ch %in% c(">", "]")) depth <- depth - 1
    if (ch == ";" && depth == 0) {
      stmt <- c(stmt, paste(buf, collapse = ""))
      buf <- c()
    } else {
      buf <- c(buf, ch)
    }
  }
  stmt <- trimws(stmt)
  stmt[stmt != ""]
}

parse_stan_decl <- function(stmt) {
  # Remove bracketed sections (<...> and [...]) to isolate "type-ish prefix + names"
  stripped <- gsub("<[^>]*>", "", stmt)
  stripped <- gsub("\\[[^\\]]*\\]", "", stripped)
  stripped <- trimws(stripped)
  
  tokens <- strsplit(stripped, "\\s+")[[1]]
  if (length(tokens) < 2) return(list())
  
  # Names are the last whitespace-token-group (may have commas)
  varnames_part <- tokens[length(tokens)]
  varnames <- trimws(strsplit(varnames_part, ",")[[1]])
  
  # base_type is whatever's left after removing brackets and the varnames part
  type_tokens <- tokens[1:(length(tokens) - 1)]
  base_type <- paste(type_tokens, collapse = " ")
  
  list(
    names = varnames,
    base_type = base_type
  )
}

parse_stan_data_block <- function(data_block_text) {
  stmts <- split_stan_statements(data_block_text)
  results <- list()
  for (s in stmts) {
    parsed <- parse_stan_decl(s)
    for (vn in parsed$names) {
      results[[vn]] <- list(
        base_type = parsed$base_type
      )
    }
  }
  results
}

check_stan_data <- function(data_block_text, data_list) {
  expected <- parse_stan_data_block(data_block_text)
  
  expected_names <- names(expected)
  actual_names <- names(data_list)
  
  missing <- setdiff(expected_names, actual_names)

  good <- NULL
  fail <- NULL
  variabs <- intersect(expected_names, actual_names)
  for(i in variabs){
    if(expected[[i]]$base_type == "int" |expected[[i]]$base_type == "array int"){
      if(is.integer(data_list[[i]])| is.numeric(data_list[[i]])){
      good <- c(good, i)
    }else{
      fail <- c(fail,i)
    }
  }
    if(substr(expected[[i]]$base_type,1,6)== "matrix"){
      if(is.matrix(data_list[[i]])){
        good <- c(good, i)
      }else{
        fail <- c(fail,i)
      }
    }
    
    
  }
   
  list(
    missing_in_list = missing,
    good = good, 
    fail= fail)
} 

get_data_prep <- function(irt_model, age_model, prior_knowlegde = NULL, parameter_fixed = NULL){
  is_2pl  <- irt_model %in% c("2PLnorm", "3PLnorm", "Testlet2PLnorm","HO2PLnorm-C",
                              "HO2PLnorm-P","HO2PLnorm-F")
  is_multi <- irt_model %in% c("Testlet2PLnorm", "Testlet1PLnorm",
                               "HO1PLnorm-C","HO2PLnorm-C",
                               "HO1PLnorm-P","HO2PLnorm-P",
                               "HO2PLnorm-F", "HO1PLnorm-F")
  is_flex_multi <- irt_model %in% c("HO1PLnorm-P","HO2PLnorm-P",
                                    "HO2PLnorm-F", "HO1PLnorm-F")
  
  p2 <- create_block_data(irt_model = irt_model, age_model = age_model,
                          is_multi = is_multi, is_2pl = is_2pl, is_flex_multi = is_flex_multi,
                          prior_knowlegde = prior_knowlegde, parameter_fixed  =  parameter_fixed )
  return(p2)
}


test_that("data_prep transforms data in correct stan format", {
  skip_on_cran()
for(irtmod in irt_models){
  temp <- all_combi[[irtmod]]
  if(irtmod %in% c("1PLnorm","2PLnorm","3PLnorm")){
    nitem <- 50
    rawdata <-  data_1PL
  }
  if(irtmod %in% c("GRMnorm")){
    nitem <- 40
    rawdata <- data_GRM$response
  }
  if(irtmod %in% c("Testlet1PLnorm","Testlet2PLnorm",
                   "HO1PLnorm-C", "HO2PLnorm-C",
                   "HO1PLnorm-P", "HO2PLnorm-P",
                   "HO1PLnorm-F", "HO2PLnorm-F")){
    nitem <- 90
    rawdata <- data_T$response
  }
  for(agemod in age_models){
    for(set in seq_len(length(temp))){
      if(any(substr(temp[[set]],1, 1) == "P")){
        prior_names <- substr(temp[[set]][which(substr(temp[[set]],1, 1) == "P")],3,10) 
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
          
        }}else{prior <- NULL; prior_names <- NULL}
      if(any(substr(temp[[set]],1, 1) == "F")){
        fixed_names <- substr(temp[[set]][which(substr(temp[[set]],1, 1) == "F")],3,10) 
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
     
      data_temp <- data_prep(raw_data = rawdata, age_variable = "age", item_variables = seq_len(nitem),
                               irt_model = irtmod, age_model = agemod, 
                               poly_mean = 4, poly_sd = 3, 
                               K_item =  data_GRM$m_categories, itemD = data_T$itemD,
                               prior_knowlegde = prior, parameter_fixed = fixed)
      
      stan_code <- get_data_prep(irt_model = irtmod, age_model = agemod, 
                                 prior_knowlegde = prior_names, parameter_fixed = fixed_names)
      test <- check_stan_data(data_block_text = stan_code, data_list = data_temp)
      expect_equal(test$missing_in_list,  character(), 
                   info = paste("For the combination", irtmod, "_", agemod, "_P", 
                               paste(prior_names, collapse = "_"), "_F", 
                               paste(fixed_names, collapse = "_"), 
                               "missing input for", paste(test$missing_in_list, collapse = "_")))
      expect_null(test$fail, 
                  info = paste("For the combination", irtmod, "_", agemod, "_P", 
                               paste(prior_names, collapse = "_"), "_F", 
                               paste(fixed_names, collapse = "_"), "_imput_",
                                     paste(test$missing_in_list, collapse = "_"), 
                               "is not in the correct format"))
      
    }
  }
}
})

