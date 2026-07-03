create_block_data <- function(irt_model, age_model, is_2pl, is_multi, is_flex_multi, prior_knowledge = NULL, parameter_fixed = NULL){
  # data
  general_data <- "data {int<lower=1> I;
  int<lower=1> J;
  int<lower=1> N;
  array[N] int<lower=1, upper=I> ii;
  array[N] int<lower=1, upper=J> jj;"

  if(irt_model == "GRMnorm"){
    general_data <- paste(general_data, "int<lower=2> K_max;
    array[I] int<lower=2, upper=K_max> K_item;
    array[N] int<lower=0, upper=K_max> y;")
  }else{
    general_data <- paste(general_data, "array[N] int<lower=0, upper=1> y;")
  }
  if(is_multi){
    general_data <- paste(general_data,"int<lower=1> U;
int<lower = 1> D;
array[N] int<lower=1, upper=D> dd;
array[N] int JD_loadings;")
  }


  if(!is.null(prior_knowledge)){
    # I think ,"lambda","unique_sigma" does not make sense
    dim_tab <- data.frame(param = c("alpha","beta", "gamma", "nu"),
                          dim = c("I","I","I","D"))
    for(j in prior_knowledge){
      if(j == "b"){
        general_data <- paste(general_data, paste0("array[I] ordered[K_max - 1] ", j,"_prior_mu;
                                                   array[I] ordered[K_max - 1] ", j,"_prior_sd;"))
      }
      else if (j %in% dim_tab$param) {
        d <- dim_tab$dim[dim_tab$param == j]
        general_data <- paste(general_data,
                              paste0("vector[",d,"] ", j, "_prior_mu;
                                      vector[",d,"] ", j, "_prior_sd;"))
      }
    }
  }

    if(is_flex_multi){
      age_data_chunk <- switch(age_model,
                               poly = "
      int S;
      int W;
      int M;
      matrix[J, S] age;
      matrix[J, W] age2;
      matrix[J, M] age3;
    ",
                               splines = "
      int S_fixed;
      int K1;
      int W_fixed;
      int K2;
      int M_fixed;
      int K3;
      matrix[J, S_fixed] X1;
      matrix[J, K1]      Z1;
      matrix[J, W_fixed] X2;
      matrix[J, K2]      Z2;
      matrix[J, M_fixed] X3;
      matrix[J, K3]      Z3;
    "
      )
    }else{
      age_data_chunk <- switch(age_model,
                               poly = "
      int S;
      int W;
      matrix[J, S] age;
      matrix[J, W] age2;
    ",
                               splines = "
      int S_fixed;
      int K1;
      int W_fixed;
      int K2;
      matrix[J, S_fixed] X1;
      matrix[J, K1]      Z1;
      matrix[J, W_fixed] X2;
      matrix[J, K2]      Z2;
    "
      )
    }


    general_data <- paste(general_data, age_data_chunk)

    if(!is.null(parameter_fixed)){
      if(irt_model != "GRMnorm" && "beta" %in% parameter_fixed){
        general_data <- paste(general_data,"vector[I] beta;")

      }
      if(irt_model == "GRMnorm" && "b" %in% parameter_fixed){
          general_data <- paste(general_data, "array[I] ordered[K_max - 1] b;")
      }

      if((is_2pl | irt_model == "GRMnorm") && "alpha" %in% parameter_fixed){
        general_data <- paste(general_data,"vector<lower=1e-6>[I] alpha;")
      }

      if(irt_model == "3PLnorm" && "gamma" %in% parameter_fixed){
        general_data <- paste(general_data, "vector<lower=0, upper = 1>[I] gamma;")
      }
    }

  return(general_data)
}

create_block_param <- function(irt_model, age_model = "poly", is_2pl, is_multi, is_flex_multi,parameter_fixed  = NULL){
HO_F <- (irt_model ==  "HO1PLnorm-F"|  irt_model  ==  "HO2PLnorm-F")
HO_P <- (irt_model ==  "HO1PLnorm-P"|  irt_model  ==  "HO2PLnorm-P")
HO_C <- (irt_model == "HO1PLnorm-C" | irt_model  == "HO2PLnorm-C")
uni_age <- (!HO_P & !HO_C & !HO_F)
  block_param <- "}
  parameters {
  vector[J] z_theta;"

  if(age_model == "poly"){
    if(uni_age){
      block_param <- paste(block_param,
                           "vector[S] p1_theta;
     vector[W] p2_theta;")
    }else{
    block_param <- paste(block_param,
                         "array[D] vector[S] p1_theta;
                          array[D] vector[W] p2_theta;")
  }
  if(HO_P){
    block_param <- paste(block_param,
                         "vector[M] p3_theta;")
  }
  if(HO_F){
    block_param <- paste(block_param,
                         "array[D] vector[M] p3_theta;")
  }
  }


  if(age_model == "splines"){
    if(uni_age){
      block_param <- paste(block_param,
                         "vector[S_fixed] l1;
    vector[K1]      u1;
    vector[W_fixed] l2;
    vector[K2]      u2;
    real<lower=1e-4>   tau1;
    real<lower=1e-4>   tau2;")
    }else{
      block_param <- paste(block_param,
                         "array[D] vector[S_fixed] l1;
    array[D] vector[K1]      u1;
    array[D] vector[W_fixed] l2;
    array[D] vector[K2]      u2;
    vector<lower=1e-4>[D]    tau1;
    vector<lower=1e-4>[D]   tau2;")
    }
    if(HO_P){
      block_param <- paste(block_param,
                           "vector[M_fixed] l3;
                            vector[K3]      u3;
                            real<lower=1e-4>   tau3;")
    }
    if(HO_F){
      block_param <- paste(block_param,
                           "array[D] vector[M_fixed] l3;
                            array[D] vector[K3]      u3;
                            vector<lower=1e-4>[D]   tau3;")
    }
  }

  if(irt_model != "GRMnorm" & !any(parameter_fixed == "beta")){
    block_param <- paste(block_param, "vector[I] beta;")

  }
  if(irt_model == "GRMnorm"){
    if( !any(parameter_fixed == "b")){
      block_param <- paste(block_param, "array[I] ordered[K_max - 1] b;")
    }
    if( !any(parameter_fixed == "alpha")){
    block_param <- paste(block_param, "vector<lower=1e-6>[I] alpha;")
    }
  }

  if(is_2pl & !any(parameter_fixed == "alpha")){
    block_param <- paste(block_param, "vector<lower=1e-6>[I] alpha;")
  }

  if(irt_model == "3PLnorm" & !any(parameter_fixed == "gamma")){
    block_param <- paste(block_param, "vector<lower=0, upper = 0.9>[I] gamma;")
  }

if(is_multi  && !HO_F){

  block_param <- paste(block_param, "vector<lower=0, upper=1>[D] nu;")
}
  if(is_multi){
    block_param <- paste(block_param, "vector[U] unique;")
  }
  return(block_param)

}

create_block_transparam <- function(irt_model, age_model = "poly", is_2pl, is_multi, is_flex_multi){

  if(!is_multi |irt_model == "Testlet2PLnorm" |irt_model == "Testlet1PLnorm"){
    block_transparam <- "}
  transformed parameters{
  vector[J] theta_mean;
  vector[J] theta_sigma;"

    if(!is_multi){
      block_transparam <- paste(block_transparam,
                                " vector[J] theta;")
    }

    if(is_multi){
      block_transparam <- paste(block_transparam, "vector[N] substraits;")
    }
    if(age_model == "poly"){
      block_transparam <- paste(block_transparam,
                                "theta_mean = age * p1_theta;
                       theta_sigma = exp(age2 * p2_theta);")
    }
    if(age_model == "splines"){
      block_transparam <- paste(block_transparam,
                                "theta_mean  =  X1 * l1 + Z1 * u1;
                       theta_sigma = exp(X2 * l2 + Z2 * u2);")
    }

    if(!is_multi){
      block_transparam <- paste(block_transparam,
                                "theta =  theta_mean + theta_sigma .* z_theta;")
    }
    if(irt_model != "GRMnorm"){
      block_transparam <- paste(block_transparam, "")
    }

    if(is_multi){
      block_transparam <- paste(block_transparam,
                                "for (n in 1:N)
  substraits[n] = theta_mean[jj[n]] + z_theta[jj[n]] *nu[dd[n]] * theta_sigma[jj[n]]  + unique[JD_loadings[n]]*sqrt(1 - nu[dd[n]]^2 )*theta_sigma[jj[n]];")
    }

  }


  if(is_flex_multi |irt_model ==  "HO1PLnorm-C" |irt_model ==  "HO2PLnorm-C"){

    block_transparam <-
      "}
      transformed parameters {
      vector[U] theta_mean;
      vector<lower=0>[U] theta_sigma;
      vector[N] substraits;"


    if(irt_model ==  "HO1PLnorm-P" |irt_model ==  "HO2PLnorm-P"){
      block_transparam <- paste(block_transparam, "vector<lower=0, upper=1>[J] lambda;")
    }
    if(irt_model == "HO1PLnorm-F" | irt_model == "HO2PLnorm-F"){
      block_transparam <- paste(block_transparam,"vector[U] unique_sigma;")
    }

    block_transparam <- paste(block_transparam,
                              "for (d in 1:D) {
        int start_idx = (d - 1) * J + 1;
        int end_idx = d * J;")


    if(age_model == "poly"){
      if(irt_model == "HO1PLnorm-F" | irt_model == "HO2PLnorm-F"){
        block_transparam <- paste(block_transparam,
                                  "unique_sigma[start_idx:end_idx] = log1p_exp(age3 * p3_theta[d]);")
      }
      block_transparam <- paste(block_transparam,
                                "theta_mean[start_idx:end_idx] = age * p1_theta[d];
                                theta_sigma[start_idx:end_idx] = log1p_exp(age2 * p2_theta[d]);}")

      if(irt_model ==  "HO1PLnorm-P" |irt_model ==  "HO2PLnorm-P"){
        block_transparam <- paste(block_transparam, "lambda = inv_logit(age3 * p3_theta);")
      }
    }
    if(age_model == "splines"){
      if(irt_model == "HO1PLnorm-F" | irt_model == "HO2PLnorm-F"){
        block_transparam <- paste(block_transparam,
                                  "unique_sigma[start_idx:end_idx] = exp(X3 * l3[d] + Z3 * u3[d]);")
      }
      block_transparam <- paste(block_transparam,
                                "theta_mean[start_idx:end_idx] = X1 * l1[d] + Z1 * u1[d];
                                theta_sigma[start_idx:end_idx] = exp(X2 * l2[d] + Z2 * u2[d]);}")

      if(irt_model ==  "HO1PLnorm-P" |irt_model ==  "HO2PLnorm-P"){
        block_transparam <- paste(block_transparam, "lambda = inv_logit(X3 * l3 + Z3 * u3);")
      }

    }



    if(irt_model ==  "HO1PLnorm-C" |irt_model ==  "HO2PLnorm-C"){
      block_transparam <- paste(block_transparam,
                                "substraits = theta_mean[JD_loadings] +
          z_theta[jj] .* nu[dd] .* theta_sigma[JD_loadings] +
          unique[JD_loadings] .* sqrt(1 - nu[dd]^2) .* theta_sigma[JD_loadings];")

    }
    if(irt_model ==  "HO1PLnorm-P" |irt_model ==  "HO2PLnorm-P"){
      block_transparam <- paste(block_transparam, "substraits = theta_mean[JD_loadings] +
                    z_theta[jj] .* nu[dd] .*lambda[jj].* theta_sigma[JD_loadings] +
                    unique[JD_loadings] .* sqrt(1 - (nu[dd].*lambda[jj])^2) .* theta_sigma[JD_loadings];")
    }
    if(irt_model ==  "HO1PLnorm-F" |irt_model ==  "HO2PLnorm-F"){
      block_transparam <- paste(block_transparam, "substraits = theta_mean[JD_loadings] + z_theta[jj] .* (theta_sigma[JD_loadings]) + unique[JD_loadings].*(unique_sigma[JD_loadings]);")
    }
  }

  return(block_transparam)

}

create_block_model <- function(irt_model, age_model = "poly", is_2pl, is_multi, is_flex_multi, prior_knowledge = NULL,parameter_fixed  = FALSE){

  if(!is_flex_multi & irt_model !=  "HO1PLnorm-C" &irt_model !=  "HO2PLnorm-C"){
    if(age_model == "poly"){
      block_model <- "}model {
    p1_theta ~ normal(0,3);
    p2_theta ~ normal(0,1);"
    }
    if(age_model == "splines"){
      block_model <- "}model {
    l1 ~ normal(0, 5);
    l2 ~ normal(0, 5);
    tau1 ~ inv_gamma(2, 1);
    tau2 ~ inv_gamma(2, 1);
    u1 ~ normal(0, tau1);
    u2 ~ normal(0, tau2);"
    }
  }
  if( irt_model ==  "HO1PLnorm-C" |irt_model ==  "HO2PLnorm-C"| irt_model ==  "HO1PLnorm-P" |irt_model ==  "HO2PLnorm-P"){
    if(age_model == "poly"){
      block_model <- "}model {
      for (d in 1:D) {
    p1_theta[d] ~ normal(0, 3);
    p2_theta[d] ~ normal(0, 1);}"
    }
    if(age_model == "splines"){
      block_model <- "} model {
    tau1 ~ inv_gamma(2, 1);
    tau2 ~ inv_gamma(2, 1);
     for (d in 1:D) {
    l1[d] ~ normal(0, 10);
    l2[d] ~ normal(0, 10);
    u1[d] ~ normal(0, tau1[d]);
    u2[d] ~ normal(0, tau2[d]);}"
    }
  }
  if( irt_model ==  "HO1PLnorm-P" |irt_model ==  "HO2PLnorm-P"){
    if(age_model == "poly"){
      block_model <- paste(block_model, "p3_theta ~ normal(0, 1);")
    }
    if(age_model == "splines"){
      block_model <- paste(block_model,
                           "l3 ~ normal(0, 5);
                           tau3 ~ inv_gamma(2, 1);
                           u3 ~ normal(0, tau3);")
    }
    }

  if( irt_model ==  "HO1PLnorm-F" |irt_model ==  "HO2PLnorm-F"){
    if(age_model == "poly"){
      block_model <- "}model {
      for (d in 1:D) {
    p1_theta[d] ~ normal(0, 3);
    p2_theta[d] ~ normal(0, 1);
    p3_theta[d] ~ normal(0, 1);}"
    }
    if(age_model == "splines"){
      block_model <- "} model {
    tau1 ~ inv_gamma(2, 1);
    tau2 ~ inv_gamma(2, 1);
    tau3 ~ inv_gamma(2, 1);
     for (d in 1:D) {
    l1[d] ~ normal(0, 10);
    l2[d] ~ normal(0, 10);
    l3[d] ~ normal(0, 10);
    u1[d] ~ normal(0, tau1[d]);
    u2[d] ~ normal(0, tau1[d]);
    u3[d] ~ normal(0, tau2[d]);}"
    }
  }

  block_model <- paste(block_model,"z_theta ~ std_normal();")

  if(is.null(prior_knowledge) & is.null(parameter_fixed)){
    if(irt_model != "GRMnorm"){
      block_model <- paste(block_model, "beta ~ normal(0, 2);")
    }
    if(is_2pl){
      block_model <- paste(block_model,"alpha ~ lognormal(0, 0.5);")
    }
    if(irt_model == "3PLnorm"){
      block_model <- paste(block_model,"gamma ~ beta(1,10);
                             for (n in 1:N) {
      y[n] ~ logit_3pl(theta[jj[n]], alpha[ii[n]], beta[ii[n]], gamma[ii[n]]);
    }")
    }
  }

  if(!is.null(prior_knowledge)){
    if(irt_model != "GRMnorm"){
      if(any(prior_knowledge == "beta")){
        block_model <- paste(block_model, "beta ~ normal(beta_prior_mu, beta_prior_sd);")
      }
      else{
        if(!any(parameter_fixed == "beta")){
          block_model <- paste(block_model, "beta ~ normal(0, 2);")
        }
      }
    }
    if(is_2pl){
      if(any(prior_knowledge == "alpha")){
        block_model <- paste(block_model,"alpha ~ normal(alpha_prior_mu, alpha_prior_sd);")

      }else{
        if(!any(parameter_fixed == "alpha")){
        block_model <- paste(block_model,"alpha ~ lognormal(0, 0.5);")
}
      }
    }
    if(irt_model == "3PLnorm"){
      if(!any(parameter_fixed == "gamma")){
        if(any(prior_knowledge == "gamma")){
          block_model <- paste(block_model,"gamma ~ normal(gamma_prior_mu, gamma_prior_sd);")
        }else{
          block_model <- paste(block_model,"gamma ~ beta(1,10);")
        }
      }
      block_model <- paste(block_model,"for (n in 1:N) {
      y[n] ~ logit_3pl(theta[jj[n]], alpha[ii[n]], beta[ii[n]], gamma[ii[n]]);
                             }")
    }
  }




  if(irt_model == "2PLnorm"){
    block_model <- paste(block_model," y   ~ bernoulli_logit(alpha[ii] .* (theta[jj] - beta[ii]));")
  }
  if(irt_model == "1PLnorm"){
    block_model <- paste(block_model," y   ~ bernoulli_logit((theta[jj] - beta[ii]));")
  }
  if(irt_model == "GRMnorm"){
    if(!any(parameter_fixed == "alpha")){
      if(any(names(prior_knowledge) == "alpha")){
        block_model <- paste(block_model,"alpha ~ normal(alpha_prior_mu, alpha_prior_sd);
                             for (i in 1:I)")
      }else{
        block_model <- paste(block_model,
                             "alpha ~ lognormal(0, 0.5);")
      }
    }
   if(!any(parameter_fixed == "b")){
     block_model <- paste(block_model, "for (i in 1:I)")
      if(any(names(prior_knowledge) == "b")){
        block_model <- paste(block_model,"b[i] ~ normal(b_prior_mu[i], b_prior_sd[i]);")
      }else{
        block_model <- paste(block_model,"b[i] ~ normal(0, 1);")
      }
   }
    block_model <- paste(block_model,
                         "for (n in 1:N) {
      int Ki = K_item[ii[n]];

        if (y[n] == 1) {
          // P(Y = 1)
          target += log1m_inv_logit(
            alpha[ii[n]] * (theta[jj[n]] - b[ii[n]][1])
          );

        } else if (y[n] == Ki) {
          // P(Y = K_i)
          target += log_inv_logit(
            alpha[ii[n]] * (theta[jj[n]] - b[ii[n]][Ki - 1])
          );

        } else {
          // P(Y = k), k = 2,...,K_i-1
          target += log(
            inv_logit(alpha[ii[n]] * (theta[jj[n]] - b[ii[n]][y[n] - 1])) -
            inv_logit(alpha[ii[n]] * (theta[jj[n]] - b[ii[n]][y[n]])));
        }
      }")


  }
  if(is_multi){
    block_model <- paste(block_model,"unique ~ normal(0,1);")
    if(irt_model !=  "HO1PLnorm-F" & irt_model !=  "HO2PLnorm-F"){
      if(any(names(prior_knowledge) == "nu")){
        block_model <- paste(block_model,"nu ~ rnorm(nu_prior_mu,nu_prior_sd);")
      }else{
        block_model <- paste(block_model,"nu ~ beta(2,2);")
      }
    }

    if(!is_2pl){
      block_model <- paste(block_model,"y ~ bernoulli_logit((substraits - beta[ii]));")
    }

    if(is_2pl){
      block_model <- paste(block_model,"y ~ bernoulli_logit(alpha[ii].*(substraits - beta[ii]));")

    }
  }
  return(block_model)
}

create_block_generatequantities <- function(irt_model, is_2pl, is_multi, is_flex_multi){

  block_quantities <-
    "} generated quantities {
    array[N] int y_rep;
    vector[N] log_lik;"

  if(!is_multi){


    if(irt_model == "1PLnorm" |irt_model == "2PLnorm"|irt_model == "3PLnorm"){
      block_quantities <- paste(block_quantities ,"array[J] real theta_rep;
    theta_rep = normal_rng(theta_mean, theta_sigma);")
    }


    if(irt_model == "1PLnorm"){
      block_quantities <- paste(block_quantities ,
                                "y_rep = bernoulli_rng(inv_logit((to_vector(theta_rep[jj]) - beta[ii])));
  for(n in 1:N)
    log_lik[n] = bernoulli_logit_lpmf(y[n] | (theta[jj[n]] - beta[ii[n]]));
  }")
    }
    if(irt_model == "2PLnorm"){
      block_quantities <- paste(block_quantities ,
                                "y_rep = bernoulli_rng(inv_logit(alpha[ii].*(to_vector(theta_rep[jj]) - beta[ii])));
    for(n in 1:N){
      log_lik[n] = bernoulli_logit_lpmf(y[n] | alpha[ii[n]] * (theta[jj[n]] - beta[ii[n]]));
    }}")
    }
    if(irt_model == "3PLnorm"){
      block_quantities <- paste(block_quantities ,
                                "for(n in 1:N){
    real prob = logit_3pl_prob(theta_rep[jj[n]], alpha[ii[n]], beta[ii[n]], gamma[ii[n]]);
    y_rep[n] = bernoulli_rng(prob);
    log_lik[n] = logit_3pl_lpmf(y[n] | theta[jj[n]], alpha[ii[n]], beta[ii[n]], gamma[ii[n]]);
  }}")
    }

    if(irt_model == "GRMnorm"){
      block_quantities <- paste(block_quantities ,
                                "for (n in 1:N) {
    int Ki = K_item[ii[n]];
    vector[Ki] p;
    vector[Ki] eta;

      for (k in 1:(Ki - 1))
        eta[k] = inv_logit(alpha[ii[n]] * (theta[jj[n]] - b[ii[n]][k]));

      p[1] = 1 - eta[1];
       if (Ki > 2){
               for (k in 2:(Ki - 1))
        p[k] = eta[k - 1] - eta[k];
       }
      p[Ki] = eta[Ki - 1];

    log_lik[n] = categorical_lpmf(y[n] | p);

    y_rep[n] = categorical_rng(p);
  }}" )
    }
  }

  if(is_multi){
    block_quantities <- paste(block_quantities , "array[J] real theta_rep;
    array[U] real unique_rep;
    for (j in 1:J)
      theta_rep[j] = normal_rng(0,1);
    for (u in 1:U)
      unique_rep[u] = normal_rng(0,1);")


    if(irt_model == "HO1PLnorm-F"){
      block_quantities <- paste(block_quantities ,
                                "for(n in 1:N){
            y_rep[n] = bernoulli_rng(inv_logit((theta_mean[JD_loadings[n]] +
            theta_rep[jj[n]]* (theta_sigma[JD_loadings[n]])  +
            unique_rep[JD_loadings[n]]*(unique_sigma[JD_loadings[n]]) - beta[ii[n]])));
      log_lik[n] = bernoulli_logit_lpmf(y[n] | (theta_mean[JD_loadings[n]] +
                    z_theta[jj[n]] * (theta_sigma[JD_loadings[n]]) +
                    unique[JD_loadings[n]] * (unique_sigma[JD_loadings[n]]) - beta[ii[n]]));
    }
}")}
    if(irt_model == "HO2PLnorm-F"){
      block_quantities <- paste(block_quantities ,
                                "for(n in 1:N){
            y_rep[n] = bernoulli_rng(inv_logit(alpha[ii[n]]*(theta_mean[JD_loadings[n]] +
            theta_rep[jj[n]]* (theta_sigma[JD_loadings[n]])  +
            unique_rep[JD_loadings[n]]*(unique_sigma[JD_loadings[n]]) - beta[ii[n]])));
      log_lik[n] = bernoulli_logit_lpmf(y[n] | alpha[ii[n]] * (theta_mean[JD_loadings[n]] +
                    z_theta[jj[n]] * (theta_sigma[JD_loadings[n]]) +
                    unique[JD_loadings[n]] * (unique_sigma[JD_loadings[n]]) - beta[ii[n]]));
    }
}")
    }

    if(irt_model == "HO1PLnorm-P"){
      block_quantities <- paste(block_quantities ,"
      for(n in 1:N){
        y_rep[n] = bernoulli_rng(inv_logit((theta_mean[JD_loadings[n]] +
                                                           theta_rep[jj[n]] *nu[dd[n]] * lambda[jj[n]]*theta_sigma[JD_loadings[n]]  +
                                                           unique_rep[JD_loadings[n]]*sqrt(1 - (nu[dd[n]]*lambda[jj[n]])^2 )*theta_sigma[JD_loadings[n]] - beta[ii[n]])));
        log_lik[n] = bernoulli_logit_lpmf(y[n] |(theta_mean[JD_loadings[n]] +
                                                                   z_theta[jj[n]] * nu[dd[n]] *lambda[jj[n]]* theta_sigma[JD_loadings[n]] +
                                                                   unique[JD_loadings[n]] * sqrt(1 - (nu[dd[n]]*lambda[jj[n]])^2) * theta_sigma[JD_loadings[n]] - beta[ii[n]]));
      }}")
    }
    if(irt_model == "HO2PLnorm-P"){
      block_quantities <- paste(block_quantities ,"
      for(n in 1:N){
        y_rep[n] = bernoulli_rng(inv_logit(alpha[ii[n]]*(theta_mean[JD_loadings[n]] +
                                                           theta_rep[jj[n]] *nu[dd[n]] * lambda[jj[n]]*theta_sigma[JD_loadings[n]]  +
                                                           unique_rep[JD_loadings[n]]*sqrt(1 - (nu[dd[n]]*lambda[jj[n]])^2 )*theta_sigma[JD_loadings[n]] - beta[ii[n]])));
        log_lik[n] = bernoulli_logit_lpmf(y[n] | alpha[ii[n]] * (theta_mean[JD_loadings[n]] +
                                                                   z_theta[jj[n]] * nu[dd[n]] *lambda[jj[n]]* theta_sigma[JD_loadings[n]] +
                                                                   unique[JD_loadings[n]] * sqrt(1 - (nu[dd[n]]*lambda[jj[n]])^2) * theta_sigma[JD_loadings[n]] - beta[ii[n]]));
      }}")
    }

    if(irt_model == "HO1PLnorm-C"){
      block_quantities <- paste(block_quantities ,"
        for(n in 1:N){
            y_rep[n] = bernoulli_rng(inv_logit((theta_mean[JD_loadings[n]] +
            theta_rep[jj[n]] *nu[dd[n]] * theta_sigma[JD_loadings[n]]  +
            unique_rep[JD_loadings[n]]*sqrt(1 - nu[dd[n]]^2 )*theta_sigma[JD_loadings[n]] - beta[ii[n]])));
      log_lik[n] = bernoulli_logit_lpmf(y[n] |(theta_mean[JD_loadings[n]] +
                    z_theta[jj[n]] * nu[dd[n]] * theta_sigma[JD_loadings[n]] +
                    unique[JD_loadings[n]] * sqrt(1 - nu[dd[n]]^2) * theta_sigma[JD_loadings[n]] - beta[ii[n]]));
    }
}")
    }
    if(irt_model == "HO2PLnorm-C"){
      block_quantities <- paste(block_quantities ,"
        for(n in 1:N){
            y_rep[n] = bernoulli_rng(inv_logit(alpha[ii[n]]*(theta_mean[JD_loadings[n]] +
            theta_rep[jj[n]] *nu[dd[n]] * theta_sigma[JD_loadings[n]]  +
            unique_rep[JD_loadings[n]]*sqrt(1 - nu[dd[n]]^2 )*theta_sigma[JD_loadings[n]] - beta[ii[n]])));
      log_lik[n] = bernoulli_logit_lpmf(y[n] | alpha[ii[n]] * (theta_mean[JD_loadings[n]] +
                    z_theta[jj[n]] * nu[dd[n]] * theta_sigma[JD_loadings[n]] +
                    unique[JD_loadings[n]] * sqrt(1 - nu[dd[n]]^2) * theta_sigma[JD_loadings[n]] - beta[ii[n]]));
    }
}")
    }

    if(irt_model == "Testlet1PLnorm"){
      block_quantities <- paste(block_quantities ,"for (n in 1:N){
    y_rep[n] = bernoulli_rng(inv_logit((theta_mean[jj[n]] +
    theta_rep[jj[n]] *nu[dd[n]] * theta_sigma[jj[n]]  +
    unique_rep[JD_loadings[n]]*sqrt(1 - nu[dd[n]]^2 )*theta_sigma[jj[n]] - beta[ii[n]])));

    log_lik[n] = bernoulli_logit_lpmf(y[n] |(theta_mean[jj[n]] +
    z_theta[jj[n]] *nu[dd[n]] * theta_sigma[jj[n]]  +
    unique[JD_loadings[n]]*sqrt(1 - nu[dd[n]]^2 )*theta_sigma[jj[n]] - beta[ii[n]]));
  }
}")
    }
    if(irt_model == "Testlet2PLnorm"){
      block_quantities <- paste(block_quantities ,"for (n in 1:N){
    y_rep[n] = bernoulli_rng(inv_logit(alpha[ii[n]]*(theta_mean[jj[n]] +
    theta_rep[jj[n]] *nu[dd[n]] * theta_sigma[jj[n]]  +
    unique_rep[JD_loadings[n]]*sqrt(1 - nu[dd[n]]^2 )*theta_sigma[jj[n]] - beta[ii[n]])));

    log_lik[n] = bernoulli_logit_lpmf(y[n] | alpha[ii[n]]*(theta_mean[jj[n]] +
    z_theta[jj[n]] *nu[dd[n]] * theta_sigma[jj[n]]  +
    unique[JD_loadings[n]]*sqrt(1 - nu[dd[n]]^2 )*theta_sigma[jj[n]] - beta[ii[n]]));
  }
}")
    }

  }

  return( block_quantities)

}

create_stan_data <- function(irt_model, age_model, prior_knowledge = NULL, parameter_fixed  = NULL){
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
                          prior_knowledge = prior_knowledge, parameter_fixed  =  parameter_fixed )
  p3 <- create_block_param(irt_model = irt_model, age_model = age_model,
                           is_multi = is_multi, is_2pl = is_2pl, is_flex_multi = is_flex_multi,
                           parameter_fixed  =  parameter_fixed)
  p4 <- create_block_transparam(irt_model = irt_model, age_model = age_model,
                                is_multi = is_multi, is_2pl = is_2pl, is_flex_multi = is_flex_multi)
  p5 <- create_block_model(irt_model = irt_model, age_model = age_model,
                           is_multi = is_multi, is_2pl = is_2pl, is_flex_multi = is_flex_multi,
                           prior_knowledge = prior_knowledge, parameter_fixed  =  parameter_fixed )
  p6 <- create_block_generatequantities(irt_model = irt_model,
                                        is_multi = is_multi, is_2pl = is_2pl,
                                        is_flex_multi = is_flex_multi)

  full_code <- paste(p2,p3,p4,p5,p6)
  if(irt_model == "3PLnorm"){
    full_code <- paste("functions {
  real logit_3pl_lpmf(int y, real theta, real alpha, real beta, real c) {
    real eta = alpha * (theta - beta);
    if (y == 1) {
      return log_sum_exp(log(c), log1m(c) + log_inv_logit(eta));
    } else {
      return log1m(c) + log1m_inv_logit(eta);
    }
  }
  real logit_3pl_prob(real theta, real alpha, real beta, real c) {
    real eta = alpha * (theta - beta);
    return c + (1 - c) * inv_logit(eta);
  }
}", full_code)
  }
  return(full_code)
}

.model_cache <- new.env(parent = emptyenv())

get_cache_dir <- function() {
  # tools::R_user_dir is the idiomatic R package cache location (R >= 4.0)
  dir <- tools::R_user_dir("IRTnorm", which = "cache")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  dir
}

get_compiled_model <- function(age_model, irt_model,prior_knowledge = NULL,parameter_fixed  = NULL) {
  key <- paste(age_model, irt_model, sep = "_")

  if(!is.null(prior_knowledge)){
    key <- paste(key,paste0("P_",paste(prior_knowledge,collapse = "")), sep = "_")
  }
  if(!is.null(parameter_fixed)){
    key <- paste(key, paste0("F_",paste(parameter_fixed,collapse = "")) , sep = "_")
  }
  # 1. Check in-memory cache first (fastest, within-session)
  if (exists(key, envir = .model_cache)) {
    return(.model_cache[[key]])
  }

  # 2. Check disk cache (persists across sessions)
  rds_path <- file.path(get_cache_dir(), paste0(key, ".rds"))
  if (file.exists(rds_path)) {
    message("Loading cached model from disk: ", key)
    model <- readRDS(rds_path)
    .model_cache[[key]] <- model   # populate in-memory cache too
    return(model)
  }

  # 3. Compile fresh, then save to both caches
  message("Compiling Stan model for: ", key)
  stan_code <- create_stan_data(age_model = age_model, irt_model = irt_model,
                                prior_knowledge = prior_knowledge,
                                parameter_fixed = parameter_fixed)
  tmp <- tempfile(fileext = ".stan")
  writeLines(stan_code, tmp)

  # --- your compilation step here ---
  model <- cmdstanr::cmdstan_model(tmp,quiet = TRUE)


  saveRDS(model, rds_path)           # persist to disk
  .model_cache[[key]] <- model       # store in memory
  return(model)
}


