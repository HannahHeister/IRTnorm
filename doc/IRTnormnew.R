## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----loadlibrary, echo=TRUE, include=TRUE, warning=FALSE, message = FALSE-----
library(IRTnorm)
library(bayesplot)
data("response_data")

## ----fit_IRTmod, echo=TRUE, include=TRUE, eval = FALSE------------------------
#  model_list <- list( A = list(irt_model = "1PLnorm", age_model = "poly",
#                                         poly_mean = 4, poly_sd = 3),
#                      B = list(irt_model = "2PLnorm", age_model = "poly",
#                                         poly_mean = 4, poly_sd = 3))
# 
# compare_IRTnorm_models(model_specifications = model_list,
#                                               raw_data = response_data,
#                                               age_variable = "age",
#                                               int_variables = 1:50)

## ----load_modelcomparee,  echo=FALSE, include= TRUE, eval = TRUE--------------
data("different_IRTmodels")
different_IRTmodels

## ----fit_poly, echo=TRUE, include=TRUE, eval = FALSE--------------------------
# my_grid <- expand.grid(mu = 1:4, sd = 0:3)
# different_poly_2PL <- compare_age_models(method = "grid",
#                                          grid = my_grid,
#                                          raw_data = response_data,
#                                          age_variable = "age",
#                                          int_variables = 1:50,
#                                          include_splines = FALSE)
# 

## ----load_modelcompare,  echo=FALSE, include= TRUE----------------------------
data("different_poly_2PL")

## ----load_modelcompare2,  echo=TRUE, include= TRUE, eval = TRUE,fig.width=6,fig.height=4, fig.align = "center"----
visual_polynomial_selection(different_poly_2PL, perspective = "ELPD", highlightranks = 3)

## ----load_modelcompare3,  echo=TRUE, include= TRUE, eval = TRUE,fig.width=6,fig.height=4, fig.align = "center"----
 visual_polynomial_selection(different_poly_2PL, perspective = "trajectory", highlightranks = 3)

## ----data_prep2, echo=TRUE, include=TRUE, eval = FALSE------------------------
# info <- estimate_norming_model(raw_data = response_data,
#                        age_variable = "age",
#                        int_variables = 1:50,
#                        irt_model  = "2PLnorm",
#                        age_model = "poly",
#                        poly_mean = 2,
#                        poly_sd = 2,
#                        seed = 1234,
#                        iter_warmup = 2000,
#                        iter_sampling =2000,
#                        extract_stan_object = TRUE)
# 

## ----load_info,  echo=FALSE, include= TRUE------------------------------------
data("info")

## ----showfit02, echo=TRUE, include=TRUE,eval = TRUE---------------------------
info$conv_info 

## ----showfit03, echo=TRUE, include=TRUE ,eval = FALSE-------------------------
# conv_plot <- bayesplot::mcmc_trace(info$model_stan_obj$draws(c("p1_theta")))
# #mcmc_trace(mod_poly$draws(c("p2_theta")))
# #mcmc_trace(mod_poly$draws(c("alpha")))
# #mcmc_trace(mod_poly$draws(c("beta")))

## ----showfit03_show, echo=TRUE, include=TRUE, fig.width=6,fig.height=4, fig.align = "center" ,eval = FALSE----
# data("conv_plot")
# conv_plot

## ----NC1, echo=TRUE, include=TRUE ,eval = FALSE-------------------------------
# NC_info <- NC(fit_info = info)
# NC_visual(NC_info = NC_info)

## ----NC2,  echo=FALSE, include= TRUE, fig.width=6,fig.height=4, fig.align = "center"----
#data("plot_NC")
#plot_NC
knitr::include_graphics("figures/plot_NC.png")

## ----ageMC, echo=TRUE, include=TRUE ,eval = FALSE-----------------------------
# NC_info_age <- NC(fit_info = info,
#                       age_range = seq(-2,2,1))

## ----NC4,  echo=FALSE, include= TRUE, eval = FALSE----------------------------
# data("NC_info_age")
# NC_visual(NC_info = NC_info_age[[1]])

## ----showfit2_NC, echo=FALSE, include=TRUE, fig.width=6,fig.height=4, fig.align = "center"----
knitr::include_graphics("figures/plot_NCage.png")


## ----showfit4_ITEM, echo=TRUE, include=TRUE, fig.show="hold",  fig.width=6,fig.height=4, fig.align = "center"----
item_person_fit(fit_info =info, parameter = "item")

## ----showfit4_2, echo=TRUE, include=TRUE, fig.show="hold", fig.width=6,fig.height=4, fig.align = "center"----
item_person_fit(fit_info =info, parameter = "person")

## ----rawsjpw, echo=TRUE, include=TRUE ,eval = FALSE---------------------------
# quantile_curves(fit_info = info, perspective = "rawscore",
#                             probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

## ----showfit3_raw, echo=TRUE, include=TRUE, fig.width=6,fig.height=4, fig.align = "center"----
knitr::include_graphics("figures/plot_rawscore.png")


## ----showfit5_consist, echo=TRUE, include=TRUE, fig.width=6,fig.height=4, fig.align = "center"----
quantile_curves(fit_info = info, perspective = "ability", 
                            probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

## ----showfit6_2, echo=TRUE, include=TRUE, fig.width=6,fig.height=4, fig.align = "center"----
normscore_dist(fit_info= info, age_range = seq(-2,2,1), binwidth = 0.5)

## ----showfit6, echo=TRUE, include=TRUE, fig.width=6,fig.height=3, fig.align = "center"----
effect_age(fit_info = info)

## ----newobs, echo=TRUE, include=TRUE, eval = FALSE----------------------------
# y <-  c(1,1,0,0,1,0,1,1,0,0,1,0,0,1,1,0,0,1,1,0,1,1,1,0,1,0,1,1,1,1,1,1,0,0,0,
#         1,0,1,1,0,0,0,0,1,1,1,0,1,0,0)
# age <- -1

## ----getinfo_norming, echo=TRUE, include=TRUE, eval = FALSE-------------------
# info_newperson <- norm_new_individual(new_y = y, new_age= age,
#                                      fit_info = info)
# # information about the norm score estimate of the individuals norm score
# info_newperson$summary

## ----printnorming,  echo=FALSE, include= TRUE---------------------------------
data("int_info_newperson")
int_info_newperson

