# get environment variables
# if you want to see commands in output file
options(echo=TRUE) 
# import args
args=(commandArgs(TRUE))
iter <- eval(parse(text = args[1]))
print(iter)



parameter_grid <- expand.grid(
  seed = 1:1000,
  sample_size = c(1000),
  skew = c('low', 'medium', 'high')
)


# load libraries
library(drtmle)
library(SuperLearner)
library(earth)
library(ranger)
library(xgboost)
library(nnet)
library(VGAM)
library(cplm)
library(pscl)
library(mgcv)

# Super Learner library
# stage-1 & stage-2 are simple since we only care about one-stage super learner
stage_1.library <- c("SL.glm")
stage_2.library <- c("SL.glm")
# For one-stage super learner
SL.library <- c("SL.glm", "SL.glmnet", "SL.earth", "SL.zip", "SL.tobit", "SL.ipredbagg",
                "SL.svm", "SL.ranger", "SL.xgboost", "SL.nnet")




#################################################################################################################
# Generate Data 
#################################################################################################################


# source in functions 
source("/projects/dbenkes/ziyue/topic_2/ate/make_data.R")
source("/projects/dbenkes/ziyue/topic_2/ate/HuberSL.R")
source("/projects/dbenkes/ziyue/topic_2/ate/SL_estimators.R")



# do your simulation for row iter of parameter_grid
# make data based on parameter_grid[iter, ]
# treating w1 as the treatment variable (0/1)

# data
data <- createData(n = parameter_grid$sample_size[iter],
                   skew = parameter_grid$skew[iter]
                   )
  
# zero percentage
zero_percent <- (sum(data$y==0) / dim(data)[1]) * 100

# outlier porportion
outlier_prop <- mean(data$y >= quantile(data$y,probs = c(0.75)) + 
                       1.5*(quantile(data$y,probs = c(0.75))-quantile(data$y,probs = c(0.25))))*100



#################################################################################################################
# Truth & Estimation of ATE
#################################################################################################################


################################################# True ATE ##################################################

# calculate true ATE numerically

n_true <- 1e6
data_1 <- createData.1(n = n_true, skew = parameter_grid$skew[iter])
data_0 <- createData.0(n = n_true, skew = parameter_grid$skew[iter])

# true ATE
true_est <- mean(data_1$y - data_0$y)


########################################## Naive Estimator ################################################

# naive ATE estimate
naive_est <- mean(data$y[data$A==1]) - mean(data$y[data$A==0])




#========================================================================#
# For standard super learner and huber loss-based super learner, 
# we only change the Outcome Regression (OR) formula for them and 
# keep the rest (g, gr, Qr) the same for comparison
#========================================================================#

############################################## Super Learner ##############################################

# Fit drtmle with super learner
sl_fit <- drtmle(W = data[,-c(1,2,12)], 
                 A = data$A, 
                 Y = data$y, 
                 family = gaussian(),
                 SL_Q = SL.library, # Only change this argument
                 glm_g = ".",
                 glm_Qr = "gn", 
                 glm_gr = "Qn",
                 a_0 = c(1,0), # specify treatment order
                 stratify = FALSE,
                 returnModels = TRUE)

# ATE estimate
standard_SL_est <- sl_fit$tmle$est[1] - sl_fit$tmle$est[2]


###################################### Huber loss-based Super Learner #####################################

# Arbitrary user-specified regressions in drtmle
# get max y
max_y <- max(data$y)

# create robustification parameter lambda
lam_candidate <- c(0.00001, 0.000025, 0.00005, 0.000075, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 
         0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1, 2.5, 
         5, 7.5, 10, 25, 50, 75, 100, 250, 500, 1000, 2500, 5000, 10000)

lam <- lam_candidate[lam_candidate <= max_y/10000]

# fit a two-stage huber loss-based superlearner for outcome regression outside of drtmle
huber_sl_fit <- HuberSL(Y = data$y, 
                        X = data[,-c(1,12)], 
                        library.2stage =list(stage1=stage_1.library,
                                             stage2=stage_2.library),
                        library.1stage = SL.library,
                        lambda = lam,
                        cvControl = list(V = 10))

# One-stage huber loss-based super learner with partial-CV
# predicted values setting w1 = 1
pred_1.partial <- predict.onestage.huber_SL.partial_CV(huber_sl_fit, newdata = data.frame(A = 1, data[,-c(1,2,12)]))
# predicted values setting w1 = 0
pred_0.partial <- predict.onestage.huber_SL.partial_CV(huber_sl_fit, newdata = data.frame(A = 0, data[,-c(1,2,12)]))

# One-stage huber loss-based super learner with nested-CV
# predicted values setting w1 = 1
pred_1.nested <- predict.onestage.huber_SL.nested_CV(huber_sl_fit, newdata = data.frame(A = 1, data[,-c(1,2,12)]))
# predicted values setting w1 = 0
pred_0.nested <- predict.onestage.huber_SL.nested_CV(huber_sl_fit, newdata = data.frame(A = 0, data[,-c(1,2,12)]))

# generate Qn list
Qn10.partial <- list(
  # first entry is predicted values setting w1 = 1
  pred_1.partial$pred,
  # second entry is predicted values setting w1 = 0
  pred_0.partial$pred
)

# generate Qn list
Qn10.nested <- list(
  # first entry is predicted values setting w1 = 1
  pred_1.nested$pred,
  # second entry is predicted values setting w1 = 0
  pred_0.nested$pred
)

# pass this list to drtmle to avoid internal estimation of 
# outcome regression (note propensity score and reduced-dimension
# regressions are still estimated internally)

huber_fit.partial <- drtmle(W = data[,-c(1,2,12)], 
                            A = data$A, 
                            Y = data$y, 
                            family = gaussian(),
                            Qn = Qn10.partial, # Only change this argument
                            glm_g = ".", 
                            glm_Qr = "gn", 
                            glm_gr = "Qn", 
                            # crucial to set a_0 to match Qn's construction!
                            a_0 = c(1,0))

huber_fit.nested <- drtmle(W = data[,-c(1,2,12)], 
                           A = data$A, 
                           Y = data$y, 
                           family = gaussian(),
                           Qn = Qn10.nested, # Only change this argument
                           glm_g = ".", 
                           glm_Qr = "gn", 
                           glm_gr = "Qn", 
                           # crucial to set a_0 to match Qn's construction!
                           a_0 = c(1,0))


# ATE estimate
huber_SL_est.partial <- huber_fit.partial$tmle$est[1] - huber_fit.partial$tmle$est[2]

huber_SL_est.nested <- huber_fit.nested$tmle$est[1] - huber_fit.nested$tmle$est[2]


#################################################################################################################
# Store ATE estimates
#################################################################################################################


# Store results
results <- c(
  ############# Data Summary ############
  zero_percent,
  outlier_prop,
  
  ########################### Truth ##########################
  true_est,
  
  ####################### Naive estimator ####################
  naive_est,
  
  ################## Standard Super Learner ##################
  standard_SL_est,
  
  ############## Huber loss-based super learner ##############
  huber_SL_est.partial,
  huber_SL_est.nested
)


# save output
save(results, file=paste0("/projects/dbenkes/ziyue/topic_2/ate/results/result_n=",parameter_grid$sample_size[i],
                         "_outlier=",parameter_grid$skew[i],
                         "_seed=", parameter_grid$seed[i],".RData"))
