################################################################################
# MEPS data analysis
################################################################################

# load libraries
library(quadprog)
library(SuperLearner)
library(pscl) 
library(VGAM) 
library(cplm)
library(caret) 
library(randomForest) 
library(e1071) 
library(moments) 
library(flexsurv) 
library(quantreg) 
library(sandwich)
library(CVXR)
library(scales)
library(nloptr)
library(spaMM)

# source functions 
source("/projects/dbenkes/ziyue/topic_2/MEPS/meps_code.R")

# Loading MEPS train_data
train <- get(load(paste0("/projects/dbenkes/ziyue/topic_2/MEPS/train.RData")))
test <- get(load(paste0("/projects/dbenkes/ziyue/topic_2/MEPS//test.RData")))


###################################################################################
# Fit Huber loss-based Super Learner 
###################################################################################
# create threshods
lam <- c(0.0001, 0.0005, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1, 2.5, 5, 7.5, 10, 
         12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35)
         
lamname <- paste('lambda =',10000*lam)


# calculate zero percentage
zero_train <- (sum(train$y==0) / dim(train)[1]) * 100
zero_test <- (sum(test$y==0) / dim(test)[1]) * 100

# calculate outlier proportion
prop_train <- mean(train$TOTEXP >= quantile(train$TOTEXP,probs = c(0.75)) + 
                     1.5*(quantile(train$TOTEXP,probs = c(0.75))-quantile(train$TOTEXP,probs = c(0.25))))*100

prop_test <- mean(test$TOTEXP >= quantile(test$TOTEXP,probs = c(0.75)) + 
                     1.5*(quantile(test$TOTEXP,probs = c(0.75))-quantile(test$TOTEXP,probs = c(0.25))))*100


# fit two-stage superlearner with different loss function    
twostage.fit <- HuberSL(Y = train$TOTEXP, 
                        X = train[,-c(1,21)], 
                        newX = test[,-c(1,21)],
                        library.2stage = list(stage1=c("SL.glm","SL.glmnet","SL.knn","SL.rpart"),
                                              stage2=c("SL.logOLS.smear", "SL.gammaIdentityGLM","SL.glmnet","SL.rf.caret1")),
                        library.1stage = c("SL.lm", "SL.glmnet", "SL.rf.caret1"),
                        lambda = lam,
                        cvControl = list(V = 10))

###################################################################################
# Obtain predictions 
###################################################################################
  
  # Two-stage
  # Discrete learner (MSE vs. HUBER)
  discrete.LS.two <- cost.fit$SL.predict$`two-stage discrete learner`$'Standard (MSE)'
  # partial-CV
  discrete.HUBER.two.partial <- cost.fit$SL.predict$`two-stage discrete learner`$'Huber-partial CV'
  # nested-CV
  discrete.HUBER.two.nested <- cost.fit$SL.predict$`two-stage discrete learner`$'Huber-nested CV'
  
  # Super Learner (MSE vs. HUBER)
  SL.LS.two <- cost.fit$SL.predict$`two-stage super learner`$'Standard (MSE)'
  # partial-CV
  SL.HUBER.two.partial <- cost.fit$SL.predict$`two-stage super learner`$'Huber-partial CV'
  # nested-CV
  SL.HUBER.two.nested <- cost.fit$SL.predict$`two-stage super learner`$'Huber-nested CV'
  
  # One-stage
  # Discrete learner (MSE vs. HUBER)
  discrete.LS.one <- cost.fit$SL.predict$`one-stage discrete learner`$'Standard (MSE)'
  # partial-CV
  discrete.HUBER.one.partial <- cost.fit$SL.predict$`one-stage discrete learner`$'Huber-partial CV'
  # nested-CV
  discrete.HUBER.one.nested <- cost.fit$SL.predict$`one-stage discrete learner`$'Huber-nested CV'
  
  # Super Learner (MSE vs. HUBER)
  SL.LS.one <- cost.fit$SL.predict$`one-stage super learner`$'Standard (MSE)'
  # partial-CV
  SL.HUBER.one.partial <- cost.fit$SL.predict$`one-stage super learner`$'Huber-partial CV'
  # nested-CV
  SL.HUBER.one.nested <- cost.fit$SL.predict$`one-stage super learner`$'Huber-nested CV'


###################################################################################
## get prediction performance (MSE vs. HUBER)
###################################################################################

  # MSE
  mse <- c(mean((test$TOTEXP - discrete.LS.two)^2),
           mean((test$TOTEXP - discrete.HUBER.two.partial)^2),
           mean((test$TOTEXP - discrete.HUBER.two.nested)^2),
           
           mean((test$TOTEXP - SL.LS.two)^2),
           mean((test$TOTEXP - SL.HUBER.two.partial)^2),
           mean((test$TOTEXP - SL.HUBER.two.nested)^2),
           
           mean((test$TOTEXP - discrete.LS.one)^2),
           mean((test$TOTEXP - discrete.HUBER.one.partial)^2),
           mean((test$TOTEXP - discrete.HUBER.one.nested)^2),
           
           mean((test$TOTEXP - SL.LS.one)^2),
           mean((test$TOTEXP - SL.HUBER.one.partial)^2),
           mean((test$TOTEXP - SL.HUBER.one.nested)^2))
  
  # MAE
  mAe <- c(mean(abs(test$TOTEXP - discrete.LS.two)),
           mean(abs(test$TOTEXP - discrete.HUBER.two.partial)),
           mean(abs(test$TOTEXP - discrete.HUBER.two.nested)),
           
           mean(abs(test$TOTEXP - SL.LS.two)),
           mean(abs(test$TOTEXP - SL.HUBER.two.partial)),
           mean(abs(test$TOTEXP - SL.HUBER.two.nested)),
           
           mean(abs(test$TOTEXP - discrete.LS.one)),
           mean(abs(test$TOTEXP - discrete.HUBER.one.partial)),
           mean(abs(test$TOTEXP - discrete.HUBER.one.nested)),
           
           mean(abs(test$TOTEXP - SL.LS.one)),
           mean(abs(test$TOTEXP - SL.HUBER.one.partial)),
           mean(abs(test$TOTEXP - SL.HUBER.one.nested))
         )
  # R sqaure
  Rsq <-  1 - mse/(var(test$TOTEXP))
  
  mse <- c(zero_train, zero_test, prop_train, prop_test, mse)
  mae <- c(zero_train, zero_test, prop_train, prop_test, mae)
  Rsq <- c(zero_train, zero_test, prop_train, prop_test, Rsq)
  
  # save output
  save(mse, file=paste0("/projects/dbenkes/ziyue/topic_2/MEPS/MSE_n=",parameter_grid$sample_size[iter],
                        "_skew=",parameter_grid$skew[iter],
                        "_seed=", parameter_grid$seed[iter],".RData"))
  
  save(mae, file=paste0("/projects/dbenkes/ziyue/topic_2/MEPS/MAE_n=",parameter_grid$sample_size[iter],
                        "_skew=",parameter_grid$skew[iter],
                        "_seed=", parameter_grid$seed[iter],".RData"))
  save(Rsq, file=paste0("/projects/dbenkes/ziyue/topic_2/MEPS/Rsq_n=",parameter_grid$sample_size[iter],
                        "_skew=",parameter_grid$skew[iter],
                        "_seed=", parameter_grid$seed[iter],".RData"))





