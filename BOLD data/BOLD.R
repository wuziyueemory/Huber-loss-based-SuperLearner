################################################################################
# BOLD data analysis
################################################################################

# load libraries
library(quadprog)
library(SuperLearner)
library(randomForest) 
library(e1071) 
library(CVXR)
library(scales)
library(nloptr)
library(spaMM)
library(glmnet)
library(caret)
library(cplm)

# source functions 
source("/projects/dbenkes/ziyue/topic_2/BOLD/SL_estimators.R")
source("/projects/dbenkes/ziyue/topic_2/BOLD/HuberSL.R")

# load dataset
train <- read.csv("/projects/dbenkes/ziyue/topic_2/BOLD/train.csv")

###################################################################################
# Fit Huber loss-based Super Learner 
###################################################################################
# 10-fold cross validation
K=10
flds <- createFolds(seq(1,nrow(train)), k = K, list = TRUE, returnTrain = FALSE)

# create threshods
lam <- c(0.00001, 0.000025, 0.00005, 0.000075, 0.0001, 0.00025, 0.0005, 0.00075,
         0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 
         0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3)
         
lamname <- paste('lambda =',10000*lam)

# calculate outlier proportion
prop_train <- mean(train$TOTEXP >= quantile(train$TOTEXP,probs = c(0.75)) + 
                     1.5*(quantile(train$TOTEXP,probs = c(0.75))-quantile(train$TOTEXP,probs = c(0.25))))*100

# Generate empty matirx to store predictions
discrete.LS.one <- rep(NA, nrow(train))
discrete.HUBER.one.partial <- rep(NA, nrow(train))
discrete.HUBER.one.nested <- rep(NA, nrow(train))
SL.LS.one <- rep(NA, nrow(train))
SL.HUBER.one.partial <- rep(NA, nrow(train))
SL.HUBER.one.nested <- rep(NA, nrow(train))

# Conduct 10-fold cross-validation model fit
for (i in 1:K){
  # create train & test set
  flags <- unlist(flds[i])
  train_temp <- train[-flags, ]
  test_temp <- train[flags,]
  
  # fit Huber loss-based super learner
  bold.fit <- onestage.HuberSL(Y = train_temp$total_RVU_12moafter, 
                               X = train_temp[,-c(1,21)], 
                               newX = test_temp[,-c(1,21)],
                               SL.library = c("SL.lm","SL.glmnet","SL.gammaIdentityGLM","SL.logOLS.smear","SL.gammaLogGLM",
                                              "SL.svm","SL.randomForest", "SL.xgboost"),
                               lambda = lam,
                               cvControl = list(V = 10))

###################################################################################
# Obtain predictions 
###################################################################################
  # One-stage
  # Discrete learner (MSE vs. HUBER)
  discrete.LS.one[flags] <- bold.fit$SL.predict$`one-stage discrete learner`$'Standard (MSE)'
  # partial-CV
  discrete.HUBER.one.partial[flags] <- bold.fit$SL.predict$`one-stage discrete learner`$'Huber-partial CV'
  # nested-CV
  discrete.HUBER.one.nested[flags] <- bold.fit$SL.predict$`one-stage discrete learner`$'Huber-nested CV'
  
  # Super Learner (MSE vs. HUBER)
  SL.LS.one[flags] <- bold.fit$SL.predict$`one-stage super learner`$'Standard (MSE)'
  # partial-CV
  SL.HUBER.one.partial[flags] <- bold.fit$SL.predict$`one-stage super learner`$'Huber-partial CV'
  # nested-CV
  SL.HUBER.one.nested[flags] <- bold.fit$SL.predict$`one-stage super learner`$'Huber-nested CV'
}


###################################################################################
## get prediction performance (MSE vs. HUBER)
###################################################################################

  # MSE
  mse <- c(mean((test$total_RVU_12moafter - discrete.LS.one)^2),
           mean((test$total_RVU_12moafter - discrete.HUBER.one.partial)^2),
           mean((test$total_RVU_12moafter - discrete.HUBER.one.nested)^2),
           
           mean((test$total_RVU_12moafter - SL.LS.one)^2),
           mean((test$total_RVU_12moafter - SL.HUBER.one.partial)^2),
           mean((test$total_RVU_12moafter - SL.HUBER.one.nested)^2))
  
  # MAE
  mae <- c(mean(abs(test$total_RVU_12moafter - discrete.LS.one)),
           mean(abs(test$total_RVU_12moafter - discrete.HUBER.one.partial)),
           mean(abs(test$total_RVU_12moafter - discrete.HUBER.one.nested)),
           
           mean(abs(test$total_RVU_12moafter - SL.LS.one)),
           mean(abs(test$total_RVU_12moafter - SL.HUBER.one.partial)),
           mean(abs(test$total_RVU_12moafter - SL.HUBER.one.nested))
         )
  
# R sqaure
  Rsq <-  1 - mse/(var(test$total_RVU_12moafter))
  
  mse <- c(prop_train, mse)
  mae <- c(prop_train, mae)
  Rsq <- c(prop_train, Rsq)
  
  # save output
  save(mse, file=paste0("/projects/dbenkes/ziyue/topic_2/BOLD/MSE.RData"))
  
  save(mae, file=paste0("/projects/dbenkes/ziyue/topic_2/BOLD/MAE.RData"))
                        
  save(Rsq, file=paste0("/projects/dbenkes/ziyue/topic_2/BOLD/Rsq.RData"))


