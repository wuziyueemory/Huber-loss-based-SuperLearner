# get environment variables
# if you want to see commands in output file
options(echo=TRUE) 
# import args
args=(commandArgs(TRUE))
iter <- eval(parse(text = args[1]))
print(iter)

# simulation parameters
parameter_grid <- expand.grid(
  seed = 1:1000,
  skew = c('low', 'medium', 'high),
  sample_size = c(250, 500, 1000, 2000)
)

################## execute job ##################
  
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
  
  # source in functions 
  source("/projects/dbenkes/ziyue/topic_2/cost/createData.R")
  source("/projects/dbenkes/ziyue/topic_2/cost/Estimators.R")
  source("/projects/dbenkes/ziyue/topic_2/cost/HuberSL.R")

  # create candidate values for robustification parameter lambda
  lam <- c(0.00001, 0.000025, 0.00005, 0.000075, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 
           0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1, 2.5, 
           5, 7.5, 10, 25, 50, 75, 100)

  lamname <- paste('lambda =',10000*lam)
  
  # do your simulation for row iter of parameter_grid
  # make training & testing data based on parameter_grid[iter,]
  train <- createData(n = parameter_grid$sample_size[iter],
                      skew = parameter_grid$skew[iter]
  )
  
  test <- createData(n = 10000,
                     skew = parameter_grid$skew[iter]
  )
  
  # zero percentage
  zero_train <- (sum(train$y==0) / dim(train)[1]) * 100
  zero_test <- (sum(test$y==0) / dim(test)[1]) * 100

  # outlier proportion
  prop_train <- mean(train$y >= quantile(train$y,probs = c(0.75)) + 
                       1.5*(quantile(train$y,probs = c(0.75))-quantile(train$y,probs = c(0.25))))*100
  
  prop_test <- mean(test$y >= quantile(test$y,probs = c(0.75)) + 
                      1.5*(quantile(test$y,probs = c(0.75))-quantile(test$y,probs = c(0.25))))*100
  
  # fit two-stage superlearner with different loss function (MSE vs. HUBER)
  cost.fit <- HuberSL(Y = train$y, X = train[,-c(1,12)], newX = test[,-c(1,12)], 
                          library.2stage =list(stage1=c("SL.glm","SL.knn","SL.randomForest","SL.glmnet"),
                                               stage2=c("SL.gammaLogGLM","SL.glmnet","SL.lm","SL.randomForest")),
                          library.1stage = c("SL.lm","SL.glmnet","SL.svm","SL.randomForest"),
                          lambda = lam,
                          cvControl = list(V = 10))
  
  ## Obtain predictions
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
  



  ## get performance evaluation (MSE vs. HUBER)
  # MSE
  mse <- c(mean((test$y - discrete.LS.two)^2),
           mean((test$y - discrete.HUBER.two.partial)^2),
           mean((test$y - discrete.HUBER.two.nested)^2),
           
           mean((test$y - SL.LS.two)^2),
           mean((test$y - SL.HUBER.two.partial)^2),
           mean((test$y - SL.HUBER.two.nested)^2),
           
           mean((test$y - discrete.LS.one)^2),
           mean((test$y - discrete.HUBER.one.partial)^2),
           mean((test$y - discrete.HUBER.one.nested)^2),
           
           mean((test$y - SL.LS.one)^2),
           mean((test$y - SL.HUBER.one.partial)^2),
           mean((test$y - SL.HUBER.one.nested)^2))
  
  # MAE
  mAe <- c(mean(abs(test$y - discrete.LS.two)),
           mean(abs(test$y - discrete.HUBER.two.partial)),
           mean(abs(test$y - discrete.HUBER.two.nested)),
           
           mean(abs(test$y - SL.LS.two)),
           mean(abs(test$y - SL.HUBER.two.partial)),
           mean(abs(test$y - SL.HUBER.two.nested)),
           
           mean(abs(test$y - discrete.LS.one)),
           mean(abs(test$y - discrete.HUBER.one.partial)),
           mean(abs(test$y - discrete.HUBER.one.nested)),
           
           mean(abs(test$y - SL.LS.one)),
           mean(abs(test$y - SL.HUBER.one.partial)),
           mean(abs(test$y - SL.HUBER.one.nested))
         )
  # R sqaure
  Rsq <-  1 - mse/(var(test$y))

  
  mse <- c(zero_train, zero_test, prop_train, prop_test, mse)
  mae <- c(zero_train, zero_test, prop_train, prop_test, mae)
  Rsq <- c(zero_train, zero_test, prop_train, prop_test, Rsq)

  # save output
  save(mse, file=paste0("/projects/dbenkes/ziyue/topic_2/cost/results/MSE_n=",parameter_grid$sample_size[iter],
                        "_skew=",parameter_grid$skew[iter],
                        "_seed=", parameter_grid$seed[iter],".RData"))
  
  save(mae, file=paste0("/projects/dbenkes/ziyue/topic_2/cost/results/MAE_n=",parameter_grid$sample_size[iter],
                        "_skew=",parameter_grid$skew[iter],
                        "_seed=", parameter_grid$seed[iter],".RData"))

  save(Rsq, file=paste0("/projects/dbenkes/ziyue/topic_2/cost/results/Rsq_n=",parameter_grid$sample_size[iter],
                        "_skew=",parameter_grid$skew[iter],
                        "_seed=", parameter_grid$seed[iter],".RData"))
