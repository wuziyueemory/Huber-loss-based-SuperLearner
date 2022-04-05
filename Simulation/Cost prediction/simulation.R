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
  source("/Users/ziyuewu//Desktop/project 2/simulation_2/simulation_code.R")

  # create candidate values for robustification parameter lambda
  lam <- c(0.00001, 0.000025, 0.00005, 0.000075, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 
           0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1, 2.5, 
           5, 7.5, 10, 25, 50, 75, 100)
  lamname <- paste('lambda =',10000*lam)
  
  # do your simulation for row iter of parameter_grid
  # make training & testing data based on parameter_grid[iter, ]
  for (iter in 1161:1170){
  train <- createData(n = parameter_grid$sample_size[iter],
                      skew = parameter_grid$skew[iter]
  )
  
  test <- createData(n = 10000,
                     skew = parameter_grid$skew[iter]
  )
  
  prop_train <- mean(train$y >= quantile(train$y,probs = c(0.75)) + 
                       1.5*(quantile(train$y,probs = c(0.75))-quantile(train$y,probs = c(0.25))))*100
  
  prop_test <- mean(test$y >= quantile(test$y,probs = c(0.75)) + 
                      1.5*(quantile(test$y,probs = c(0.75))-quantile(test$y,probs = c(0.25))))*100
  
  # fit two-stage superlearner with different loss function    
  twostage.fit <- twostageSL.HUBER3(Y = train$y, X = train[,-c(1,12)], newX = test[,-c(1,12)],
                                    library.2stage =list(stage1=c("SL.glm","SL.knn","SL.randomForest","SL.glmnet"),
                                                         stage2=c("SL.gammaIdentityGLM","SL.glmnet","SL.lm",
                                                                  "SL.randomForest")),
                                    library.1stage = c("SL.lm","SL.glmnet","SL.svm","SL.randomForest"),
                                    lambda = lam,
                                    cvControl = list(V = 10))
  
  # Two-stage
  # Discrete prediction (LS vs. optimal HUBER)
  discrete.LS.two <- twostage.fit$discrete.predict$twostage$`Square Loss`
  # Fake
  discrete.HUBER.two.fake <- twostage.fit$discrete.predict$twostage$`Huber Loss-optimal lambda`$fake
  # True
  discrete.HUBER.two.true <- twostage.fit$discrete.predict$twostage$`Huber Loss-optimal lambda`$true
  
  # SL prediction (LS vs. optimal HUBER)
  SL.LS.two <- twostage.fit$SL.predict$twostage$`Square Loss`
  # Fake
  SL.HUBER.two.fake <- twostage.fit$SL.predict$twostage$`Huber Loss-optimal lambda`$fake
  # True
  SL.HUBER.two.true <- twostage.fit$SL.predict$twostage$`Huber Loss-optimal lambda`$true
  
  # One-stage
  # Discrete prediction (LS vs. optimal HUBER)
  discrete.LS.one <- twostage.fit$discrete.predict$onestage$`Square Loss`
  # Fake
  discrete.HUBER.one.fake <- twostage.fit$discrete.predict$onestage$`Huber Loss-optimal lambda`$fake
  # True
  discrete.HUBER.one.true <- twostage.fit$discrete.predict$onestage$`Huber Loss-optimal lambda`$true
  
  # SL prediction (LS vs. optimal HUBER)
  SL.LS.one <- twostage.fit$SL.predict$onestage$`Square Loss`
  # Fake
  SL.HUBER.one.fake <- twostage.fit$SL.predict$onestage$`Huber Loss-optimal lambda`$fake
  # True
  SL.HUBER.one.true <- twostage.fit$SL.predict$onestage$`Huber Loss-optimal lambda`$true
  
  ## get prediction performance (LS vs. optimal HUBER)
  # MSE
  mse <- c(mean((test$y - discrete.LS.two)^2),
           mean((test$y - discrete.HUBER.two.fake)^2),
           mean((test$y - discrete.HUBER.two.true)^2),
           
           mean((test$y - SL.LS.two)^2),
           mean((test$y - SL.HUBER.two.fake)^2),
           mean((test$y - SL.HUBER.two.true)^2),
           
           mean((test$y - discrete.LS.one)^2),
           mean((test$y - discrete.HUBER.one.fake)^2),
           mean((test$y - discrete.HUBER.one.true)^2),
           
           mean((test$y - SL.LS.one)^2),
           mean((test$y - SL.HUBER.one.fake)^2),
           mean((test$y - SL.HUBER.one.true)^2))
  
  # R sqaure
  Rsq <-  1 - mse/(2*var(test$y))
  
  mse <- c(mse,prop_train,prop_test)
  Rsq <- c(Rsq,prop_train,prop_test)
  
  ## generate optimal lambda for Huber Loss from training
  # two discrete
  optimal.lambda.two.discrete.fake <- twostage.fit$huber.optimal.lambda$twostage$`Discrete Learner`$fake
  optimal.lambda.two.discrete.true <- twostage.fit$huber.optimal.lambda$twostage$`Discrete Learner`$true
  # two SL
  optimal.lambda.two.SL.fake <- twostage.fit$huber.optimal.lambda$twostage$SuperLearner$fake
  optimal.lambda.two.SL.true <- twostage.fit$huber.optimal.lambda$twostage$SuperLearner$true
  # one discrete
  optimal.lambda.one.discrete.fake <- twostage.fit$huber.optimal.lambda$onestage$`Discrete Learner`$fake
  optimal.lambda.one.discrete.true <- twostage.fit$huber.optimal.lambda$onestage$`Discrete Learner`$true
  # one SL
  optimal.lambda.one.SL.fake <- twostage.fit$huber.optimal.lambda$onestage$SuperLearner$fake
  optimal.lambda.one.SL.true <- twostage.fit$huber.optimal.lambda$onestage$SuperLearner$true
  
  mse <- c(mse,optimal.lambda.two.discrete.fake,optimal.lambda.two.discrete.true,
           optimal.lambda.two.SL.fake,optimal.lambda.two.SL.true,
           optimal.lambda.one.discrete.fake,optimal.lambda.one.discrete.true,
           optimal.lambda.one.SL.fake,optimal.lambda.one.SL.true)
  
  Rsq <- c(Rsq,optimal.lambda.two.discrete.fake,optimal.lambda.two.discrete.true,
           optimal.lambda.two.SL.fake,optimal.lambda.two.SL.true,
           optimal.lambda.one.discrete.fake,optimal.lambda.one.discrete.true,
           optimal.lambda.one.SL.fake,optimal.lambda.one.SL.true)
  
  ## generate optimal lambda for Huber Loss from testing
  two.discrete <- apply(twostage.fit$discrete.predict$twostage$`Huber Loss`, 2, function(x) mean((x-test$y)^2))
  test.two.discrete <- lamname[which.min(two.discrete)]
  two.SL <- apply(twostage.fit$SL.predict$twostage$`Huber Loss`, 2, function(x) mean((x-test$y)^2))
  test.two.SL <- lamname[which.min(two.SL)]
  one.discrete <- apply(twostage.fit$discrete.predict$onestage$`Huber Loss`, 2, function(x) mean((x-test$y)^2))
  test.one.discrete <- lamname[which.min(one.discrete)]
  one.SL <- apply(twostage.fit$SL.predict$onestage$`Huber Loss`, 2, function(x) mean((x-test$y)^2))
  test.one.SL <- lamname[which.min(one.SL)]
  
  mse <- c(mse,test.two.discrete,test.two.SL,test.one.discrete,test.one.SL)
  Rsq <- c(Rsq,test.two.discrete,test.two.SL,test.one.discrete,test.one.SL)
  
  ## generate MSE for optimal lambda for Huber Loss from testing
  mse.test.two.discrete <- as.numeric(two.discrete[which.min(two.discrete)])
  mse.test.two.SL <- as.numeric(two.SL[which.min(two.SL)])
  mse.test.one.discrete <- as.numeric(one.discrete[which.min(one.discrete)])
  mse.test.one.SL <- as.numeric(one.SL[which.min(one.SL)])
  
  Rsq.test.two.discrete <- 1 - mse.test.two.discrete/(2*var(test$y))
  Rsq.test.two.SL <- 1 - mse.test.two.SL/(2*var(test$y))
  Rsq.test.one.discrete <- 1 - mse.test.one.discrete/(2*var(test$y))
  Rsq.test.one.SL <- 1 - mse.test.one.SL/(2*var(test$y))
  
  mse <- c(mse,mse.test.two.discrete,mse.test.two.SL,mse.test.one.discrete,mse.test.one.SL)
  Rsq <- c(Rsq,Rsq.test.two.discrete,Rsq.test.two.SL,Rsq.test.one.discrete,Rsq.test.one.SL)
  
  # save output
  save(mse, file=paste0("/Users/ziyuewu//Desktop/project 2/simulation_2/MSE_n=",parameter_grid$sample_size[iter],
                        "_skew=",parameter_grid$skew[iter],
                        "_seed=", parameter_grid$seed[iter],".RData"))
  
  save(Rsq, file=paste0("/Users/ziyuewu//Desktop/project 2/simulation_2/Rsq_n=",parameter_grid$sample_size[iter],
                        "_skew=",parameter_grid$skew[iter],
                        "_seed=", parameter_grid$seed[iter],".RData"))
}
