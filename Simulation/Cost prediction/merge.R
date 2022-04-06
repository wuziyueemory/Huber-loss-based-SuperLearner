
# simulation parameters
parameter_grid <- expand.grid(
  seed = 1:1000,
  skew = c('low', 'medium', 'high),
  sample_size = c(250, 500, 1000, 2000)
)

# algorithm name
algo.name <- c("Zero percentage: train",
               "Zero percentage: test",
               "Outlier proportion: train",
               "Outlier proportion: test",
               "Two-stage discrete learner: MSE",
               "Two-stage discrete learner: Huber partial-CV",
               "Two-stage discrete learner: Huber nested-CV",
               "Two-stage super learner: MSE",
               "Two-stage super learner: Huber partial-CV",
               "Two-stage super learner: Huber nested-CV",
               "One-stage discrete learner: MSE",
               "One-stage discrete learner: Huber partial-CV",
               "One-stage discrete learner: Huber nested-CV",
               "One-stage super learner: MSE",
               "One-stage super learner: Huber partial-CV",
               "One-stage super learner: Huber nested-CV"
)

algo.num <- length(algo.name)

# MSE
mse_result <- matrix(NA, nrow = nrow(parameter_grid), ncol = algo.num)
for(i in 1:nrow(parameter_grid)){
  tmp <- get(load(paste0("/projects/dbenkes/ziyue/topic_2/cost/results/MSE_n=",parameter_grid$sample_size[i],
                         "_skew=",parameter_grid$skew[i],
                         "_seed=", parameter_grid$seed[i],".RData")))
  mse_result[i,] <- tmp
}
mse_result <- cbind(parameter_grid, mse_result)
colnames(mse_result)[4:33] <- algo.name

# MAE
mae_result <- matrix(NA, nrow = nrow(parameter_grid), ncol = algo.num)
for(i in 1:nrow(parameter_grid)){
  tmp <- get(load(paste0("/projects/dbenkes/ziyue/topic_2/cost/results/MAE_n=",parameter_grid$sample_size[i],
                         "_skew=",parameter_grid$skew[i],
                         "_seed=", parameter_grid$seed[i],".RData")))
  mae_result[i,] <- tmp
}
mae_result <- cbind(parameter_grid, mae_result)
colnames(mae_result)[4:33] <- algo.name

# R square
Rsq_result <- matrix(NA, nrow = nrow(parameter_grid), ncol = algo.num)
for(i in 1:nrow(parameter_grid)){
  tmp <- get(load(paste0("/projects/dbenkes/ziyue/topic_2/cost/results/Rsq_n=",parameter_grid$sample_size[i],
                         "_skew=",parameter_grid$skew[i],
                         "_seed=", parameter_grid$seed[i],".RData")))
  Rsq_result[i,] <- tmp
}
Rsq_result <- cbind(parameter_grid,Rsq_result)
colnames(Rsq_result)[4:33] <- algo.name

# save result object
save(mse_result, 
     file=paste0("/projects/dbenkes/ziyue/topic_2/cost/MSE_",parameter_grid$sample_size[1],
                 "_",parameter_grid$skew[1],".RData"))
save(mae_result, 
     file=paste0("/projects/dbenkes/ziyue/topic_2/cost/MAE_",parameter_grid$sample_size[1],
                 "_",parameter_grid$skew[1],".RData"))
save(Rsq_result, 
     file=paste0("/projects/dbenkes/ziyue/topic_2/cost/Rsq_",parameter_grid$sample_size[1],
                 "_",parameter_grid$skew[1],".RData"))
