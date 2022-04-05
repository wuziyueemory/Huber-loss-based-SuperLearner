# ----------------------------------
# Merge results
# ----------------------------------
options(echo=TRUE) # if you want to see commands in output file

## simulation parameters
parameter_grid <- expand.grid(
  seed = 1:1000,
  sample_size = c(1000),
  skew = c('low', 'medium', 'high')
)

## result name
result.name <- c(
  ############# Data Summary ############
  "zero percentage",
  "outlier proportion",
  
  ########################### Truth ##########################
  "true ATE",
  
  ####################### Naive estimator ####################
  "naive ATE",
  
  ################## Standard Super Learner ##################
  "standard SL ATE",
  
  ############## Huber loss-based super learner ##############
  "huber SL ATE_partial",
  "huber SL ATE_nested"
)


result.num <- length(result.name)

## Store results
overall <- matrix(NA, nrow = nrow(parameter_grid), ncol = result.num)

for(i in 1:nrow(parameter_grid)){
  tmp <- get(load(paste0("/projects/dbenkes/ziyue/topic_2/ate/results/result_n=",parameter_grid$sample_size[i],
                         "_outlier=",parameter_grid$skew[i],
                         "_seed=", parameter_grid$seed[i],".RData")))
  overall[i,] <- tmp
}
overall <- cbind(parameter_grid,overall)
colnames(overall)[4:10] <- result.name



## Save result object
save(overall, 
     file=paste0("/projects/dbenkes/ziyue/topic_2/ate/n=",parameter_grid$sample_size[1],".RData"))
