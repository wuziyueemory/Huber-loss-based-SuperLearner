__README for simulation code__



# Overview 

This folder includes the materials for Monte Carlo simulation regarding the evaluation of Huber loss-based super learner in terms of average treatment effect (ATE) estimation. True ATE was calculated numerically in each setting. We consider two type of estimators for ATE, including a unadjusted sample mean difference estimator and a targeted minimum loss estimators (TMLE). A TMLE of the ATE requires two inputs -- an estimate of the probability of the intervention (the so-called propensity score) and an estimate of the conditional mean of the costs (the so-called outcome regression). Each TMLE considered in our simulation used a main term logistic regression for the propensity score. We considered three different approaches to modeling the outcome regression: a standard super learner, the partial-CV Huber loss-based super learner, and the nested-CV Huber loss-based super learner. Each super learner consisted of the same ten candidate algorithms and was based on ten-fold cross-validation. 

Specifically, this folder includes the code for running ATE estimations `simulation.R`, functions of Huber loss-based super learner `HuberSL.R`, functions for generating simulated data `make_data.R`, functions for candidate algorithms in the super learner library `SL_estimators.R`, functions for merging simulation results `merge.R`, and shell script for executing simulation `simulation.sh` and merging results `merge.sh`.

The code for this simulation were run on a Linux system with a Slurm Workload Manager. The scripts can be executed with different commands passed in from the system environment. We have also included a shell script (`simulation.sh`) that shows the workflow of executing the R scripts. However, note that some of the options for batching the jobs are native to the host system they were executed on and thus will error if executed on other systems. To read more about this work flow see [this page](https://github.com/FredHutch/slurm-examples/tree/master/centipede). 


# Questions

There is a substantial amount of code associated with this project and
there was a significant amount of computational burden in executing the
production size jobs. The code submitted along with the manuscript needs 
to be modified in some places in order to ease its use on different systems. 
I have not debugged the code across all different systems and I 
cannot guarantee that the code will run error-free on new systems. If you come 
across any issues, [please reach out to me by email](ziyue.wu@emory.edu) 
and I am happy to help sort them out. 
