__README for code used in Monte Carlo simulation about cost prediction__



# Overview 

This folder includes the materials for Monte Carlo simulation regarding the evaluation of Huber loss-based super learner in terms of cost prediction 
under different data generating settings. In particular, for a fixed zero-percentage (∼35%), we were interested in comparing performance of 
Huber loss-based super learner under varying sample sizes and varying outlier proportions. We considered settings with low (∼3.5%), medium(∼10%), 
and high (∼20%) outlier proportions and compared performance of methods at sample sizes of 250, 500, 1000, and 2000.

We considered four candidate algorithms for the super learner: ordinary least squares (OLS), lasso, Support Vector Machine (SVM) and Random Forest. 
These four algorithms were used to define six super learners: standard discrete learner (based on optimizing MSE), standard super learner 
ensemble (based on optimizing MSE), partial-CV Huber loss-based discrete learner, nested-CV Huber loss-based discrete learner, 
partial-CV Huber loss-based super learner ensemble, and nested-CV Huber loss-based super learner ensemble. Cost predictions were evaluated using MSE,
MAE and R^2.

Specifically, this folder includes the code for running cost prediction `simulation.R`, functions of Huber loss-based super learner `HuberSL.R`, 
functions for generating simulated data `createData.R`, functions for candidate algorithms in the super learner library `SL_estimators.R`, 
functions for merging simulation results `merge.R`, and shell script for executing simulation `simulation.sh` and merging results `merge.sh`.

The code for this simulation were run on a Linux system with a Slurm Workload Manager. The scripts can be executed with different commands passed 
in from the system environment. We have also included a shell script (`simulation.sh`) that shows the workflow of executing the R scripts. However, 
note that some of the options for batching the jobs are native to the host system they were executed on and thus will error if executed on other 
systems. To read more about this work flow see [this page](https://github.com/FredHutch/slurm-examples/tree/master/centipede). 


# Questions

There is a substantial amount of code associated with this project and
there was a significant amount of computational burden in executing the
production size jobs. The code submitted along with the manuscript needs 
to be modified in some places in order to ease its use on different systems. 
I have not debugged the code across all different systems and I 
cannot guarantee that the code will run error-free on new systems. If you come 
across any issues, [please reach out to me by email](ziyue.wu@emory.edu) 
and I am happy to help sort them out. 
