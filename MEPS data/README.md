__README for MEPS data analysis code__

# Overview 
In this folder we includes codes for conducting real data analysis with data from the Medical Expenditure Panel Survey (MEPS).
MEPS is a national survey on the financing and use of medical care of families and individuals, 
their medical providers, and employers across the United States. We used the longitudinal data of MEPS 2016-2017, with 2016 MEPS used as a 
training sample and the same individuals in 2017 MEPS used as a testing sample. We developed a prediction model for total annual healthcare expenditures based on demographics, medical conditions, and insurance characteristics.

We have also included a shell script (`meps.sh`) that shows the 
workflow of executing the R scripts `meps.R`. However, note that some of 
the options for batching the jobs are native to the host system 
they were executed on and thus will error if executed on other 
systems. To read more about this work flow see 
[this page](https://github.com/FredHutch/slurm-examples/tree/master/centipede). 

This folder includes the training dataset `train.RData` which contains data from MEPS 2016,
testing dataset `test.RData` which contains data from MEPS 2017, R script for MEPS 
data analysis `meps.R`, functions for Huber loss-based super learner `HuberSL.R`, functions 
for candidate estimators `SL_estimators.R` and shell script for executing MEPS data analysis 
R scripts `meps.sh`.
