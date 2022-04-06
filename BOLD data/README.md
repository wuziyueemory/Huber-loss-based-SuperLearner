__README for BOLD data analysis__


# Overview 
In this folder we include codes for conducting real data analysis with data from the Back Pain Outcomes using Longitudinal Data (BOLD).
BOLD is a large, community-based registry of patients aged 65 years and older who presented with primary care visits for a 
new episode of back pain from March 2011 to March 2013. Expenditures in BOLD data were calculated as total relative value units (RVUs), 
a measure of value used in the US Medicare reimbursement formula for physician services. We used the total annual RVUs one year after index 
as the outcome in this study and the model was developed based on covariates from baseline patient questionnaires and EHR. 

We have also included a shell script (`bold.sh`) that shows the 
workflow of executing the R scripts `BOLD.R`. However, note that some of 
the options for batching the jobs are native to the host system 
they were executed on and thus will error if executed on other 
systems. To read more about this work flow see 
[this page](https://github.com/FredHutch/slurm-examples/tree/master/centipede). 

This folder includes the training dataset `train.csv` which contains data from BOLD (served
as both training and validation set through 10-fold cross-validation), R script for BOLD data analysis `BOLD.R`, 
functions for Huber loss-based super learner `HuberSL.R`, functions 
for candidate algorithms in super learner library `SL_estimators.R` and shell script for executing 
BOLD data analysis R scripts `BOLD.sh`.


