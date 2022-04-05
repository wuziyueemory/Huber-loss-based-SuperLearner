#!/bin/bash

echo "#SBATCH --mail-user=zwu56@emory.edu." >> script.out
echo "#SBATCH --partition=short-cpu" >> script.out
echo "#SBATCH --job-name=merge" >> script.out
echo "#SBATCH --output=/projects/dbenkes/ziyue/topic_2/ate/merge.out" >> script.out
echo "#SBATCH --error=/projects/dbenkes/ziyue/topic_2/ate/merge.err" >> script.out

module purge
module load R/4.0.2
srun R --vanilla  < /projects/dbenkes/ziyue/topic_2/ate/merge.R
