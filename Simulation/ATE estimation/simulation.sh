#!/bin/bash

##################### Change these constants ##############################
path_read="/projects/dbenkes/ziyue/topic_2/ate/"  # load directory

# performing 1000 simulations

# if path directory doesn't exist, make it
[ ! -d ${path_read}results ] && mkdir ${path_save}results
[ ! -d ${path_read}out ] && mkdir ${path_read}out
[ ! -d ${path_read}err ] && mkdir ${path_read}err



for ((i = 1; i <= 3000; i ++ ))
do

echo "#!/bin/bash" >> script.out
echo "#SBATCH --partition=short-cpu" >> script.out
echo "#SBATCH --job-name=ate_${i}" >> script.out
echo "#SBATCH --mem=16G" >> script.out
echo "#SBATCH --error=${path_read}err/1000_${i}.err" >> script.out
echo "#SBATCH --output=${path_read}out/1000_${i}.out" >> script.out
echo "#SBATCH --mail-user=zwu56@emory.edu." >> script.out
echo "module purge" >> script.out
echo "module load R/4.0.3" >> script.out
echo "R CMD BATCH \"--args i=${i}\" ${path_read}sim.R" >> script.out

chmod +x script.out
sbatch ./script.out
rm -rf script.out

done

