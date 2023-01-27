#!/bin/bash
#SBATCH --job-name=singlecore_job 
#SBATCH --ntasks=1
#SBATCH --mem=400G
#SBATCH --requeue

module load R
R CMD BATCH ./tidy_collated.R
