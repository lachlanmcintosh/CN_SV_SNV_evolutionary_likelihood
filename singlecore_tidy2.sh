#!/bin/bash
#SBATCH --job-name=singlecore_job 
#SBATCH --ntasks=1
#SBATCH --mem=1500G
#SBATCH --requeue

module load R
R CMD BATCH ./tidy_collated2.R
