#!/bin/bash
#SBATCH --job-name=singlecore_job 
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --requeue

source ~/anaconda3/etc/profile.d/conda.sh
conda activate sage

sage collate.sage ${1} ${2} ${3} ${4} ${5}

