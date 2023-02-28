#!/bin/bash
#SBATCH --job-name=singlecore_job 
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --requeue
#SBATCH --qos=bonus


source ~/anaconda3/etc/profile.d/conda.sh
conda activate sage

python run_simulation_and_analysis2.py $1

