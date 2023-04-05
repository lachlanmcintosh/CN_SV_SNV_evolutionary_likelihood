#!/bin/bash
#SBATCH --job-name=singlecore_job
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --requeue
#SBATCH --qos=bonus
#SBATCH --array=1-<MAX_INDEX>
#SBATCH --output=slurm-%j_%a.out

source ~/anaconda3/etc/profile.d/conda.sh
conda activate sage

python run_simulation_and_analysis2.py $SLURM_ARRAY_TASK_ID


