#!/bin/bash
#SBATCH --job-name=singlecore_job 
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --requeue
#SBATCH --qos=bonus


source ~/anaconda3/etc/profile.d/conda.sh
conda activate sage

#sage everything_h5py.sage ${1} ${2} ${3} ${4} ${5}
sage everything.sage ${1} ${2} ${3} ${4} ${5}

