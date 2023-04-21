#!/bin/bash
#SBATCH --job-name=singlecore_job 
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --requeue

sage matrix_gen.sage ${1}

