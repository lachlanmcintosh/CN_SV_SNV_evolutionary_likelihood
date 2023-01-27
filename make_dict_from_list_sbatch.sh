#!/bin/bash
#SBATCH --job-name=singlecore_job 
#SBATCH --ntasks=1
#SBATCH --mem=100G


python make_dict_from_list.py

