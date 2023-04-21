#!/bin/bash
#SBATCH --job-name=log_transform
#SBATCH --output=log_transform_%j.log
#SBATCH --error=log_transform_%j.log
#SBATCH --mem=100GB

python log_concatenated_csv_files.py
