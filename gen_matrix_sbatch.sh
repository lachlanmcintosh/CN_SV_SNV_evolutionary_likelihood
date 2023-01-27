#!/bin/bash
#SBATCH --job-name=singlecore_job 
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --requeue

sbatch gen_matrix_sbatch_solo.sh #sage matrix_gen.sage &

source ~/anaconda3/etc/profile.d/conda.sh
conda activate sage

# n,l 40,20

n=32 # this is the max copy number allowed in the model
l=128 # this is the maximum mnumber of paths in the model
p="c32_p128_v1" # the path breakdown

for i in 1 2 3 4 5 6 7 8
do
  ./run_GD.sh $n
  cd precomputed
  ./run_allpaths.sh $n $l $p
  cd ..
  ./run_collate.sh $n $l $p
  sleep 3600
done
