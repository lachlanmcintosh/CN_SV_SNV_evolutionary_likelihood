#!/bin/bash

sbatch touch_all.sh

atatime=40
reps=1
total_jobs=$((atatime * reps))

# Submit the array job
sbatch --array=1-${total_jobs} ./simulation_singlecore.sh

# Sleep between repetitions if required
for rep in $(seq 2 $reps)
do
  sleep 7200
  start_index=$((atatime * (rep - 1) + 1))
  end_index=$((atatime * rep))
  sbatch --array=${start_index}-${end_index} ./simulation_singlecore.sh
  #sbatch --array=$((atatime * (rep - 1) + 1))-${total_jobs} ./simulation_singlecore.sh
done
