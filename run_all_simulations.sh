#!/bin/bash

atatime=400
reps=40
for rep in $(seq 1 $reps)
do
  for i in $(seq 0 $atatime)
  do
    j=$(( i+atatime*rep ))
    echo $j
    sbatch ./simulation_singlecore.sh $j
  done
  sleep 5000
done
