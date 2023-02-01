#!/bin/bash

for rep in {0..10}
do
  for i in {0..100}
  do
    sbatch ./simulation_singlecore.sh 
  done
  sleep 60
done

