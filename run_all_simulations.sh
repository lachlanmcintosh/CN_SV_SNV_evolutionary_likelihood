#!/bin/bash

for rep in {0..0}
do
  for i in {0..10}
  do
    sbatch ./simulation_singlecore.sh $i
  done
  sleep 60
done

