#!/bin/bash

for rep in {1..4}
do
  for i in {0..400}
  do
    j=$(( $i*$rep ))
    echo $j
    sbatch ./simulation_singlecore.sh $j
  done
  sleep 5000 
done

