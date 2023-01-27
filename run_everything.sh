#!/bin/bash

sbatch gen_matrix_sbatch_solo.sh #sage matrix_gen.sage &

n=33 # this is the max copy number allowed in the model
l=128 # this is the maximum mnumber of paths in the model
p="p128_v3" # the path breakdown

for rep in {0..1}
do
  for i in {0..100}
  do
    for j in {0..100}
    do
      if [ $(($i+$j)) -lt 101 ]
      then
        FILE="collated_"$i"_"$j"_"$n"_"$p".csv" 
        if [ -f $FILE ]
        then
          echo $FILE
          echo "exists."
        else
          sbatch ./everything_singlecore_33.sh $i $j $n $l $p
        fi
      fi
      echo $i
      echo $j
    done  
    sleep 60 
  done
done

