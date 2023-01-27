#!/bin/bash

#sbatch gen_matrix_sbatch_solo.sh #sage matrix_gen.sage &

n=129 # this is the max copy number allowed in the model
l=128 # this is the maximum mnumber of paths in the model
p="p128_v3" # the path breakdown

for rep in {0..0}
do
  for i in {0..100}
  do
    for j in {0..100}
    do
      if [ $(($i+$j)) -lt 101 ]
      then
        FILE="precomputed/collated_"$i"_"$j"_"$n"_"$p".csv" 
        if ! [ -f $FILE ]
        then
          #echo $FILE
          #echo "exists."
        #else
          echo $FILE
          echo "exists."
          sbatch ./everything_singlecore.sh $i $j $n $l $p
        fi
      fi
      echo $i
      echo $j
    done  
    #sleep 10
  done
done

