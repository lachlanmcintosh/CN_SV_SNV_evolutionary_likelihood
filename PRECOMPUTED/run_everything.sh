#!/bin/bash

max_path_length=8

sage make_combos_adjacent.sage "$max_path_length"
sbatch gen_matrix_sbatch_solo.sh "$max_path_length"

max_CN=$((2 ** $max_path_length)) # this is the max copy number allowed in the model
path_length=8 # this is the maximum mnumber of paths in the model
path_description="p"$max_path_length"_v4" # the path breakdown

for rep in {0..0}
do
  for p_up in {0..100}
  do
    for p_down in {0..100}
    do
      if [ $(($p_up+$p_down)) -lt 101 ]
      then
        FILE="MATRICES/subbed_u"$p_up"_d"$p_down"_"$path_description".csv"
        if ! [ -f $FILE ]
        then
          echo $FILE
          echo "doesn't exist!"
          sbatch ./everything_singlecore.sh $p_up $p_down $max_CN $path_length $path_description
        fi
      fi
      echo $p_up
      echo $p_down
    done
    sleep 300
  done
done

