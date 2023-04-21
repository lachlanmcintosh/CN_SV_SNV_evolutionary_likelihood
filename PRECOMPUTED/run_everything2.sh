#!/bin/bash

max_path_length=8

#sage make_combos_adjacent.sage "$max_path_length"
#sbatch gen_matrix_sbatch_solo.sh "$max_path_length"

max_CN=$((2 ** $max_path_length))
path_length=8
path_description="p"$max_path_length"_v4"

max_jobs=400

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
          
          num_jobs=$(squeue -u "$USER" | awk 'NR > 1' | wc -l)
          while [ $num_jobs -ge $max_jobs ]
          do
            sleep 10
            num_jobs=$(squeue -u "$USER" | awk 'NR > 1' | wc -l)
          done
          
          sbatch ./everything_singlecore.sh $p_up $p_down $max_CN $path_length $path_description
        fi
      fi
      echo $p_up
      echo $p_down
    done
  done
done

python check_all_files_exist.py


