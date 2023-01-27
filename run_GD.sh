#!/bin/bash

n=$1
for i in {0..100}
do
  for j in {0..100}
  do
    FILE="precomputed/pre_mat"$n"_u"$i"_d"$j".sobj" 
    #echo $FILE
    if [ -f $FILE ]
    then
      echo "exists."
    else
      if [ $(($i+$j)) -lt 101 ]
      then
        sbatch ./singlecore.sh $i $j $n
      fi
      echo $i
      echo $j
    fi
  done  
done

