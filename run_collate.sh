#!/bin/bash

n=$1
l=$2
p=$3

for i in {0..100}
do
  for j in {0..100}
  do
    FILE="collate_"$i"_"$j"_"$n"_"$p".csv" 
    echo $FILE
    if [ -f $FILE ]
    then
      echo "exists."
    else
      if [ $(($i+$j)) -lt 101 ]
      then
        sbatch ./singlecore_collate.sh $i $j $n $l $p
      fi
      echo $i
      echo $j
    fi
  done  
done

