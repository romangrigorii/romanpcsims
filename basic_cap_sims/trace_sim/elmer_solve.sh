#!/bin/bash

echo $1

for entry in "UNV"/*/
do
  echo "$entry"
  cd "$entry"
  result=${PWD##*/}
  echo "Solving $result.sif"
  ElmerSolver "$result.sif"
  cd ..
  cd ..
done
