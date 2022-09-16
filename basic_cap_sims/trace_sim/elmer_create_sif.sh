#!/bin/bash

echo $1

for entry in "UNV"/*
do
  echo "$entry"
  ElmerGrid 8 2 "$entry"
done
