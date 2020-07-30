#!/usr/bin/bash

# for i in {${1}..${2}..${3}}
for (( i = ${1}; i <= $2; i += $3 ))
do
  echo Starting with $i sub-regions...
  Rscript run_mob.R fit us -i 2000 --rand-sample-subnat=$i -o speed_test &> speed_${i}.log
  echo end.
done
