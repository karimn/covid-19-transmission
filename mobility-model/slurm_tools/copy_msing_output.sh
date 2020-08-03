#/usr/bin/env bash

slurm_tools_dir=`dirname $0`
first_job_id=$1
last_job_id=$2
destination=$3

for f in `seq ${first_job_id} ${last_job_id}`
do
  fit_file=multi_$f_mob.RData
  results_file=multi_$f_mob_results.RData

  if [ -e $fit_file ]
  then
    cp $fit_file $destination
  fi

  if [ -e $results_file ]
  then
    cp $results_file $destination
  fi
done
