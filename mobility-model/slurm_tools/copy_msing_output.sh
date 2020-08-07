#/usr/bin/env bash

slurm_tools_dir=`dirname $0`
first_job_id=$1
last_job_id=$2
destination=$3

for f in `seq ${first_job_id} ${last_job_id}`
do
  fit_file=multi_${f}_mob.RData
  results_file=multi_${f}_mob_results.rds

  if [ -e $fit_file ]
  then
    cp $fit_file $destination --verbose
  else 
    echo "File ${fit_file} not found"
  fi

  if [ -e $results_file ]
  then
    cp $results_file $destination --verbose
  else
    echo "File ${results_file} not found"
  fi
done
