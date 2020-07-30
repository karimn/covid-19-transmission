#/usr/bin/env bash

slurm_tools_dir=`dirname $0`
first_job_id=$1
last_job_id=$2

for f in `seq ${first_job_id} ${last_job_id}`
do
  log_file=mob_$f.log
  if [ -e $log_file ]
  then
    ${slurm_tools_dir}/log_parser.awk $log_file
  fi
done
