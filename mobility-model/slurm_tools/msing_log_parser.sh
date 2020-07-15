#/usr/bin/env bash

#printf "job_id\ttask_id\tdiv_tran\tlow_bfmi\tRhat\tess_bulk\tess_tail\titer1\titer2\titer3\titer4\tstatus\tcountry_name\n"

slurm_tools_dir=`dirname $0`

for f in $* 
do 
  ${slurm_tools_dir}/log_parser.awk mob_$f.log
done
