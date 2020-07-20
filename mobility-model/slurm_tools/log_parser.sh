#/usr/bin/env bash

slurm_tools_dir=`dirname $0`

for f in $* 
do 
  ${slurm_tools_dir}/log_parser.awk $f 
done
