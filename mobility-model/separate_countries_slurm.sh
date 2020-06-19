#!/usr/bin/env bash

#SBATCH -n 4                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-04:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared           # serial_requeue   # Partition to submit to
#SBATCH --mem=8000 # 8GB    # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../temp/log/mob_%A_%a.log  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e ../temp/log/mob_%A_%a.log  # File to which STDERR will be written, %j inserts jobid

module load gcc/9.2.0-fasrc01 R_core/3.6.3-fasrc01

#if [ ! -d "../temp/log/job_%A" ]
#then
#  mkdir "../temp/log/job_%A"
#fi

Rscript run_mob.R fit ${SLURM_ARRAY_TASK_ID} -i 2000 -o "{all_country_codes}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_mob" --hyperparam=separate_hyperparam.yaml --show-script-options

