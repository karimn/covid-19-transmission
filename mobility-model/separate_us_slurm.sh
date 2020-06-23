#!/usr/bin/env bash

#SBATCH -n 4                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-14:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared           # serial_requeue   # Partition to submit to
#SBATCH --mem=16000 # 16GB    # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../temp/log/mob_us_%j.log  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e ../temp/log/mob_us_%j.log  # File to which STDERR will be written, %j inserts jobid

module load gcc/9.2.0-fasrc01 R_core/3.6.3-fasrc01

RUN_SUFFIX=mob

Rscript run_mob.R fit us -i 4000 -o "us_${SLURM_JOB_ID}_${RUN_SUFFIX}" --hyperparam=separate_hyperparam.yaml --show-script-options

