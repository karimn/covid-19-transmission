#!/usr/bin/env bash

#SBATCH -n 4                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-12:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared           # serial_requeue   # Partition to submit to
#SBATCH --mem=12000 # 12GB    # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../temp/log/mob_%j.log  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e ../temp/log/mob_%j.log  # File to which STDERR will be written, %j inserts jobid

module load gcc/9.2.0-fasrc01 R_core/3.6.3-fasrc01

run_type=$1 # fit vs prior prediction
run_suffix=$2
iter=$3
output_args="-o multi_${SLURM_JOB_ID}_${run_suffix} --output-dir=${SCRATCH}/kremer_lab/karimn/mob_results/joint_runs"
hierarch_args="--complete-pooling=trend"
mob_model_type="--mobility-model-type=${4}"
epi_config="--epidemic-cutoff=${5}"
stan_controls="--adapt-delta=0.99"

shift 5
countries=$@

Rscript run_mob.R $run_type $countries -i $iter --hyperparam=joint_hyperparam.yaml --job-id=${SLURM_JOB_ID} --show-script-options --include-param-trend $output_args $hierarch_args $mob_model_type $epi_config $stan_controls --fixed-ifr
