#!/bin/bash

COUNTRIES=`Rscript get_nat_with_subnat.R`
#COUNTRIES=28,43,32

echo "Running countries: ${COUNTRIES}"

sbatch --array=$COUNTRIES separate_countries_slurm.sh 
