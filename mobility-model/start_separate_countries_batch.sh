#!/bin/bash

COUNTRIES=`Rscript get_nat_with_subnat.R`

sbatch separate_countries_slurm.sh --array=$COUNTRIES
