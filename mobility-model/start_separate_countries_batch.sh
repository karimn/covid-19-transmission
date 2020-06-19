#!/usr/bin/env bash
#
# Usage: start_separate_countries_batch.sh [<countries>]
#

PATH=../docopts:$PATH

source ../docopts/docopts.sh --auto "$@"

if [ ! -z "${ARGS[<countries>]}" ] ; then
  COUNTRIES="${ARGS[<countries>]}"
else
  COUNTRIES=`Rscript get_nat_with_subnat.R`
fi

echo "Running countries: ${COUNTRIES}"

sbatch --array=$COUNTRIES separate_countries_slurm.sh 
