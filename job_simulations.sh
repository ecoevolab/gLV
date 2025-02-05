#!/bin/env bash
#$ -N Exp01-run_01
#$ -o $JOB_NAME.log
#$ -e $JOB_NAME.error
#$ -cwd
#$ -l h_rt=30:00:00   # Set a runtime limit of 30 hours
#$ -pe openmp 10
#$ -l h_rss=30G       # Request 15GB of memory per core

# Load R module and code
module load r/4.4.0
Rscript /mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/Experiment01/sims-run_01.R


