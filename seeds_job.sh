#!/bin/env bash
#$ -N seeds
#$ -o ./Logs_errors/$JOB_NAME.log
#$ -e ./Logs_errors/$JOB_NAME.error
#$ -cwd
#$ -l h_rt=02:00:00   # Set a runtime limit of 2 hours
#$ -pe openmp 1
#$ -l h_rss=15G        # Request 2GB of memory per core

# Load the R module (if necessary)
module load r/4.4.0

# Run the R script
Rscript Seeds_generator.R

# Combine seeds data
cd /mnt/atgc-d3/sur/users/mrivera/Simulations
cat ./Seeds_01.tsv ./Seeds_02.tsv > ./combined_seeds.tsv

