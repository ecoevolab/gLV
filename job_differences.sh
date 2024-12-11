#!/bin/env bash
#$ -N out_differences
#$ -o ./Logs_errors/$JOB_NAME.log
#$ -e ./Logs_errors/$JOB_NAME.error
#$ -cwd
#$ -l h_rt=15:00:00   # Set a runtime limit of 15 hours
#$ -pe openmp 1
#$ -l h_rss=15G        # Request 15GB of memory per core

# Load the R module (if necessary)
module load r/4.4.0

# Run the R script
Rscript Diff_cluster.R
