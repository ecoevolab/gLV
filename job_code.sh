#!/bin/env bash
#$ -N Sim7
#$ -o ../Logs_errors/$JOB_NAME.log
#$ -e ../Logs_errors/$JOB_NAME.error
#$ -cwd
#$ -l h_rt=15:00:00   # Set a runtime limit of 15 hours
#$ -pe openmp 1
#$ -l h_rss=30G       # Request 15GB of memory per core

# Print the start time
# echo "Job start time: $(date)"

# Load the R module (if necessary)
module load r/4.4.0

# Run the R script
Rscript runner.R

# Print the end time
# echo "Job end time: $(date)"

