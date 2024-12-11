#!/bin/env bash
#$ -N cleaner
#$ -o ./Logs_errors/$JOB_NAME.log
#$ -e ./Logs_errors/$JOB_NAME.error
#$ -cwd
#$ -l h_rt=02:00:00   # Set a runtime limit of 2 hours
#$ -pe openmp 1
#$ -l h_rss=15G        # Request 2GB of memory per core

# Clean directories
find /mnt/atgc-d3/sur/users/mrivera/Outputs -type f -exec rm {} +
find /mnt/atgc-d3/sur/users/mrivera/Parameters -type f -exec rm {} +
find /mnt/atgc-d3/sur/users/mrivera/Scan -type f -exec rm {} +
find /mnt/atgc-d3/sur/users/mrivera/Logs_errors -type f -exec rm {} +
