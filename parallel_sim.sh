#!/bin/env bash

# Define output directories for logs and errors
LOG_DIR="/mnt/atgc-d3/sur/users/mrivera/Logs_errors"

# Define the job name prefix and the R script to run
JOB_PREFIX="Simulations"

# Iterate through different job submissions (varying the job name or internal parameters in R)
for i in {1..10}; do

    # Submit the job with qsub
    qsub -N "${JOB_PREFIX}_${i}" \
         -o "${LOG_DIR}/${JOB_PREFIX}_${i}.log" \
         -e "${LOG_DIR}/${JOB_PREFIX}_${i}.error" \
         -cwd \
         -l h_rt=15:00:00 \
         -l h_rss=30G <<EOF
#!/bin/env bash
# Print the start time
echo "Job start time: $(date)"

# Load the R module (if necessary)
module load r/4.4.0

# Run the R script without external parameters
Rscript Cluster_code.R

# Print the end time
echo "Job end time: $(date)"
EOF

done

