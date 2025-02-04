
#--------------------------------Create symbolic links--------------------------------#

# Define source and target directories
source_parent <- "/mnt/atgc-d3/sur/users/mrivera/testing/test02"
# source_parent <- "/home/rivera/Cluster/testing/test02"

target_parent <- "/mnt/atgc-d3/sur/users/mrivera/testing/tmp-test"
# target_parent <- "/home/rivera/Cluster/testing/Tolerances"

# Get all worker directories
worker_dirs <- list.dirs(source_parent, recursive = FALSE, full.names = TRUE)

# Define a function to process each worker directory
process_worker <- function(worker) {
  # Get all files in the worker directory recursively
  files <- list.files(worker, recursive = TRUE, full.names = TRUE)

  # Iterate over each file
  for (file_path in files) {
    # Extract relative path
    relative_path <- sub(paste0("^", source_parent, "/worker_\\d+/"), "", file_path)

    # Construct the target path
    target_path <- file.path(target_parent, relative_path)

    # Ensure the target directory exists
    target_dir <- dirname(target_path)
    if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)

    # Create symbolic link
    file.symlink(file_path, target_path)
  }

  msg <- cat("The symbolic links of worker ", worker, " were completed...\n")
  return(msg)
}

library(parallel)
num_cores <- detectCores() - 1  # Use one less than the total number of cores
cat("The number of cores that will be used are: ", num_cores, "\n")

# Use mclapply for parallel processing (for Unix-like systems)
mclapply(worker_dirs, process_worker, mc.cores = num_cores)
