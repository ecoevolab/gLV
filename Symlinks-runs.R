
#--------------------------------Create symbolic links--------------------------------#

# Define source and target directories
source_parent <- "/mnt/atgc-d3/sur/users/mrivera/testing/repetition"

target_parent <- "/mnt/atgc-d3/sur/users/mrivera/testing/Sims_Exp02/Combination-tols"


# Get all worker directories
worker_dirs <- list.dirs(source_parent, recursive = FALSE, full.names = TRUE)

process_worker <- function(worker) {
  
  # Get all files in the worker directory recursively
  files <- list.files(worker, recursive = TRUE, full.names = TRUE)
  
  relative_paths <- sub(paste0("^", source_parent, "/worker_\\d+/"), "", files)
  
  # Construct target paths
  target_paths <- file.path(target_parent, relative_paths)
  
  # Get unique target directories
  target_dirs <- unique(dirname(target_paths))
  
  # Create target directories only if they don't exist
  for (dir in target_dirs) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  }
  
  # Create symbolic links in a vectorized way
  mapply(file.symlink, files, target_paths)
  
  
 cat("The symbolic links of worker ", worker, " were completed...\n")
}

lapply(worker_dirs, process_worker)
