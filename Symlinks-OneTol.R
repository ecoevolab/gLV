
#--------------------------------Create symbolic links--------------------------------#
# Define source and target directories
source_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/D13M02Y25"

target_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/D13M02Y25-outs"

# Create target directories if they don't exist
if (!dir.exists(target_path)) dir.create(target_path, recursive = TRUE) 

# Get all worker directories
worker_dirs <- list.dirs(source_path, recursive = FALSE, full.names = TRUE)

lapply(worker_dirs, function(worker) {
  
  cat("Starting symbolic link of worker ", worker, " ...\n")
  
  # Get all files in the worker directory recursively
  files <- list.files(worker, recursive = TRUE, full.names = TRUE)
  
  relative_paths <- sub(paste0("^", source_path, "/worker_\\d+/"), "", files)
  
  # Construct target paths
  target_paths <- file.path(target_path, relative_paths)
  
  # Get unique target directories
  target_dirs <- unique(dirname(target_paths))
  
  # Create symbolic links in a vectorized way
  mapply(file.symlink, files, target_paths)
  cat("Worker ", worker, " completed ...\n")
})

 



                 