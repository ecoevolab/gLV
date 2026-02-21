
#--------------------------------Create symbolic links--------------------------------#
# Define source and target directories
source_path <- c("/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Experiment-01/D10M02Y24_01",
                 "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Experiment-01/D10M02Y24-RAW")

target_path <- c("/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Unified/PropsUnified-D10M02Y24",
                 "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Unified/RawUnified-D10M02Y24")

process_links <- function(source_path, target_path) {
  
  # Get all worker directories
  worker_dirs <- list.dirs(source_path, recursive = FALSE, full.names = TRUE)
  
  lapply(worker_dirs, function(worker) {
    # Get all files in the worker directory recursively
    files <- list.files(worker, recursive = TRUE, full.names = TRUE)
    
    relative_paths <- sub(paste0("^", source_path, "/worker_\\d+/"), "", files)
    
    # Construct target paths
    target_paths <- file.path(target_path, relative_paths)
    
    # Get unique target directories
    target_dirs <- unique(dirname(target_paths))
    
    # Create target directories if they don't exist
    lapply(target_dirs, function(dir) if (!dir.exists(dir)) dir.create(dir, recursive = TRUE))
    
    # Create symbolic links in a vectorized way
    mapply(file.symlink, files, target_paths)
    # cat("The symbolic links of worker ", worker, " were completed...\n")
    print(paste("Mapping worker:", source_path, "->", target_path))
  })
  
  cat( rep("-", 20), "\nSymbolic links of directory ", source_path, " COMPLETED\n", rep("-", 20))
}


# Apply the function with correct index mapping
Map(process_links, source_path, target_path)

# Map(function(src, tgt) {
#   print(paste("Mapping:", src, "→", tgt))
# }, source_path, target_path)


# Define source and target directories
source_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Experiment-01/D10M02Y24_01"

target_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Unified/PropsUnified-D10M02Y24"
                 