#' Generate symbolic links of simulation files
#'
#' This function creates symbolic links to the simulation files located in each worker directory and consolidates them into a single target directory. Each worker corresponds to a core utilized during the parallelization of simulations.
#'
#' @param source_path Character string: The path to the root directory containing all worker directories where the simulation files are stored.
#' @param target_path Character string: The path where the symbolic links will be created, i.e., the target directory to unify the simulation files.
#'
#' @return This function does not return a value. It performs the creation of symbolic links in the target directory.
#' 
#' @details 
#' The function iterates through each worker directory located in `source_path`, identifies the simulation files, and generates symbolic links to those files in the `target_path`. This process helps in organizing simulation outputs from different parallel workers into a single location for easy access and further processing.
#'
#' @examples
#' # Example of usage
#' generate_symlinks_OneTol("path/to/worker/directories", "path/to/target/directory")
#'
#' @export
generate_symlinks <- function(source_path, target_path){
  # Create target directories if they don't exist
  if (!dir.exists(target_path)) dir.create(target_path, recursive = TRUE) 
  
  # Get all worker directories
  worker_dirs <- dir(source_path, recursive = FALSE, full.names = TRUE,  pattern = "^worker_.*")
  
  lapply(worker_dirs, function(worker) {
    
    message("Starting symbolic link of ", basename(worker), " ...\n")
    
    # Get all files in the worker directory recursively
    files <- list.files(worker, recursive = TRUE, full.names = TRUE)
    
    relative_paths <- sub(paste0("^", source_path, "/worker_\\d+/"), "", files)
    
    # Construct target paths
    target_paths <- file.path(target_path, relative_paths)
    
    # Get unique target directories and create them if missing
    dirs_to_create <- unique(dirname(target_paths))
    sapply(dirs_to_create, function(d) {
      if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    })
    
    # Create symbolic links, checking if they don't already exist
    invisible(mapply(function(src, tgt) {
      if (!file.exists(tgt)) {
          file.symlink(src, tgt)
        cat("Created symlink for: ", basename(tgt), "\n")
      } else {
        cat("Symlink already exists for: ", basename(tgt), "\n")
      }
      return(NULL)
    }, files, target_paths))
    
    message("Finished symbolic link of ", basename(worker), " ...\n")
  })
  
  return()
}





                 