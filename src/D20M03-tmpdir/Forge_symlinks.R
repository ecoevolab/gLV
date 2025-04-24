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



generate_symlinks <- function(workers_path, outs_path, exts_path){
  
  # Create target directories if they don't exist
  if (!dir.exists(outs_path)) dir.create(outs_path, recursive = TRUE) 
  if (!dir.exists(exts_path)) dir.create(exts_path, recursive = TRUE) 
  
  # Get all worker directories
  worker_dirs <- dir(workers_path, recursive = FALSE, full.names = TRUE,  pattern = "^worker_.*")
  
  lapply(worker_dirs, function(worker) {
    
    # worker <- worker_dirs[1]
    message("Starting symbolic link of ", basename(worker), " ...\n")
    
    # Get all files in the worker directory recursively
    outputs_files <- list.files(path = worker, pattern = "^O_.*\\.tsv$", full.names = TRUE)
    extinct_files <- list.files(path = worker, pattern = "^Ext_.*\\.rds$", full.names = TRUE)
    
    # Create symlinks in the corresponding directories
    file.symlink(outputs_files, outs_path)
    file.symlink(extinct_files, exts_path)
    
    message("Finished symbolic link of ", basename(worker), " ...\n")
  })
  
  return()
  
}




                 