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
gen_syml <- function(worker_path, exp_dir) {

  cat(rep("-", 30), "\n")
  cat("Starting symbolic link of ", basename(worker_path), " ...\n")
    
  # Define file patterns 
  file_patterns <- c(
    "^O_.*\\.feather$",           # Output files
    "^E_.*-Info\\.feather$",      # Info extension files
    "^E_.*-Sp.*\\.feather$"   # Spec extension files
  )

  # Define target directories
  target_dirs <- c(
    file.path(exp_dir, "Outputs"),
    file.path(exp_dir, paste0("Exts-", "Info")),
    file.path(exp_dir,  paste0("Exts-", "Specs"))
  )


  # Find files and create symlinks in one operation
  mapply(function(pattern, target_dir) {
    # Find files matching the pattern
    files <- list.files(worker_path, recursive = TRUE, full.names = TRUE, pattern = pattern)
    
    # Skip if no files found
    if (length(files) == 0) return(invisible(NULL))
    
    # Create target directory if needed
    if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
    
    # Create target paths
    target_paths <- file.path(target_dir, basename(files))
    
    # Create symlinks for files that don't already have them
    needs_link <- !file.exists(target_paths)
    
    if (any(needs_link)) {
      # Create symlinks only where needed
      mapply(file.symlink, files[needs_link], target_paths[needs_link])
      cat("Created", sum(needs_link), "symlinks for pattern:", pattern, "\n")
    }
    
    # Report existing symlinks
    if (any(!needs_link)) {
      cat("Skipped", sum(!needs_link), "existing symlinks for pattern:", pattern, "\n")
    }
    
    invisible(NULL)
  }, file_patterns, target_dirs, SIMPLIFY = FALSE)

  cat("Ending symbolic link of ", basename(worker_path), " ...\n")
}





                 