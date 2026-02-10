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
gen_syml <- function(source_dir) {

  cat(rep("-", 30), "\n")
  cat("Starting symbolic link of ", source_dir, "\n") 
  exp_dir = dirname(source_dir)

  # Define file patterns 
  to_link <- list(
    list(pattern = "RawOutput_*.feather", target = file.path(exp_dir, "Outputs")),            # Output files
    list(pattern = "E_*-S*.feather", target = file.path(exp_dir, "ExtOutputs")),               # extinctions-output
    list(pattern =  "ExtSummary_*.feather", target = file.path(exp_dir,  "ExtSummaries")),    # extinctions-summary
    list(pattern =  "A_*.feather", target = file.path(exp_dir,  "Interactions")),             # interactions matrix      
    list(pattern =  "Topology_*.feather", target = file.path(exp_dir,  "Topologies"))         # Network-topology
  )
  
  for (dir in to_link) {
    # Look for files
    mcmd <- sprintf('find "%s" -type f -name "%s" ', source_dir, dir$pattern)               # Command
    src_files = system(mcmd, intern = TRUE)                                                 # Full-paths
    if (length(src_files)==0) {
      cat("No files found for pattern:", dir$pattern, "\n")
      next
    }
    # Create dirctory if it doesnt exist
    if (!(dir.exists(dir$target))){
      dir.create(dir$target) 
      cat("Directory", dir$target,  "created: \n", sep = " ")
    }
    # Generate-new-paths
    target_paths = file.path(dir$target, basename(src_files))
    file.rename(src_files, target_paths)
    cat(">> Generated", length(target_paths), "for pattern:", dir$pattern, "\n", sep=" ")
  }
}
