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
gen_syml <- function(source_dir, tgt_dir) {

  cat(rep("-", 30), "\n")
  cat("Starting symbolic link of ", source_dir, " ...\n") 
  
  # tictoc::tic("Starting symbolic links...\n")
  # Define file patterns 
  to_link <- list(
    list(pattern = "O_*.feather", target = file.path(tgt_dir, "raw-ODEs")),                 # Output files
    list(pattern = "E_*-S*.feather", target = file.path(tgt_dir, "exts-outs")),             # extinctions-output
    list(pattern =  "E_*-Info.feather", target = file.path(tgt_dir,  "exts-info"))          # extinctions-info
  )
  dir = to_link[[2]]
  for (dir in to_link) {
    if (!(dir.exists(dir$target))){
      dir.create(dir$target) 
      cat("Directory", dir$target,  "created: \n", sep = " ")
    }
    
    #====================== Source-files ======================
    mcmd <- sprintf('find "%s" -type f -name "%s" ', source_dir, dir$pattern)               # Command
    src_files = system(mcmd, intern = TRUE)
    src_files = list.files(source_dir, full.names = TRUE, recursive = TRUE, pattern = "^E_.*-S.*\\.feather$")
    src_basenames <- basename(src_files)                                                    # source-file-names
    
    #====================== Already-existing-files ======================
    #mcmd <- sprintf('ls "%s" ', dir$target)                    # Command
    #tgt_files = system(mcmd, intern = TRUE)                       # files with symlink already
    #tgt_basenames <- basename(tgt_files)                          # target-file-names

    # Get files missing symlink
    #missing_files <- src_files[!src_basenames %in% tgt_basenames]
    target_paths = file.path(dir$target, src_basenames)
    file.symlink(src_files, target_paths)
    commands <- sprintf('ln -sf "%s" "%s"', src_files, target_paths)
    command_string <- paste(commands, collapse = "\n")
    system("bash", input = command_string)


  }
  # tictoc::toc()
}
