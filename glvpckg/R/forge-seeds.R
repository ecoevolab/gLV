#' Generate seeds for parameter generation
#'
#' This function generates `n` random seeds within a specified range, subsamples them, and selects one seed to generate each item for the \link{forge_data} function.
#'
#' @param n Number of possible seeds to generate.
#' @param min Minimum possible seed value (inclusive).
#' @param max Maximum possible seed value (inclusive).
#' @param wd Working directory where the seeds file will be saved.
#'  
#' @import utils
#' @importFrom random randomNumbers
#' 
#' @return A data frame containing the generated seeds and their corresponding sampled values saved as a file in the specified working directory.
#'
#' @examples
#' # Example usage
#' wd = "~/Documents/LAB_ECO/Simulations"
#' forge_seeds(n = 200, min = 2, max = 1000, wd)
#' # Check the specified working directory for the seeds file
#' 
#' @export


forge_seeds <- function(n, min, max, wd) {
    
    # Attach libraries
    requireNamespace("utils")
    if (!requireNamespace("random", quietly = TRUE)) {
      stop("The 'random' package is required but not installed.")
    }
    
  
    # Validate input
    if (min >= max) {
      stop("'min' must be less than 'max'.")
    }
    
    if (n <= 0) {
      stop("'n' must be a positive number.")
    }
    
    # Generate seeds
    seeds <- random::randomNumbers(n, min, max)
    
    # Check for existing file and create a unique name with 2-digit numbering
    i <- 1
    seeds_file_name <- sprintf("Seeds_%02d.tsv", i)  # Start with Seeds_01.tsv
    seeds_path <- file.path(wd, seeds_file_name)
    
    # Increment the number if the file already exists
    while (file.exists(seeds_path)) {
      i <- i + 1
      seeds_file_name <- sprintf("Seeds_%02d.tsv", i)
      seeds_path <- file.path(wd, seeds_file_name)
    }
    
    # Create file and save seeds
    file.create(seeds_path)
    write.table(seeds, file = seeds_path, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    cat("Seeds generated and saved on path: ", seeds_path)
}

