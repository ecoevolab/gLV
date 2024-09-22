#' forge_seeds Function
#'
#' This function generates `n` number of random seeds, subsamples them, and selects one to generate each item for the `init_data` function.
#'
#' @param n Number of possible seeds generated.
#' @param min Minimum possible seed value.
#' @param max Maximum possible seed value.
#' @param wd Working directory where the seeds file will be saved.
#'  
#' @import utils
#' @importFrom random randomNumbers
#' 
#' @examples
#' # Example usage
#' wd = "~/Documents/LAB_ECO/"
#' forge_seeds(n = 200, min = 2, max = 1000, wd)
#' 
#' @export
#' 

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
    
    # Generate path to save the seeds
    seeds_path <- file.path(wd, "Seeds.tsv")
    file.create(seeds_path) # Create file
    write.table(seeds, file = seeds_path, sep = "\t", row.names = FALSE, col.names = FALSE) # Save
    
    cat("Seeds generated and saved on path: ", seeds_path)
}

