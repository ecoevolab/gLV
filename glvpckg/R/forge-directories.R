#' forge_directories function
#'
#' This function creates the necessary directories for simulations in the specified working directory.
#'
#' @param wd Character. The path to the working directory where the directories will be created.
#'
#' @details
#' The function constructs directories used for simulations, including "Outputs", "Parameters", and "Scan".
#' 
#' @return Prints messages to the console indicating which directories were created or if they already exist.
#' 
#' @export
#'
#' @examples
#' # Example usage
#' wd <- "~/Documents/LAB_ECO"
#' forge_directories(wd)

forge_directories <- function(wd) {
  
  # Function to check if directories are generated
  detect <- function(directory) {
    
    # Check if the directory exists
    if (!dir.exists(directory)) {
      # If it doesn't exist, create it
      dir.create(directory)
      message("Directory '", directory, "' created.")
    } else {
      message("Directory '", directory, "' already exists. No action taken.")
    }
  }
  
  # Create Outputs directory
  out_path <- file.path(wd, "Outputs")
  detect(out_path)
  
  # Create Parameters directory
  params_path <- file.path(wd, "Parameters")
  detect(params_path)
  
  # Create Scan directory
  scan_path <- file.path(wd, "Scan")
  detect(scan_path)
}


