#' Generate required directories for the package's functions
#'
#' This function creates the necessary directories for simulations in the specified working directory.
#'
#' @param wd Character. The path to the working directory where the directories will be created.
#'
#' @details
#' The function creates directories used for simulations, including: 
#'
#' \itemize{
#'   \item \code{\strong{Outputs} } Directory where the simulation matrix of population changes over time is saved.
#'   \item \code{\strong{Parameters} }: Directory for saving the parameters used to generate the data. Parameters can be saved either by line using \link{params_line_saver}, or by seed using \link{params_seed_saver}.
#'   \item \code{ \strong{Scan} }: Directory where the results of the steady state search algorithms are saved, including global stable generations (using \link{SS_find_and_save_all}) or by individual species (using \link{individual_SS_find_and_save}).
#'   \item  \code{ \strong{Differences} }: Directory where the log-transformed differences \eqn{\log(t+1) - \log(t)} of the output matrix are saved.
#' }
#' 
#' @return Prints messages to the console indicating which directories were created or if they already existed.
#' 
#' @examples
#' # Example usage
#' wd <- "~/Documents/LAB_ECO/Simulations"
#' forge_directories(wd)
#' 
#' @export



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
  
  # Check if parent directory exist
  detect(wd)
  
  # Create Outputs directory
  out_path <- file.path(wd, "Outputs")
  detect(out_path)
  
  # Create Parameters directory
  params_path <- file.path(wd, "Parameters")
  detect(params_path)
  
  # Create Scan directory
  scan_path <- file.path(wd, "Scan")
  detect(scan_path)
  
  # Create Differences directory
  diff_path <- file.path(wd, "Differences")
  detect(diff_path)
}


