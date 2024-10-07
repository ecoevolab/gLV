#' Search for the best tolerance
#'
#' This function generates the necessary parameters for a simulation and compares different tolerance values to identify the best one for finding steady states.
#'
#' @param paramSettings A list containing the following elements:
#'   \itemize{
#'     \item \code{C0}: Numeric. The probability of no interaction between species (i.e., 0 interaction).
#'     \item \code{CN}: Numeric. The probability of negative interaction (<0) between species.
#'     \item \code{Diag_val}: Numeric. The diagonal values used in the interaction matrix.
#'     \item \code{tolerance}: Numeric. Tolerance value for determining the steady state. The tolerance is transformed in \code{individual_diff_SS} using \code{tolerance^2} to define the threshold for stability.
#'     \item \code{times}: Numeric. The number of generations to simulate.
#'   }
#' @param seeds_path Character. The path to the seeds file used for the simulation.
#' @param wd Character. Working directory where the results will be saved.
#' @param individual Logical. If TRUE, the steady state search will be conducted per species rather than globally. 
#'
#' @return A character string containing a unique ID to be used for saving files.
#'
#' @examples
#' # Example usage:
#' 
#' wd = "~/Documents/LAB_ECO/Simulations"
#' # forge_seeds(n = 200, min = 2, max = 1000, wd)
#' seeds_path <- file.path(wd, "Seeds.tsv" )
#' params <- init_data(N_species = 2, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#' 
#' # Run simulation
#' times <- 100  # Define the number of generations
#' output <- run_simulation(N_species = 2, params = params, times = times)
#' 
#' # Generate unique ID
#' uniqueID <- forge_id(wd)
#' 
#'  # Function
#' forge_tolerance(output, uniqueID, initial_tolerance = 0.05, individual = TRUE)
#' 
#' @export

forge_tolerance <- function(output, uniqueID, initial_tolerance, individual) {
  
    # Set initial tolerance
    tolerance <- initial_tolerance
    counter <- 1
    
    repeat {
      
      if (individual){
        # Calculate and save individual steady states
        individual_SS_find_and_save(uniqueID, output, tolerance, wd) 
        cat("Individual steady state generations search done and saved with tolerance ", tolerance, "\n")
        
      } else {
        # Calculate and save all steady states
        capture.output(all_SS_find_and_save(uniqueID, output, tolerance, wd) )
        cat("ALL steady state generations search done and saved with tolerance ", tolerance, "\n")
      }
      
      # Decrease tolerance by 1 order of magnitude
      tolerance <- tolerance / 10
      counter <- counter + 1
      
      if (counter >= 7) {
        break
      }
    }
    
    # Print message
    message("Steady Search completed with ID: ", uniqueID)
  
}
