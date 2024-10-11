#' Search for the Best Tolerance
#'
#' This function generates the necessary parameters for a simulation and compares different tolerance values to identify the best one for finding steady states.
#'
#' @param output Matrix. A matrix of simulation results where each row represents a different species and each column represents a specific time step. The values in the matrix indicate the abundance of each species at each time step.
#' @param uniqueID Character. A unique identifier for the simulation run, typically used for tracking and organizing results.
#' @param initial_tolerance Numeric. A tolerance value that determines the threshold for considering a species to be at steady state. 
#' @param individual Logical. If TRUE, the steady state search will be conducted per species rather than globally. 
#' @param wd Character. The path to the working directory where files will be saved. 
#' 
#' @return A character string containing a unique ID to be used for saving files related to the simulation results.
#'
#' @examples
#' # Example usage:
#' 
#' wd = "~/Documents/LAB_ECO/Simulations"
#' seeds_path <- file.path(wd, "Seeds.tsv" )
#' params <- forge_data(N_species = 2, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#' 
#' # Run simulation
#' output <- run_simulation(N_species = 2, params = params, times = 100)
#' 
#' # Generate unique ID
#' uniqueID <- forge_id(wd)
#' 
#'  # Function
#' forge_tolerance(output, uniqueID, initial_tolerance = 0.05, individual = TRUE, wd)
#' 
#' @export

forge_tolerance <- function(output, uniqueID, initial_tolerance, individual, wd) {
  
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
