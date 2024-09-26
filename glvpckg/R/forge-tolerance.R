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
#' wd <- "~/Documents/LAB_ECO/testing"
#' seeds_path <- file.path("~/Documents/LAB_ECO", "Seeds.tsv")
#' 
#' paramSettings <- list(N_species = 2, 
#'                       C0 = 0.45, 
#'                       CN = 0.2, 
#'                       Diag_val = -0.5, # Parameters for data generation
#'                       tolerance = 0.05, # Tolerance used for Steady State search
#'                       times = 100) # Number of generations   
#' individual = TRUE
#' 
#' # Function
#' forge_tolerance(paramSettings, seeds_path, wd, individual)
#' @export

forge_tolerance <- function(paramSettings, seeds_path, wd, individual) {
  
  # Set initial tolerance
  tolerance <- paramSettings$tolerance
  counter <- 1
  
  # Generate unique ID for this run
  uniqueID <- forge_id(wd)
    
  # Generate directories once before the loop
  suppressMessages( forge_directories(wd) )
    
  # Pre-generate the simulation parameters that do not change
  params <- init_data(
    N_species = paramSettings$N_species, 
    seeds_path = seeds_path, 
    C0 = paramSettings$C0, 
    CN = paramSettings$CN, 
    Diag_val = paramSettings$Diag_val
  )
  
  # Save parameters and output data
  params_seed_saver(
      N_species = paramSettings$N_species,
      C0 = paramSettings$C0,
      CN = paramSettings$CN,
      Diag_val = paramSettings$Diag_val,
      params = params,
      uniqueID = uniqueID,
      wd = wd)
  
   # Run simulation and save it
    output <- run_simulation(N_species = paramSettings$N_species, params = params, times = paramSettings$times)
    output_saver(output, uniqueID, wd) 
    
 repeat {
    
    if (individual){
      # Calculate and save individual steady states
      capture.output(individual_SS_find_and_save(uniqueID, output, tolerance, wd) )
      message("Individual steady state generations search done and saved with tolerance ", tolerance)
    } else {
      # Calculate and save all steady states
      capture.output(SS_find_and_save_all(uniqueID, output, tolerance, wd) )
      message("ALL steady state generations search done and saved with tolerance ", tolerance)
    }
    
    # Decrease tolerance by 1 order of magnitude
    tolerance <- tolerance / 10
    counter <- counter + 1
    
    if (counter >= 7) {
      break
    }
  }
  
  message("Steady Search completed...")
  
  return(uniqueID)
}
