#' individual_SS_find_and_save Function
#'
#' This function searches for steady states in a simulation by analyzing both methods available, using rolling variance with a window size of 10% and differences^2.
#'
#' @param uniqueID Character: A unique identifier for the simulation run, typically generated using \code{forge_id}.
#' @param output Matrix: Simulation results where rows represent species and columns represent time steps. Each element in the matrix represents the abundance of a species at a given time step.
#' @param tolerance Numeric: Tolerance value for determining the steady state. The tolerance is transformed on \code{individual_diff_SS} using \code{tolerance^2} to define the threshold for stability.
#' @param wd Character: Working directory where the results will be saved. The directory should exist prior to running the function.
#'
#' @return A message indicating where the table was saved. The stable generations are saved using the specified unique ID.
#'
#' @details The function calculates the start of the steady state using both method available.
#'
#' @examples
#' # Example usage:
#'
#' wd = "~/Documents/LAB_ECO"
#' seeds_path <- file.path(wd, "Seeds.tsv")
#' params <- init_data(N_species = 2, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#' 
#' # Run simulation
#' times <- 100  # Define the number of generations
#' output <- run_simulation(N_species = 2, params = params, times = times)
#' 
#' # Generate unique ID
#' uniqueID <- forge_id(wd)
#'
#' # Search for Individual steady state.
#' tolerance <- 0.00005
#' individual_diff_SS(uniqueID, output, tolerance, wd)
#'
#' @export

individual_SS_find_and_save <- function(uniqueID, output, tolerance, wd) {
  
  
  # Apply functions for SS searching
  result1 <- individual_diff_SS(uniqueID, output, tolerance, wd)
  result2 <- individual_rolling_window_SS(uniqueID, output, tolerance, wd)
  
  #-------------------------- Save Stable Generations -------------------------#
  RDS_path <- file.path(wd, "Scan", "SS_Individual.rds")
  
  # Create the new entry
  new_entry <- list(
    ID = uniqueID,
    Tolerance = tolerance,
    "Differences_method"= result1$Stable,
    "Rolling_variance_method" = result2$Stable
  )
  
  # Append new entry or initialize new data
  if (file.exists(RDS_path)) {
    existing_data <- readRDS(RDS_path)  # Load existing data
    updated_data <- c(existing_data, list(new_entry) )  # Append new entry
  } else {
    updated_data <- list(new_entry)  # Initialize new data list
  }
  
  saveRDS(updated_data, file = RDS_path)  # Save updated data
  
  #-----------------------------------------------------------------------#
  cat("Individual steady state generations search done and saved\n", 
          "\tWith ID:", uniqueID, "\n",
          "\tAt path:", RDS_path, "\n")
  
}







