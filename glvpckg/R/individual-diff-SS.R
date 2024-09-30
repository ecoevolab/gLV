#' Individual Steady State Search by Squared Differences
#'
#' This function searches for steady states in a simulation by analyzing the squared differences in the simulation results over time.
#'
#' @param uniqueID Character: A unique identifier for the simulation run, typically generated using \code{forge_id}.
#' @param output Matrix: Simulation results where rows represent species and columns represent time steps. Each element in the matrix represents the abundance of a species at a given time step.
#' @param tolerance Numeric: Tolerance value for determining the steady state. The tolerance is transformed using \code{tolerance^2} to define the threshold for stability.
#' @param wd Character: Working directory where the results will be saved. The directory should exist prior to running the function.
#'
#' @return Vector with the generation where the species reached the steady state. 
#'
#' @details 
#' The function calculates the squared differences between successive time steps for each species. If the squared differences fall below the transformed tolerance (\code{tolerance^2}), the species is considered to have reached a steady state.
#' 
#' Specifically, the squared difference for each species at time step `t` is computed as:
#' \deqn{D(X_t, X_{t+1}) = (X_{t+1} - X_t)^2}
#' where `D` is the squared difference and `X_t` and `X_{t+1}` are the abundances of the species at consecutive time steps.
#'
#' @examples
#' # Example usage:
#'
#' wd = "~/Documents/LAB_ECO/Simulations"
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

individual_diff_SS <- function(uniqueID, output, tolerance, wd) {
  
  specs <- nrow(output)  # Number of species
  stable_gen <- numeric(specs)  # Preallocate for stable generations
  ss_counts <- numeric(specs)  # Preallocate for counting steady state generations
  trans_tolerance <- tolerance^2  # Transform tolerance
  
  #---------------------------- Search Stability -----------------------#
  for (s in seq_len(specs)) {
    V_spec <- diff(as.numeric(output[s, ]))^2  # Get squared differences
    stable_points <- which(V_spec < trans_tolerance)  # Stable generations
    
    stable_gen[s] <- ifelse(length(stable_points) > 0, stable_points[1], NA)  # First generation where value < tolerance
  }
  
  # Naming stable generations and counts
  names(stable_gen) <- paste0("Specie", seq_len(specs))

  #--------------------- Print Messages --------------------------------------#
  # cat("Individual differences^2 steady State search done \n")
  
  return(stable_gen)  # Vector with the generation where the species reached the steady state
}


