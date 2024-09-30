#' Individual Steady State Search by Proportion of Change
#'
#' This function searches for steady states in a simulation by analyzing the \eqn{\log(t / (t+1))} of species abundance across generations.
#'
#' @param uniqueID Character: A unique identifier for the simulation run, typically generated using \code{forge_id}.
#' @param output Matrix: Simulation results where rows represent species and columns represent time steps. Each element in the matrix represents the abundance of a species at a given time step.
#' @param tolerance Numeric: Tolerance value for determining the steady state. The tolerance represents the proportional change between successive generations, calculated using \eqn{\log(t / (t+1))}.
#' @param wd Character: The working directory where the results will be saved. The directory must exist prior to running the function.
#'
#' @details
#' The function computes the logarithmic differences of species abundance to identify stable generations. 
#' It determines a steady state when the absolute difference in logarithmic values between consecutive generations falls 
#' below the specified tolerance. This calculation leverages the logarithmic property \eqn{\log\left(\frac{x}{y}\right) = \log(x) - \log(y)} to facilitate the analysis.
#' 
#' Values of zero in the \code{output} matrix are treated as \code{NA} to avoid 
#' taking the logarithm of zero. If no stable generations are found for a species, the output will reflect this 
#' with an \code{NA} entry for that species.
#'
#' @return A numeric vector indicating the generation where each species reached the steady state. 
#' If a species did not reach a steady state, the corresponding value will be \code{NA}.
#'
#' @examples
#' wd = "~/Documents/LAB_ECO/Simulations"
#' seeds_path <- file.path(wd, "Seeds.tsv" )
#' 
#' # Generate parameters
#' params <- init_data(N_species = 2, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#' 
#' # Run simulation
#' times <- 100  # Define the number of generations
#' output <- run_simulation(N_species = 2, params = params, times = times)
#' 
#' # Generate unique ID
#' uniqueID <- forge_id(wd)
#' 
#' # Example usage:
#' tolerance=0.05
#' individual_prop_SS(uniqueID, output, tolerance, wd)
#'
#' @export

individual_prop_SS <- function(uniqueID, output, tolerance, wd) {
  
  # Ensure required package is loaded
  requireNamespace("utils")
  
  #--------------------------Declare data-------------------------------#
  specs <- nrow(output)  # Number of species
  gens <- ncol(output)
  stable_gen <- numeric(specs)  # Preallocate for stable generations
  log_diff_df <- matrix(NA, nrow = specs, ncol = gens - 1)  # Preallocate log_diff matrix
  
  #---------------------------- Precompute Log Safely ------------------#
  # Replace 0 values in output with NA to avoid log(0) = -Inf
  safe_output <- ifelse(output == 0, NA, output)
  log_output <- log(safe_output)  # Apply logarithm safely
  
  #---------------------------- Search Stability ------------------------#
  for (s in seq_len(specs)) {
    V_spec <- log_output[s, ]  # Use safe log values
    log_diff <- diff(V_spec)
    log_diff_df[s, ] <- log_diff
    
    # Get stable generations
    stable_points <- which(abs(log_diff) < tolerance & !is.na(log_diff))  # Ignore NAs
    stable_gen[s] <- ifelse(length(stable_points) > 0, stable_points[1], NA)
  }
  
  # Naming stable generations
  names(stable_gen) <- paste0("Specie", seq_len(specs))
  
  #------------------------------- Save differences -----------------------------#
  # Log differences path
  diff_path <- file.path(wd, "Differences", paste0("ld", uniqueID, ".tsv"))
  
  # Convert log_diff_df to data frame for saving
  log_diff_df <- as.data.frame(log_diff_df)
  colnames(log_diff_df) <- paste0("Diff", seq_len(gens - 1))
  
  # Save log differences
  utils::write.table(log_diff_df, file = diff_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  #--------------------- Print Messages -----------------------------------------#
  # cat("Log differences steady state search done and saved\n",
  #     "\tWith ID:", uniqueID, "\n",
  #     "\tLog differences saved at path:", diff_path, "\n")
  
  return(stable_gen)
}
