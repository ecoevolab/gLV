#' Individual Steady State Search by Rolling Variance
#'
#' This function analyzes simulation results to identify steady states by calculating rolling variance of species abundances across time steps. It checks when the variance fall below a specified tolerance level, indicating that the species has stabilized.
#'
#' @param uniqueID Character: A unique identifier for the simulation run, typically used for tracking and organizing results.
#' @param output Matrix: A matrix of simulation results where each row represents a different species and each column represents a specific time step. The values in the matrix indicate the abundance of each species at each time step.
#' @param tolerance Numeric: A tolerance value that determines the threshold for considering a species to be at steady state. This value represents the maximum allowable difference in rolling means that can occur between consecutive generations for a species to be deemed stable.
#' @param wd Character: The working directory where the results will be saved. This directory must exist prior to running the function.
#'
#' @return A numeric vector \code{Stable_vec}, where each index corresponds to a species and contains the generation at which that species reaches the steady state. If a species does not reach a steady state, the corresponding value will be `NA`.
#'
#' @details 
#' The function calculates the rolling variance of the output using a window size of \eqn{generations * 0.1}. It determines the generation at which the rolling variance falls below the specified tolerance, indicating that the species has stabilized.
#' 
#' The variance for a given window size `k` is calculated using the equation:
#' 
#' \deqn{Var(X) = \frac{1}{n} \sum_{i=1}^{n} (X_i - \bar{X})^2}
#' 
#' where `n` is the size of the rolling window and \bar{X} is the mean of the values within that window.
#' 
#' For example, when analyzing the first row of the output matrix with a window size of 10, the variance is calculated using the values from columns 1 to 10. This rolling process continues for each subsequent column, moving the window forward by one column at a time. At each step, the variance is recalculated by including the new value and excluding the oldest value in the window. 
#' 
#' The process stops when there are not enough values remaining to form a complete window of size `k`, such as before the last `k` columns.
#' 
#' @importFrom zoo rollapply
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
#' tolerance <- 0.005
#' individual_rolling_window_SS(uniqueID, output, tolerance, wd)
#'
#' @export

individual_rolling_variance_SS <- function(uniqueID, output, tolerance, wd) {
  
  # Ensure the zoo package is available
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("The 'zoo' package is required but not installed.")
  }
  
  #--------------------------Declare data-------------------------------#
  specs <- nrow(output)  # Number of species
  times <- ncol(output) # Number of generations
  window_size <- max(1, round(times * 0.1))  # Calculate window size ensuring it's at least 1
  stable_gen <- numeric(specs)  # Vector to store the generations where species reached steady state

  # Function to calculate stability for one row
  find_stability <- function(row) {
    moving_var <- zoo::rollapply(row, width = window_size, FUN = var, fill = NA, align = "left")  # Calculate moving variance
    stable_points <- which(moving_var < tolerance)  # Find generations where moving variance < tolerance
    return(stable_points)
  }
  
  #-------------------------- Calculate Stable Generations -------------------------#
  for (s in seq_len(specs)) {
    
    # Directly convert the current row to numeric
    stable_points <- find_stability(as.numeric(output[s, ]))  
    stable_gen[s] <- ifelse(length(stable_points) > 0, stable_points[1], NA)  # First generation where value < tolerance
  }
  
  names(stable_gen) <- paste0("Species", seq_len(specs))

  #--------------------- Print Messages --------------------------------------#
  # cat("Rolling window variance steady State search done \n")
  
  return(stable_gen)
}

