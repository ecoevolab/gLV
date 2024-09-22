#' Rwindow_individual Function
#'
#' This function analyzes simulation results to identify steady states by calculating rolling means and checking when the differences fall below a specified tolerance.
#'
#' @param uniqueID Character: A unique identifier for the simulation run.
#' @param output Matrix: Simulation results where rows represent species and columns represent time steps.
#' @param tolerance Numeric: Tolerance value for determining the steady state.
#' @param wd Character: Working directory where the results will be saved.
#'
#' @return Numeric vector \code{Stable_vec}, where each index corresponds to a species and contains the generation at which that species reaches the steady state.
#'
#' @details The function calculates the rolling mean of the output with a window size of \code{#generations * 0.1}.
#' It then computes the differences between consecutive rolling means and determines the generation where the difference is less than the specified tolerance.
#'
#' @importFrom zoo rollmean
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
#' individual_rolling_window_SS(uniqueID, output, tolerance, wd)
#'
#' @export

individual_rolling_window_SS <- function(uniqueID, output, tolerance, wd) {

  # Ensure the zoo package is available
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("The 'zoo' package is required but not installed.")
  }

  specs <- nrow(output)  # Number of species
  times <- ncol(output)  # Number of generations
  window_size <- max(round(times * 0.1), 1)  # Calculate window size and ensure it's at least 1
  Stable_vec <- numeric(specs) # Preallocate vector

  # Function to calculate stability for one row
  find_stability <- function(row, window_size, tolerance) {
    moving_avg <- zoo::rollmean(row, window_size, fill = NA, align = "left")  # Calculate moving average
    difference <- abs( diff(moving_avg) )
    stable_gen <- which(difference < tolerance)[1]  # Find the first generation where moving average < tolerance
    return(stable_gen + 1)
  }

  # Apply the stability function to each row
  Stable_vec <- sapply(seq_len(specs), function(s) {
    row <- as.numeric(output[s, ])  # Convert the current row to numeric vector
    find_stability(row, window_size, tolerance)  # Store the result
  })

  names(Stable_vec) <- paste0("Specie", seq_len(specs))

  # Define file path for RDS
  RDS_path <- file.path(wd, "Scan", "SS_individual.rds")

  # Create the new entry
  new_entry <- list(
    ID = uniqueID,
    Method = "Rolling_window",
    Steady_State = Stable_vec
  )

  # Check if the file exists and load/append data accordingly
  if (file.exists(RDS_path)) {
    existing_data <- readRDS(RDS_path)  # Load existing data
    updated_data <- c(existing_data, list(new_entry))  # Append new entry
  } else {
    updated_data <- list(new_entry)  # Initialize new data list if no file exists
  }

  saveRDS(updated_data, file = RDS_path)  # Save updated data

  # Print messages
  cat("Steady State search done and saved\n",
      "\tWith ID:", uniqueID, "\n",
      "\tStable generations path:", RDS_path, "\n\n")

  return(Stable_vec)
}
