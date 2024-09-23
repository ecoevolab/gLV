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
  times <- ncol(output) # Number of generations
  window_size <- max(1, round(times * 0.1))  # Calculate window size ensuring it's at least 1
  stable_gen <- numeric(specs)  # Vector to store the generations where species reached steady state
  ss_counts <- numeric(specs) # Vector to store the number of generations where species are on steady state
  
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
    
    if (length(stable_points) > 1) {  # Ensure there are at least 2 points to check
      runs <- rle(diff(stable_points) == 1)  # Run-length encoding for sequential checks
      
      if (tail(runs$values)) {  # Check if the last run is TRUE (steady state at end)
        index <- sum(runs$lengths[-length(runs$lengths)])  # Sum lengths except the last
        cons_gens <- stable_points[(index + 1):length(stable_points)]  # Use vector indexing directly
        ss_counts[s] <- length(cons_gens)  # Calculate for how many generations the system is in steady state
      } else {
        message("The steady state is not at the end of the simulation. This could indicate no steady state or oscillatory behavior.")
      }
      
      stable_gen[s] <- ifelse(length(stable_points) > 0, stable_points[1], NA)  # First generation where value < tolerance
    } else {
      message("Not enough stable points to determine sequential generations.")
    }
  }
  
  names(stable_gen) <- paste0("Species", seq_len(specs))
  
  #------------------------------Create Data Frame-----------------------------#
  SS_df <- data.frame(
    ID = uniqueID,
    "#Generations" = times,
    "#Species" = specs,
    Tolerance = tolerance,
    "#Steady_start" = sum(stable_gen, na.rm = TRUE),
    "#Steady_generations" = sum(ss_counts, na.rm = TRUE),  # Avoid NA in summation
    Method = "roll_var",
    Individual = TRUE
  )
  
  tmp <- SS_df
  
  #----------------------------Save Data Frame--------------------------------#
  SS_ind_path <- file.path(wd, "Scan", "SS_Individual_RollingVariance.tsv")
  
  # Read existing table if it exists, else create new
  if (file.exists(SS_ind_path)) {
    SS_table <- read.delim(SS_ind_path, sep = "\t", header = TRUE)  # Read existing data
    SS_df <- rbind(SS_table, SS_df)  # Combine with new data
  }
  
  write.table(SS_df, file = SS_ind_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
  
  #--------------------- Print Messages --------------------------------------#
  cat("Rolling window steady State search done and saved\n",
          "\tWith ID:", uniqueID, "\n",
          "\tData Frame path:", SS_ind_path, "\n")
  
  return(list(table = tmp, # Table saved
              Stable = stable_gen # Vector with on what generation does the specie reached the steady state
              )
  )
}

