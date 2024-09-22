#' individual_diff_SS Function
#'
#' This function searches for steady states in a simulation by analyzing the squared differences in the simulation results over time.
#'
#' @param uniqueID Character: A unique identifier for the simulation run, typically generated using \code{forge_id}.
#' @param output Matrix: Simulation results where rows represent species and columns represent time steps. Each element in the matrix represents the abundance of a species at a given time step.
#' @param tolerance Numeric: Tolerance value for determining the steady state. The tolerance is transformed using \code{tolerance^2} to define the threshold for stability.
#' @param wd Character: Working directory where the results will be saved. The directory should exist prior to running the function.
#'
#' @return A message indicating where the table was saved. The stable generations are saved using the specified unique ID.
#'
#' @details The function calculates the squared differences between successive time steps for each species. If the squared differences fall below the specified tolerance (in its log-transformed form), the species is considered to have reached a steady state. The function then determines the generation at which each species reaches this steady state and saves the results to a file.
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


individual_diff_SS <- function(uniqueID, output, tolerance, wd) {
  
  specs <- nrow(output)  # Number of species
  stable_gen <- numeric(specs)  # Preallocate for stable generations
  ss_counts <- numeric(specs)  # Preallocate for counting steady state generations
  trans_tolerance <- tolerance^2  # Transform tolerance
  
  #----------------------------Search Stability-----------------------#
  for (s in seq_len(specs)) {
    
    V_spec <- diff(as.numeric(output[s, ]))^2  # Get squared differences
    stable_points <- which(V_spec < trans_tolerance)  # Stable generations
    stable_gen[s] <- stable_points[1] # First generation where value < tolerance
    
    if (length(stable_points) > 1) {  # Ensure there are at least 2 points to check
      # Check if the stable points are sequential
      runs <- rle(diff(stable_points) == 1)
      
      # Check if the last run is TRUE (stable state at end)
      if (tail(runs$values, 1)) {
        index <- sum(runs$lengths[-length(runs$lengths)])  # Sum lengths except the last
        cons_gens <- stable_points[(index + 1):length(stable_points)]  # Use vector indexing directly
        ss_counts[s] <- length(cons_gens)
      } else {
        cat("The Steady state is not at the end of the simulation.\nThis could indicate there is not a steady state or the function is oscillatory.")
      }
    } else {
      cat("Not enough stable points to determine sequential generations.\n")
    }
}
  
  # Generation where the steady state is reached
  names(stable_gen) <- paste0("Specie", seq_len(specs))
  
  # Number of Steady Generations of species
  names(ss_counts) <- paste0("#Steady_gens", seq_len(specs))
  
  
  #------------------------------Create Data Frame-----------------------------#
  SS_df <- data.frame(
    ID = uniqueID,
    "#Generations" = ncol(output),
    "#Species" = specs,
    Tolerance = tolerance,
    Transformed_tolerance = format(trans_tolerance, digits = 5, nsmall = 5),
    "#Steady_start" = sum(stable_gen, na.rm = TRUE),
    "#Steady_generations" = sum(ss_counts, na.rm = TRUE),  # Avoid NA in summation
    Method = "diff^2",
    Individual = TRUE
  )
  tmp <- SS_df
  
  #----------------------------Save Data Frame--------------------------------#
  SS_ind_path <- file.path(wd, "Scan", "SS_individual.tsv")
  
  # Read existing table if it exists, else create new
  if (file.exists(SS_ind_path)) {
    SS_table <- read.delim(SS_ind_path, sep = "\t", header = TRUE)  # Read existing data
    SS_df <- rbind(SS_table, SS_df)  # Combine with new data
  }
  
  write.table(SS_df, file = SS_ind_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
  
  #--------------------------Save Stable Generations-------------------------#
  RDS_path <- file.path(wd, "Scan", "SS_individual.rds")
  
  # Create new entry
  new_entry <- list("ID" = uniqueID,
                    "Method" = "log(diff^2)",
                    "Steady_State" = stable_gen)  # Replace with your steady state vector
  
  # Append new entry or initialize new data
  if (file.exists(RDS_path)) {
    existing_data <- readRDS(RDS_path)  # Load existing data
    updated_data <- append(existing_data, new_entry)  # Append new entry
  } else {
    updated_data <- list(new_entry)  # Initialize new data list
  }
  
  saveRDS(updated_data, file = RDS_path)  # Save updated data
  
  #---------------------Print Messages--------------------------------------#
  cat("Steady State search done and saved\n", "\t With ID:", uniqueID, "\n")
  cat("\t Data Frame path:", SS_ind_path, "\n")
  cat("\t Stable generations path:", RDS_path, "\n \n")
  
  return(tmp)
}




