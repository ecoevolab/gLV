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
#' @details The function calculates the squared differences between successive time steps for each species. If the squared differences fall below the tranformed tolerance (\code{tolerance^2}), the species is considered to have reached a steady state. 
#' 
#' The function then determines at how many generations does the specie keep on steady state at the simulation.
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
  
  #---------------------------- Search Stability -----------------------#
  for (s in seq_len(specs)) {
    V_spec <- diff(as.numeric(output[s, ]))^2  # Get squared differences
    stable_points <- which(V_spec < trans_tolerance)  # Stable generations
    
    stable_gen[s] <- ifelse(length(stable_points) > 0, stable_points[1], NA)  # First generation where value < tolerance
    
    if (length(stable_points) > 1) {  # Ensure there are at least 2 points to check
      runs <- rle(diff(stable_points) == 1)  # Check for sequential stable points
      
      # Check if the last run indicates a steady state at the end
      if (tail(runs$values, 1)) {
        index <- sum(runs$lengths[-length(runs$lengths)])  # Sum lengths except the last
        ss_counts[s] <- length(stable_points) - index  # Count stable generations
      } else {
        message("The steady state is not at the end of the simulation. This could indicate no steady state or oscillatory behavior.")
      }
    } else {
      message("Not enough stable points to determine sequential generations.")
    }
  }
  
  # Naming stable generations and counts
  names(stable_gen) <- paste0("Specie", seq_len(specs))
  names(ss_counts) <- paste0("#Steady_gens", seq_len(specs))
  
  #------------------------------ Create Data Frame -----------------------------#
  SS_df <- data.frame(
    ID = uniqueID,
    "#Generations" = ncol(output),
    "#Species" = specs,
    Tolerance = tolerance,
    Transformed_tolerance = format(trans_tolerance, digits = 5, nsmall = 5),
    "#Steady_start" = sum(stable_gen, na.rm = TRUE),
    "#Steady_generations" = sum(ss_counts, na.rm = TRUE),  # Avoid NA in summation
    Method = "diff^2",
    Individual = TRUE,
    stringsAsFactors = FALSE
  )
  
  tmp <- SS_df
  
  #---------------------------- Save Data Frame --------------------------------#
  SS_ind_path <- file.path(wd, "Scan", "SS_individual_differences.tsv")
  
  # Read existing table if it exists, else create new
  if (file.exists(SS_ind_path)) {
    SS_table <- read.delim(SS_ind_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)  # Read existing data
    SS_df <- rbind(SS_table, SS_df)  # Combine with new data
  }
  
  write.table(SS_df, file = SS_ind_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
  
  #--------------------- Print Messages --------------------------------------#
  cat("Individual differences steady State search done and saved\n", 
          "\tWith ID:", uniqueID, "\n",
          "\tData Frame path:", SS_ind_path, "\n")
  
  return(list(table = tmp,  # Table saved
              Stable = stable_gen))  # Vector with the generation where the species reached the steady state
}


