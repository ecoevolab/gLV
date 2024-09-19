#' SS_find_and_save_all Function
#'
#' This function searches for steady states in a simulation by analyzing the differences in the simulation results over time and saves the
#' generation where the steady state is reached.
#'
#' @param uniqueID Character: A unique identifier for the simulation run.
#' @param output Matrix: Simulation results where rows represent species and columns represent time steps.
#' @param tolerance Numeric: Tolerance value for determining the steady state.
#' @param wd Character: Working directory where the results will be saved.
#'
#'
#'
#' @examples
#' # Example usage:
#'
#' wd <- "~/Documents/LAB_ECO" # Working directory
#'
#' # Generate parameters
#' N_species = 2
#' C0 = 0.45
#' CN = 0.2
#' Diag_val = -0.5
#'
#' seeds_path <- file.path(wd, "Seeds.tsv" )
#' params <-  init_data(N_species, seeds_path, C0, CN, Diag_val)
#'
#' # Simulate species interactions
#' times <- 100 # Define the number of generations
#' output <- run_simulation(N_species, params = params, times = times, norm = FALSE)
#'
#' # Generate ID
#' uniqueID <- forge_ID(wd)
#'
#' # All  Steady states
#' tolerance = 0.05
#' SS_find_and_save_all(uniqueID, output, tolerance, wd)
#' @export

SS_find_and_save_all <- function(uniqueID, output, tolerance, wd) {

  # Apply functions for SS searching
  result1 <- SS_roll_window_all(output, tolerance)
  result2 <- SS_diff_means_all(output, tolerance)

  # Create data frame
  SS_df <- data.frame(
    ID = uniqueID,
    Generations_num = ncol(output),
    Species_num = nrow(output),
    Roll_window_Tolerance = tolerance,
    Diff_means_tolerance = result2$Transformed_Tolerance,
    Roll_window_Stable_generation = result1$Stable_gen,
    Diff_means_Stable_generation = result2$Stable_gen,
    Individual = FALSE,
    stringsAsFactors = FALSE
  )

  # Define the file path
  SS_path <- file.path(wd, "Scan", "SS_ALL.tsv")

  # Save data frame
  if (file.exists(SS_path)) {
    # Read existing data
    SS_table <- read.table(SS_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    # Combine new data with existing data
    SS_join <- rbind(SS_table, SS_df)
    # Save combined data
    write.table(SS_join, file = SS_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  } else {
    # Save new data
    write.table(SS_df, file = SS_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  }

  cat("Steady State search done and saved \n", "With ID", uniqueID, "\n at path:", SS_path, "\n")
}

