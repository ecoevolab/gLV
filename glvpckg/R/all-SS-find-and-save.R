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
#'  @import utils
#'
#' @examples
#' # Example usage:
#' wd = "~/Documents/LAB_ECO/Simulations"
#'
#' # Generate parameters
#' seeds_path <- file.path(wd, "Seeds.tsv" )
#' params <-  forge_data(N_species = 2, seeds_path, C0 = 0.45, CN =  0.2, Diag_val = -0.5)
#'
#' # Simulate species interactions
#' output <- run_simulation(N_species = 2, params = params, times = 100)
#'
#' # Generate ID
#' uniqueID <- forge_id(wd)
#'
#' # All  Steady states
#' tolerance = 0.05
#' all_SS_find_and_save(uniqueID, output, tolerance, wd)
#' 
#' @export

all_SS_find_and_save <- function(uniqueID, output, tolerance, wd) {
  
  # Attach package
  requireNamespace("utils")

  # Apply functions for SS searching
  result1 <- all_rolling_var_SS(output, tolerance)
  result2 <- all_prop_SS(output, tolerance)

  # Create data frame
  SS_df <- data.frame(
    ID = uniqueID,
    Generations_num = ncol(output),
    Species_num = nrow(output),
    Tolerance = tolerance,
    Roll_var_Stable_generation = result1,
    Diff_means_Stable_generation = result2,
    Individual = FALSE
  )

  # Define the file path
  SS_path <- file.path(wd, "Scan", "SS_ALL.tsv")

  # Save data frame
  if (file.exists(SS_path)) {
    # Read existing data
    SS_table <- utils::read.table(SS_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    # Combine new data with existing data
    SS_join <- rbind(SS_table, SS_df)
    # Save combined data
    utils::write.table(SS_join, file = SS_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  } else {
    # Save new data
    utils::write.table(SS_df, file = SS_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  }

  cat("Steady State search done and saved \n", 
      "\tWith ID", uniqueID, "\n",
      "\tAt path:", SS_path, "\n")
}

