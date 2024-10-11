#' Individual Steady State Search with All Methods
#'
#' This function searches for steady states in a simulation by analyzing all available methods: rolling variance with a window size of 10%, squared differences between generations, and differences of log abundances between consecutive time steps.
#'
#' @param uniqueID Character: A unique identifier for the simulation run, typically generated using \code{forge_id}.
#' @param output Matrix: Simulation results where rows represent species and columns represent time steps. Each element in the matrix represents the abundance of a species at a given time step.
#' @param tolerance Numeric: Tolerance value for determining the steady state. 
#' @param wd Character: Working directory where the results will be saved. The directory should exist prior to running the function.
#'
#' @return A message indicating where the results were saved (as an RDS file). The stable generations are stored using the specified unique ID.
#'
#' @details 
#' The function calculates the steady state of the output using three methods:
#' \itemize{
#'   \item \link{individual_raw_diff_SS}: Computes steady states by analyzing squared differences between generations \eqn{(t+1) - t} . On this method tolerance is transformed using \code{tolerance^2} to define the threshold for stability.
#'   \item \link{individual_rolling_variance_SS}: Determines steady states using rolling variance with a specified window size.
#'   \item \link{individual_prop_SS}: Analyzes differences in log abundances to identify steady states.
#' }
#'
#' This function saves simulation results in separate RDS files based on the method used. 
#' If a new simulation is conducted with an existing ID, the results will be appended to the corresponding RDS file. Note that the function does not verify whether a 
#' simulation with the same ID and tolerance has already been stored in the file.
#' 
#' @examples
#' # Example usage:
#'
#' wd = "~/Documents/LAB_ECO/Simulations"
#' seeds_path <- file.path(wd, "Seeds.tsv")
#' params <- forge_data (N_species = 2, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#' 
#' # Run simulation
#' output <- run_simulation(N_species = 2, params = params, times = 100)
#' 
#' # Generate unique ID
#' uniqueID <- forge_id(wd)
#'
#' # Search for Individual steady state.
#' tolerance <- 0.005
#' individual_SS_find_and_save(uniqueID, output, tolerance, wd)
#'
#' @export

individual_SS_find_and_save <- function(uniqueID, output, tolerance, wd) {
  
  # Apply functions for SS searching
  result1 <- individual_raw_diff_SS(uniqueID, output, tolerance, wd) 
  result2 <- individual_rolling_variance_SS(uniqueID, output, tolerance, wd) 
  result3 <- individual_prop_SS(uniqueID, output, tolerance, wd) 
  
  #-------------------------- Do list of Stable Generations -------------------------#
 
  
  # Function to save JSON file per method
  entry_per_method <- function(result, tolerance, method_name) {
    
    entry <- list(
      Tolerance = tolerance,
      Generations = result
    )
    
    nested_entry <- list(ID = uniqueID,
                     data = entry)
    
    # Generate file path
    RDS_path <- file.path(wd, "Scan", paste0(method_name, ".rds", sep = "") )
    
    # Append new entry or initialize new data
    if (file.exists(RDS_path)) {
      
        # Load existing data
        old_list <- readRDS(RDS_path)  
        
        # Append new entry
        updated_entry <- list(old_list, nested_entry)
        
        # Save appended entry
        saveRDS(updated_entry, file = RDS_path) 
        #cat("RDS file updated")
      } else {
      
      # file doesn't exist, so we save it 
      saveRDS(nested_entry, file = RDS_path) 
      #cat("RDS file generated, it didnt exist")
    }
  }
  
  entry1 <- entry_per_method(result1, tolerance, method_name = "Raw_diff")
  entry2 <- entry_per_method(result2, tolerance, method_name = "Rolling_var")
  entry3 <- entry_per_method(result3, tolerance, method_name = "Log_diff")

  
  #--------------------------Get maximum and minimum by each method-----------#
  SS_df <- data.frame(
    ID = uniqueID,
    tolerance = tolerance,
    "Max_diff" = max(result1), "#Specie" = which.max(result1),
    "Min_diff" = min(result1), "#Specie"  = which.min(result1),
    "Max_Rolling_var" = max(result2), "#Specie"  = which.max(result2),
    "Min_Rolling_var" = min(result2), "#Specie"  = which.min(result2),
    "Max_Proportions" = max(result3), "#Specie"  = which.max(result3),
    "Min_Proportions" = min(result3), "#Specie"  = which.min(result3)
  )
  
  tmp <- SS_df
  info_path <- file.path(wd, "Scan/SS_diff_table.tsv")
  
  # Read existing table if it exists, else create new
  if (file.exists(info_path)) {
    SS_table <- read.delim(info_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)  # Read existing data
    SS_df <- rbind(SS_table, SS_df)  # Combine with new data
  } else {
     write.table(SS_df, file = info_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
  }
  
  #-----------------------------------------------------------------------#
  # cat("Individual steady state generations search done and saved\n", 
  #         "\tWith ID:", uniqueID, "\n",
  #         "\tSS table saved at path:", info_path, "\n")
  
  return(tmp)
}







