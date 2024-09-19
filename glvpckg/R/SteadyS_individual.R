#' SS_diffmeansALL Function
#'
#' This function searches for steady states in a simulation by analyzing the squared differences in the simulation results over time.
#'
#' @param uniqueID Character: A unique identifier for the simulation run.
#' @param output Matrix: Simulation results where rows represent species and columns represent time steps.
#' @param tolerance Numeric: Tolerance value for determining the steady state. It is transformed using \code{log(tolerance^2)}.
#' @param params List: A list containing parameters for the simulation, including 'Interactions' and 'Growths'.
#' @param wd Character: Working directory where the results will be saved.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Stable_gens}: Numeric vector where each index corresponds to a species and contains the generation at which that species reaches the steady state.
#'   \item \code{Transformed_tolerance}: Numeric value of the transformed tolerance, used for comparison. log(tolerance^2)
#'   \item \code{Tolerance}: Numeric value of the original tolerance.
#' }
#'
#' @details The function calculates the squared differences between successive time steps for each species, applies a log transformation,
#' and determines if the differences fall below a specified tolerance. It then calculates the generation at which each species reaches a steady state.
#'
#'
#' @examples
#' # Example usage:
#'
#'wd <- "~/Documents/LAB_ECO"
#'
#' # Initial parameters
#' N_species = 2
#' C0 = 0.45
#' CN = 0.2
#' Diag_val = -0.5
#'
#' # Generate ID
#' uniqueID <- generate_uniqueID(wd)
#'
#' # Generate parameters
#' seeds_path <- file.path(wd, "Seeds.tsv" )
#' params <- generate (N_species, seeds_path, C0, CN, Diag_val)
#'
#' # Generate simulation
#' times <- 20 # Define the number of generations
#' output <- Simulate_output(N_species, params = params, times = times, norm = FALSE)
#'
#' tolerance = 0.05
#' SS_individual(uniqueID, output, tolerance, params, wd)
#'
#' @export

SS_individual <- function(uniqueID, output, tolerance, params, wd) {

  #--------------------Calculate the analytic SS----------------#
  A_inv <- solve(as.matrix(params$Interactions)) # Inverse of the matrix
  detA <- det(params$Interactions)

  if (detA != 0) {
    P <- -A_inv %*% params$Growths
    cat("There is an analytical solution for the Steady State...\n")
  } else {
    P <- numeric(nrow(output)) # No analytical solution
  }

  specs <- nrow(output) # Number of species
  Stable_vec <- numeric(specs) # Preallocate vector
  trans_tolerance <- log(tolerance^2) # Apply transformation to tolerance
  difference <- 0 # Initialize difference

  #----------------------------Search Stability-----------------------#
  for (s in seq_len(specs)) {
    V_spec <- diff(as.numeric(output[s,]))^2 # Get differences and square them
    V_spec <- log(ifelse(V_spec == 0, 10, V_spec)) # Log transformation

    Stable_col <- which(V_spec < trans_tolerance)[1] # Generation where the mean < tolerance

    if (!is.null(P) && length(P) == specs) {
      difference <- difference + (output[s, Stable_col] - P[s])
    }

    Stable_vec[s] <- Stable_col + 1 # Add vector
  }

  if (any(P < 0)) {
    difference <- "Negative numbers found on Populations vector (analytical)"
  }

  names(Stable_vec) <- paste0("Specie", seq_len(specs))

  #------------------------------Create data frame-----------------------------#
  SS_df <- data.frame(
    ID = uniqueID,
    Generations_num = ncol(output),
    Species_num = specs,
    Tolerance = tolerance,
    Transformed_tolerance = format(trans_tolerance, digits = 5, nsmall = 5),
    Expected_vs_Searched = difference,
    Individual = TRUE
  )

  #----------------------------Save Data frame--------------------------------#
  Sind_path <- file.path(wd, "Scan", "SS_individual.tsv")

  if (file.exists(Sind_path)) { # File exists
    SS_table <- read.delim(Sind_path, sep = "\t", header = TRUE) # Read existing table
    SS_df <- rbind(SS_table, SS_df) # Combine with new data
  }

  write.table(SS_df, file = Sind_path, sep = "\t", row.names = FALSE, col.names = TRUE) # Save

  #--------------------------Save stable generations-------------------------#
  # Define file path for RDS
  RDS_path <- file.path(wd, "Scan", "SS_individual.rds")

  # Create the new entry
  new_entry <- list("ID" = uniqueID,
                    "Method" = "log(diff^2)",
                    "Steady_State" = Stable_vec)  # Replace with your steady state vector

  # Check if the file exists and load/append data accordingly
  if (file.exists(RDS_path)) {
    existing_data <- readRDS(RDS_path) # Load existing data
    updated_data <- append(existing_data, new_entry) # Append new entry
  } else {
    updated_data <- list(new_entry) # Initialize new data list if no file exists
  }

  saveRDS(updated_data, file = RDS_path) # Save updated data

  #---------------------Print messages--------------------------------------#
  cat("Steady State search done and saved\n", "\t With ID:", uniqueID, "\n")
  cat ("\t Data Frame path:", Sind_path, "\n")
  cat ("\t Stable generations path:", RDS_path, "\n \n")

  return(list("Stable gens" = Stable_vec,
              "Transformed_tolerance" = trans_tolerance,
              "Tolerance" = tolerance)
  )
}



