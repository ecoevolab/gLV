#' SS_diffmeansALL Function
#'
#' This function searches for steady states in a simulation by analyzing the squared differences in the simulation results over time.
#'
#' @param uniqueID
#' @param output Matrix: Simulation results where rows represent species and columns represent time steps.
#' @param tolerance Numeric: Tolerance value for determining the steady state. It is transformed using \code{log(tolerance^2)}.
#' @param wd Character: Working directory where the results will be saved.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{Diff_mean}: Numeric vector. Column means of the squared differences between successive time steps, transformed using \code{log}.
#'     \item \code{Stable_gen}: Data frame containing the following elements:
#'       \itemize {
#'         \item \code{ID}: Character. Unique identifier for the simulation.
#'         \item \code{Generations_num}: Numeric. Number of generations in the simulation.
#'         \item \code{Species_num}: Numeric. Number of species in the simulation.
#'         \item \code{Tolerance}: Numeric. Original tolerance value used.
#'         \item \code{Transformed_tolerance}: Numeric. Tolerance value transformed using \code{log(tolerance^2)}, displayed with 5 decimal places.
#'         \item \code{Stable_generation}: Numeric or NA. Generation where the mean squared differences first fall below the transformed tolerance. \code{NA} if no such generation is found.
#'         \item \code{Individual}: Logical. Flag indicating if the result is for an individual simulation.
#'       }
#'   }
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
#' # Generate parameters
#' seeds_path <- file.path(wd, "Seeds.tsv" )
#' params <- generate (N_species, seeds_path, C0, CN, Diag_val)
#'
#' # Generate simulation
#' times <- 20 # Define the number of generations
#' output <- Simulate_output(N_species, params = params, times = times, norm = FALSE)
#'
#' tolerance = 0.05
#' SS_diffmeansALL(uniqueID, output, tolerance, wd)
#'
#' @export

SS_diffmeansALL <- function(uniqueID, output, tolerance, wd) {

  #------------------------Get differences------------------------#
  gens <- ncol(output) # Number of generations
  specs <- nrow(output) # Number of species
  og_tol <- tolerance
  tolerance <- log(tolerance^2) # Apply transformation to tolerance

  #--------------- Compute differences and square them-----------#
  Diff_output <- matrix(nrow = specs, ncol = gens - 1) # Preallocate matrix
  for (r in 1:specs) {
    v <- output[r, ]
    Diff_output[r, ] <- diff(v)^2
  }

  #---------------------------Create data frame------------------------#
  # Apply log transformation, replace 0 with 10
  ln_mat <- log(ifelse(Diff_output == 0, 10, Diff_output))

  rownames(ln_mat) <- paste("specie", 1:specs, sep = "") # Assign row names
  colnames(ln_mat) <- seq(1, gens - 1) # Assign column names

  #-------------Calculate column means and stable generation-------------------#
  Diff_mean <- colMeans(Diff_output, na.rm = TRUE) # Column means
  Stable_gen <- which(Diff_mean < tolerance)[1] # Generation where difference < tolerance

  if (is.na(Stable_gen)) {
    Stable_gen <- NA # Handle case where no generation meets the tolerance
  } else {
    Stable_gen <- Stable_gen + 1 # Adjust for 1-based index
  }

  #------------------------------Create data frame-----------------------------#
  SS_df <- data.frame(
    'ID' = uniqueID,
    "Generations_num" = gens,
    "Species_num" = specs,
    'Tolerance' = og_tol,
    "Transformed_tolerance" = format(tolerance, digits = 5, nsmall = 5),
    'Stable_generation' = Stable_gen,
    'Individual' = FALSE
  )

  #-------------------------Save data frame--------------------------------#
  SS_path <- file.path(wd, "Scan", "SS_ALL.tsv") # Use file.path for cross-platform compatibility

  if (file.exists(SS_path)) { # File exists
    SS_table <- read.delim(SS_path, sep = "\t", header = TRUE) # Read table
    SS_join <- rbind(SS_table, SS_df) # Join tables
  }

  write.table(SS_join, file = SS_path, sep = "\t", row.names = FALSE, col.names = TRUE) # Save
  cat("Steady Sate search done and Saved \n", "With ID", uniqueID, " \n at path:", SS_path, "\n")

  return(list("Diff_means" = Diff_mean,
              "Stable_df" = SS_df)
         )
}

#' Rwindow_ALL Function
#'
#' This function searches for steady states in a simulation by analyzing the differences in the simulation results over time.
#'
#' @param output Matrix: Simulation results where rows represent species and columns represent time steps.
#' @param tolerance Numeric: Tolerance value for determining the steady state.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{Tolerance_trans}: Numeric. The tolerance value.
#'     \item \code{Diff_mean}: Numeric vector. Column means of the squared differences between successive time steps, transformed using \code{log}.
#'     \item \code{Stable_gen}: Numeric. The generation (time step) at which the system first reaches the steady state, based on the tolerance threshold.
#'   }
#'
#' @examples
#' # Example usage:
#'
#------------------------------------------
#' wd <- "~/Documents/LAB_ECO"

#' # Initial parameters
#' N_species = 2
#' C0 = 0.45
#' CN = 0.2
#' Diag_val = -0.5

#' # Generate parameters
#' seeds_path <- file.path(wd, "Seeds.tsv" )
#' params <- generate (N_species, seeds_path, C0, CN, Diag_val)

#' # Generate simulation
#' times <- 20 # Define the number of generations
#' output <- Simulate_output(N_species, params = params, times = times, norm = FALSE)
#'
#' tolerance = 0.05
#' resul = Rwindow_ALL(output, tolerance)
#' @export

Rwindow_ALL <- function(output, tolerance) {

  # Ensure the zoo package is available
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("The 'zoo' package is required but not installed.")
  }

  #----------------------------Calculate Column Means-----------------#
  times <- ncol(output)  # Number of generations
  col_means <- colMeans(output, na.rm = TRUE)  # Compute column means

  #----------------------------Apply Moving Average---------------------#
  window_size <- round(times * 0.1)  # Calculate window size

  # Function to calculate moving average
  # "left" covers following rows
  moving_avg <- rollmean(col_means, window_size, fill = NA, align = "left")

  # Find stability based on moving average
  Stable_gen <- which(moving_avg < tolerance)[1]  # Find the first generation where moving average < tolerance

  return(list("Tolerance" = tolerance,
              "Stable_gen" = Stable_gen,
              "ColMeans" = col_means)
        )
}
