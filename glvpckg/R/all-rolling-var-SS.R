#' All Steady State Search by Rolling Variance
#'
#' @details
#' This function searches for steady states in a simulation by analyzing the changes in the simulation results over time using rolling variance. The steady state is detected when the variance between successive time steps falls below a given tolerance threshold.
#'
#' @param output Matrix: A matrix of simulation results where rows represent species and columns represent time steps. Typically generated with \link{run_simulation}.
#' @param tolerance Numeric: The tolerance value used to determine the steady state based on the rolling variance.
#'
#' @return \code{Stable_gen}: A numeric value representing the generation (time step) where the system first reaches the steady state, or \code{NA} if no steady state is found.
#'
#' @importFrom zoo rollapply
#' 
#' @examples
#' # Example usage:
#'
#' wd <- "~/Documents/LAB_ECO"
#'
#' # Initial parameters
#' N_species = 2
#' C0 = 0.45
#' CN = 0.2
#' Diag_val = -0.5
#'
#' # Generate parameters
#' seeds_path <- file.path(wd, "Seeds.tsv" )
#' params <-  forge_data(N_species, seeds_path, C0, CN, Diag_val)
#'
#' # Generate simulation
#' output <- run_simulation(N_species, params = params, times = 20)
#'
#' tolerance = 0.0
#' all_rolling_var_SS(output, tolerance)
#' @export

all_rolling_var_SS <- function(output, tolerance) {

  # Ensure the zoo package is available
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("The 'zoo' package is required but not installed.")
  }

  #----------------------------Calculate Column Means-----------------#
  times <- ncol(output)  # Number of generations
  col_means <- colMeans(output, na.rm = TRUE)  # Compute column means

  #----------------------------Apply Moving Average---------------------#
  window_size <- round(times * 0.1)  # Calculate window size

  # Function to calculate moving variance
  # "left" covers following rows
  moving_var <- zoo::rollapply(col_means, width = window_size, FUN = var, fill = NA, align = "left")  # Calculate moving variance

  # Find the first generation where moving variance < tolerance
  Stable_gen <- which(moving_var < tolerance)[1] 

  return(Stable_gen)
}
