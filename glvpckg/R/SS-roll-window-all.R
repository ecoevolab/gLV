#' SS_roll_window_all Function
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
#' @importFrom zoo rollmean
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
#' params <-  init_data(N_species, seeds_path, C0, CN, Diag_val)
#'
#' # Generate simulation
#' times <- 20 # Define the number of generations
#' output <- run_simulation(N_species, params = params, times = times)
#'
#' tolerance = 0.05
#' result <- SS_roll_window_all(output, tolerance)
#' @export

SS_roll_window_all <- function(output, tolerance) {

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
  moving_avg <- zoo::rollmean(col_means, window_size, fill = NA, align = "left")

  # Find stability based on moving average
  Stable_gen <- which(moving_avg < tolerance)[1]  # Find the first generation where moving average < tolerance

  return(list("Tolerance" = tolerance,
              "Stable_gen" = Stable_gen,
              "ColMeans" = col_means)
  )
}
