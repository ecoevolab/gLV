#' Steady State Search by Proportional Change in Column Means
#'  
#' @details
#' This function searches for steady states in a simulation by analyzing the changes in the logarithm of the mean simulation results over time. 
#' It leverages the logarithmic property \eqn{\log\left(\frac{x}{y}\right) = \log(x) - \log(y)} to simplify the comparison of successive time steps.
#'
#' @param output Matrix: A matrix of simulation results, where rows represent species and columns represent time steps. Typically generated with \link{run_simulation}.
#' @param tolerance Numeric: The tolerance value used to determine steady states. It compares the differences in consecutive column means over time.
#'
#' @return \code{Stable_gen}: A numeric value representing the first generation where the difference between means is less than the tolerance, or \code{NA} if no steady state is found.
#' @examples
#' # Example usage:
#' wd <- "~/Documents/LAB_ECO/Simulations"
#' seeds_path <- file.path(wd, "Seeds.tsv")
#' params <- init_data(N_species = 2, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#' 
#' # Run simulation for 100 generations
#' times <- 100  
#' output <- run_simulation(N_species = 2, params = params, times = times)
#'
#' # Search for steady states with a tolerance of 0.05
#' tolerance <- 0.05
#' all_prop_SS(output, tolerance)
#'
#' @export

all_prop_SS <- function(output, tolerance) {
  
  # Compute column means directly
  log_output <- ifelse(output <= 0, NA, log(output) )   # Apply logarithm safely
  log_mean <- colMeans(log_output)
  
  # Compute differences between consecutive means
  differ_mat <- diff(log_mean)
  
  # Find the first generation where the difference is below tolerance
  stable_gen <- which(abs(differ_mat) < tolerance)[1]
  
  return(stable_gen)
}


