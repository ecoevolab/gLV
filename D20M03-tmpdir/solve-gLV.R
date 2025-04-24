#' Solves the Generalized Lotka-Volterra (gLV) equation numerically with a fixed tolerance of `1e-06`.
#' The function utilizes the Dormand-Prince method as the solver.
#' 
#' @param times A numeric vector of times at which to simulate the ODE.
#' @param params A list containing the parameters required for the gLV model.
#'
#' @return A matrix or data frame with the solution of the gLV equations over the specified time points.
#' 
#' @examples
#' # Example usage of the function
#' result <- solve_gLV(times = 10, params = list(M = matrix(rnorm(3^2), nrow=3, ncol=3), mu = runif(3, 0, 1), x0 = runif(3, 0.001, 1)))

solve_gLV <- function(times, params) {
  
  
  # Define the equation
  glv_model <- function(t, x0, params) {
    
    r <- params$mu         # Growth rate vector
    A <- params$M          # Interaction matrix
    
    # Compute dx/dt for each species
    dx <- x0 * (r + A %*% x0)
    list(dx)
  }
  
  # Define extinction thershold function
  ext_thr <- function(t, y, params) {
    y[y < 10^-8 ] <- 0
    return(y)
  }
  
  time_seq <- seq(1, times, by = 1)  # Define the time sequence
  
  # Get solution
  results <- tryCatch(
    R.utils::withTimeout(
      deSolve::ode(y = params$x0, # Variable to evaluate
                   times = time_seq, # Sequence of times
                   func = glv_model, # gLV equation
                   parms = params, # parameters for gLV
                   method = "ode45", # Dormand-Prince method
                   rtol = 1e-06, # relative tolerance
                   atol = 1e-06,
                   events = list(func = ext_thr, time = time_seq)),
    timeout = 600), # 600 seconds cap (10 minutes)
    error = function(e) {
      message(">> Simulation failed... skipping")
      return(NULL) # Return NULL instead of an NA matrix
    })
  
  
  # Remove the first column (`time`) and transpose it so columns represent generations
  # Check for valid output and return transposed results
  if (!is.null(results) && ncol(results) > 1) {
    return(t(results[, -1]))
  } else {
    return(matrix(NA, nrow = nrow(params$M), ncol = times)) # Return properly shaped NA matrix
  }
}
