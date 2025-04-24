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
  
  # Define the gLV model
  glv_model <- function(t, x0, params) {
    x0[x0 < 1e-8] <- 0 # Ignore the effect of species with population below a threshold
    
    # dx/dt = X *(r + A * X)
    dx <- x0 * (params$mu + params$M %*% x0)
    list(dx)
  }
  
  time_seq <- seq_len(times)  # Times to simulate
  
  # Try solving the system with a timeout
  results <- tryCatch(
    R.utils::withTimeout(
      deSolve::ode(
        y = params$x0,
        times = time_seq,
        func = glv_model,
        parms = params,
        method = "ode45",
        rtol = 1e-6,
        atol = 1e-6
      ),
      timeout = 600
    ),
    error = function(e) {
      message(">> Simulation failed... skipping")
      return(NULL)
    }
  )
  
  # Process results if valid
  if (!is.null(results) && ncol(results) > 1) {
    tmp <- results[, -1] |>  # remove time column
      t() |>                 # transpose
      as.data.frame() |>     # convert to df
      # Populations that went extinct (no effect)
      dplyr::mutate(across(everything(), ~ replace(., . < 1e-8, 0))) 
    
    return(tmp)
  }
  
  # If results not valid: return NA matrix 
  matrix(NA, nrow = nrow(params$M), ncol = times)
}
