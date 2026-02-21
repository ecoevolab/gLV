#' Solve Model Using miaSim with Timeout
#'
#' The function leverages the miaSim package to simulate the model using the ode45 solver. It enforces a timeout of 600 seconds to prevent long-running simulations. Both relative and absolute tolerances for ode45 are set to 1e-06 to ensure precision. 
#'
#' @param params A list of parameters required by the gLV.
#' @param t_end End time of the simulation (aka. generations)
#' @param t_step Numeric: Interval between simulation
#' 
#'
#' @return A TreeSummarizedExperiment object containing:
#' 	- assays: A list with the transposed output matrix (species counts).
#' 	- colData: A DataFrame with the simulation time.
#' 	- metadata: A list containing migration_p and error_variance
#'
#' @seealso
#' \emph{miaSim} for model simulation and \emph{ode45} for solving ordinary differential equations.
#'
#' @examples
#' params <- list(x0 = vector of 3 species,
#'               M = M by M matrix,
#'               mu = vector of length 3 of unifor 0.001 to 1)
#'          
#' result <- solve_model_timeout(model, initial_state, params)
#' head(result)
#'
#' @export
mia.simulate <- function(params = params, n_t = n_t, t_step = t_step){
  
  # Check that matrrix is square
  if (nrow(params$M) != ncol(params$M)) {
    stop("Matrix is not n*n", call. = TRUE)
  }
  
 # Run simulation
  mia.res <- tryCatch(
    R.utils::withTimeout(miaSim::simulateGLV(n_species = nrow(params$M), 
                                             names_species = names(params$x0),
                                             A = params$M,
                                             x0 = params$x0,
                                             growth_rates = params$mu,
                                             sigma_migration = 0,
                                             epoch_p = 0,
                                             t_external_events = NULL,
                                             t_external_durations = NULL,
                                             stochastic = FALSE,
                                             migration_p = 0,
                                             error_variance = 0,
                                             norm = FALSE,
                                             t_end = n_t,
                                             t_step = 1,
                                             t_store = t_end),
                         timeout = 600), 
    error = function(e) {
      message(">> Simulation failed... skipping")
      return(NULL) # Return NULL instead of an NA matrix
    })
  
  # Ensure msim is valid before extracting counts
  if (inherits(mia.res, "TreeSummarizedExperiment")) {
    counts <- SummarizedExperiment::assay(mia.res, "counts")
    colnames(counts) <- round(mia.res$time, 1)
    return(counts)
  } else {
    return(matrix(NA, nrow = nrow(params$M), ncol = n_t)) # Return properly shaped NA matrix
  }
}
