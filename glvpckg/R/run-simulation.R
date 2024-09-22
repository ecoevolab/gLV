#' Simulate_output Function
#'
#' This function simulates species interactions over a specified number of generations using miaSim::simulateGLV package
#'
#' @param N_species Numeric: Number of species involved in the simulation.
#' @param params List: A list including the following elements:
#'   \itemize{
#'     \item Interactions: If \code{NULL}, a random interaction matrix is generated using \code{runif(n = n_species, min = 0, max = 1)}.
#'     \item Growths: Numeric vector representing the growth rates of the simulated species.
#'     \item Population: Numeric vector representing the initial abundances of the simulated species.
#'     \item Seeds: Numeric vector representing the seeds used to generate the interaction matrix, growth rates, and populations.
#'   }
#' @param times Numeric: The number of generations to simulate.
#'
#' @return The function returns the simulation results over time, typically in the form of a matrix or data frame.
#' 
#' @importFrom miaSim simulateGLV
#' @export
#'
#' @examples
#' wd <- "~/Documents/LAB_ECO"
#'
#' # Generate parameters
#' seeds_path <- file.path(wd, "Seeds.tsv")
#' params <- init_data(N_species = 2, seeds_path = seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#'
#' # Simulate species interactions
#' times <- 100  # Define the number of generations
#' output <- run_simulation(N_species = 2, params = params, times = times)

run_simulation <- function(N_species, params, times) {

  # Ensure the ids package is available
  if (!requireNamespace("miaSim", quietly = TRUE)) {
    stop("The 'miaSim' package is required but not installed.")
  }

  glvmodel <- miaSim::simulateGLV(n_species = N_species,
                                A = params$alpha, # interaction matrix
                                x0 = params$Pobl, # Initial abundances
                                growth_rates = params$r, # Growth rates
                                t_start = 0,
                                t_store = times,
                                t_end= times,
                                migration_p = 0,
                                stochastic = FALSE, # Ignorar ruido
                                norm = TRUE) # FALSE=conteo, TRUE=proporciones

  output <- glvmodel@assays@data@listData[["counts"]]

  return(output)
}




