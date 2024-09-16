#' Simulate_output Function
#'
#' This function simulates species interactions over a specified number of generations.
#'
#' @param n_species Numeric: Number of species involved in the simulation.
#' @param params List: A list including the following elements:
#'   \itemize{
#'     \item Interactions: If \code{NULL}, a random interaction matrix is generated using \code{runif(n = n_species, min = 0, max = 1)}.
#'     \item Growths: Numeric vector representing the growth rates of the simulated species.
#'     \item Population: Numeric vector representing the initial abundances of the simulated species.
#'     \item Seeds: Numeric vector representing the seeds used to generate the interaction matrix, growth rates, and populations.
#'   }
#' @param times Numeric: The number of generations to simulate.
#' @param norm Logical: Whether the output time series should be normalized to proportions (\code{norm = TRUE}) or returned as raw counts (\code{norm = FALSE}, default).
#'
#' @return The function returns the simulation results over time, typically in the form of a matrix or data frame.
#' @export
#'
#' @examples
#' # Example usage:
#' wd <- "~/Documents/LAB_ECO" # Working directory
#'
#' # Generate parameters
#' seeds_path <- file.path(wd, "Seeds.tsv" )
#' # It can also be: "~/Documents/LAB_ECO/Seeds.tsv"
#' params <- generate (N_species = 2, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#'
#' # Simulate species interactions
#' times <- 100 # Define the number of generations
#' Output <- Simulate_output(N_species = 2, params = params, times = times, norm = FALSE)

Simulate_output <- function(N_species, params, times, norm) {

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
                                norm = norm) # FALSE=conteo, TRUE=proporciones

  output <- glvmodel@assays@data@listData[["counts"]]

  return(output)
}



#' Output_save function
#'
#' This function saves the provided output data to a TSV file with a specified ID in a given working directory.
#'
#' @param output A data frame or matrix containing the simulation output to be saved.
#' @param uniqueID A character string containing a unique ID to be used for saving files. Generated previously with \code{generate_uniqueID} function.
#' @param wd Character. Path to the working directory where files will be saved. The file will be saved in the path: \code{wd/Outputs/O_uniqueID.tsv}.
#'
#' @details
#' The function constructs the file path using the working directory (`wd`), the "Outputs" subdirectory,
#' and the unique identifier (`ID`). The file is saved with a ".tsv" extension where the columns represent generations
#' and rows species.
#'
#' @return A message indicating the ID and the file path is printed to the console.
#' @export
#'
#' @examples
#'
#' wd <- "~/Documents/LAB_ECO"
#'
#' # Generate parameters
#' seeds_path <- file.path(wd, "Seeds.tsv" )
#' # It can also be: "~/Documents/LAB_ECO/Seeds.tsv"
#' params <- generate (N_species = 2, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#'
#' # Simulate species interactions
#' times <- 100 # Define the number of generations
#' Output <- Simulate_output(N_species = 2, params = params, times = times, norm = FALSE)
#'
#' # Generate unique ID
#' uniqueID <- generate_uniqueID(wd)
#'
#' # Save output
#' output_save(Output, uniqueID, wd)

output_save <- function(Output, uniqueID, wd) {

  # Output path
  out_path <- file.path(wd, "Outputs", paste0("O_", uniqueID, ".tsv") )

  # Save output data
  write.table(Output, file = out_path, sep = "\t", row.names = FALSE, col.names = TRUE)

  cat("Output saved with ID: ", uniqueID, "\n Path:", out_path)
}

