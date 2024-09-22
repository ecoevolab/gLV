#' Save Simulation Output to a File
#'
#' This function saves the provided simulation output to a TSV file with a specified unique ID in the given working directory.
#'
#' @param output A data frame or matrix containing the simulation output to be saved. The columns represent generations, and rows represent species.
#' @param uniqueID Character. A unique identifier for the file, generated using the \code{generate_uniqueID} function.
#' @param wd Character. The path to the working directory where the file will be saved. The output file will be saved in the path: \code{wd/Outputs/O_uniqueID.tsv}.
#'
#' @details
#' The function constructs the file path using the working directory (`wd`), the "Outputs" subdirectory, and the unique identifier (`uniqueID`). The output is saved in TSV format, where columns represent generations and rows represent species.
#'
#' @return A message is printed to the console indicating the file's ID and save location.
#' 
#' @import utils
#' 
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
#' output <- run_simulation(N_species = 2, params = params, times = times, norm = FALSE)
#'
#' # Generate unique ID
#' uniqueID <- forge_id(wd)
#'
#' # Save output to file
#' output_saver(output, uniqueID, wd)

output_saver <- function(output, uniqueID, wd) {
  
  # attach package
  requireNamespace("utils")

  # Output path
  out_path <- file.path(wd, "Outputs", paste0("O_", uniqueID, ".tsv") )

  # Save output data
  utils::write.table(output, file = out_path, sep = "\t", row.names = FALSE, col.names = TRUE)

  cat("Output saved with ID: ", uniqueID, "\n Path:", out_path)
}
