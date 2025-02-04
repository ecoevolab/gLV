#' Save Parameters by Lines
#'
#' This function saves simulation parameters line by line to a specified file for future reference.
#'
#' @param params List. A list generated with the \link{forge_data} function, containing the following elements:
#' \itemize{
#'   \item \code{Interactions}: A matrix representing the interaction values between species.
#'   \item \code{Growths}: A numeric vector representing the growth rates of the simulated species.
#'   \item \code{Population}: A numeric vector representing the initial abundances of the simulated species.
#'   \item \code{Seeds}: A numeric vector representing the seeds used to generate the populations, interaction matrix, and growth rates.
#' }
#' @param uniqueID Character. A unique ID string used for saving files. This ID should be generated previously using the \code{forge_ID} function.
#' @param wd Character. The path to the working directory where the file will be saved. The file will be saved at: `wd/Parameters/P_uniqueID.tsv`.
#'  
#' @return A message indicating the path where the parameters file was stored.
#'
#' @import utils
#' @examples
#' # Example usage:
#' wd <- "~/Documents/LAB_ECO/Simulations"
#' 
#' # Generate parameters
#' seeds_path <- file.path(wd, "Seeds.tsv")
#' params <- init_data(N_species = 10, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#' 
#' # Generate unique ID
#' uniqueID <- forge_id(wd)
#' 
#' # Save parameters by line
#' params_line_saver(params, uniqueID, wd)
#' 
#' @export


params_line_saver <- function(params, uniqueID, wd) {
  
  # attach package 
  requireNamespace("utils")

  # Generate unique file paths
  pms_path <- file.path(wd, "Parameters", paste0("P_", uniqueID, ".tsv") )

  # Separate parameters
  Interactions <- params$Interactions
  Growths <- params$Growths
  Population <- params$Population
  Seeds <- params$Seeds

  # Create a function to save each section with headers
  save_section <- function(title, data, pms_path) {
    con <- file(pms_path, open = "at") # Open connection in append mode
    on.exit(close(con))
    cat("\n", title, "\n", file = pms_path, append = TRUE)
    utils::write.table(data, file = pms_path, sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
    writeLines("", con) # Add a new line
  }

  # Save each section
  save_section("Interactions", Interactions, pms_path)
  save_section("Grow rates", Growths, pms_path)
  save_section("Initial populations", Population, pms_path)
  save_section("Seeds", Seeds, pms_path)

  cat("Lines saved with ID: ", uniqueID, "\n Path used:", pms_path)
}
