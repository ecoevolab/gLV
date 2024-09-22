#' Seed_saver Function
#'
#' This function saves specified parameters to a table format. The parameters include the number of species, interaction probabilities, diagonal values of the interaction matrix, and seed values used in different aspects of the simulation.
#'
#' @param N_specs Integer. Number of species involved in the simulation.
#' @param C0 Numeric. Probability of no interaction between species (i.e., zero interaction).
#' @param CN Numeric. Probability of negative interaction (<0) between species.
#' @param Diag_val Numeric vector. Diagonal values used in the interaction matrix.
#' @param params List. Parameters used in the simulation, specifically \code{params$Semilla} containing:
#' \itemize{
#'   \item \code{Semilla[1]}: Seed for initial populations.
#'   \item \code{Semilla[2]}: Seed for the interaction matrix.
#'   \item \code{Semilla[3]}: Seed for species growth rates.
#' }
#' @param uniqueID Character. A unique ID string used for saving files. This ID should be generated previously using the \code{generate_uniqueID} function.
#' @param wd Character. Path to the working directory where files will be saved. The file will be saved in the path: \code{wd/Parameters/Seeds_save.tsv}.
#'
#' @details
#' The function stores the provided parameters in a tabular format for future reference. The table includes species count, interaction probabilities,
#' matrix diagonal values, and seed values.
#'
#' @examples
#'  wd <- "~/Documents/LAB_ECO"
#'
#'  # Generate parameters for simulation
#'  seeds_path <- file.path(wd, "Seeds.tsv" )
#'  params <- init_data(N_species = 2, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#'
#'  # Generate unique ID
#'  uniqueID <- forge_ID(wd)
#'
#'  # Save parameters using Seed_saver function
#'  seed_saver(N_specs = 2,  C0 = 0.45, CN = 0.2, Diag_val, params, uniqueID, wd)
#'
#' @export

seed_saver <- function(N_species, C0, CN, Diag_val, params, uniqueID, wd) {

  # Save seeds information
  seeds_df <- data.frame(
    'ID_simulation' = uniqueID,
    "N_specs" = N_species,
    "Prob_0" = C0,
    "Prob_neg" = CN,
    "Diagonal" = Diag_val,
    "Population_seed" = params$Seeds[1],
    "Interactions_seed" = params$Seeds[2],
    "Growth_seed" = params$Seeds[3]
  )

  Ps_path <- file.path(wd, "Parameters", "Seeds_save.tsv")

  if (file.exists(Ps_path)) {
    existing_data <- read.delim(Ps_path, sep = "\t", header = TRUE)
    combined_data <- rbind(existing_data, seeds_df)
  } else {
    combined_data <- seeds_df
  }

  # Save table
  write.table(combined_data, file = Ps_path, sep = "\t", row.names = FALSE, col.names = TRUE)

  cat("Seeds saved with ID: ", uniqueID, "\n Path:", Ps_path)
}


