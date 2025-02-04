#' Save Parameters by Seed for Future Regeneration
#'
#' This function saves specified parameters in a tabular format. The parameters include the number of species, interaction probabilities, diagonal values of the interaction matrix, and seed values used in different aspects of the simulation.
#'
#' @param N_species Integer. The number of species involved in the simulation.
#' @param C0 Numeric. The probability of no interaction between species (i.e., zero interaction).
#' @param CN Numeric. The probability of negative interaction (<0) between species.
#' @param Diag_val Numeric vector. The diagonal values used in the interaction matrix.
#' @param params List. Parameters used in the simulation, specifically \code{params$Semilla} will be saved, containing:
#' \itemize{
#'   \item \code{Semilla[1]}: Seed for initial populations.
#'   \item \code{Semilla[2]}: Seed for the interaction matrix.
#'   \item \code{Semilla[3]}: Seed for species growth rates.
#' }
#' @param uniqueID Character. A unique ID string used for saving files. This ID should be generated previously using the \code{forge_ID}function.
#' @param wd Character. The path to the working directory where files will be saved. The file will be saved at: `wd/Parameters/Seeds_save.tsv`.
#'
#' @details
#' This function stores the provided parameters in a tabular format for future reference. The table includes the species count, interaction probabilities,
#' matrix diagonal values, and seed values.
#'  
#' @import utils
#' 
#' @examples
#' wd <- "~/Documents/LAB_ECO/Simulations"
#' 
#' # Generate parameters
#' seeds_path <- file.path(wd, "Seeds.tsv")
#' params <- init_data(N_species = 10, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#' 
#' # Generate unique ID
#' uniqueID <- forge_id(wd)
#'
#' # Save parameters by seeds
#' params_seed_saver(N_species = 10,  C0 = 0.45, CN = 0.2, Diag_val = -0.5, params, uniqueID, wd)
#'
#' @export

params_seed_saver <- function(N_species, C0, CN, Diag_val, params, uniqueID, wd) {
  
  # attach package 
  requireNamespace("utils")

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
    existing_data <- utils::read.delim(Ps_path, sep = "\t", header = TRUE)
    combined_data <- rbind(existing_data, seeds_df)
  } else {
    combined_data <- seeds_df
  }

  # Save table
  utils::write.table(combined_data, file = Ps_path, sep = "\t", row.names = FALSE, col.names = TRUE)

  cat("Parameters seeds saved \n\tID:", uniqueID, "\n\tPath:", Ps_path, "\n")
}


