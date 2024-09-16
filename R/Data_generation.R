#' Generate a UniqueID function
#'
#' This function generates a unique ID by checking the output directory to ensure no duplicate IDs are created.
#'
#' @param wd Character: The working directory path where output files are stored.
#'   The function will generate a unique ID based on the existing files in the "Outputs" subdirectory.
#'
#' @return A character string containing a unique ID to be used for saving files.
#' @export
#'
#' @examples
#' # Example usage:
#' wd <- "~/Documents/LAB_ECO"
#' uniqueID <- generate_uniqueID(wd)

generate_uniqueID <- function(wd) {

  # Ensure the ids package is available
  if (!requireNamespace("ids", quietly = TRUE)) {
    stop("The 'ids' package is required but not installed.")
  }

  repeat {
    directory = file.path(wd,"Outputs")
    ID <- ids::random_id(1, 3)
    all_files <- list.files(directory, full.names = TRUE) # List all files in the directory
    matching_files <- grep(ID, basename(all_files), value = TRUE) # Check if any file name contains the ID as a substring

    # Check if any file matches the pattern
    if (length(matching_files) == 0) {
      return(ID = ID) # Return the unique ID
    }
  }
}

#' Generate Function
#'
#' This function generates the parameters required for simulation, which can then be used in the \code{Simulate_output} function.
#'
#' @param N_species Numeric. Number of species involved in the simulation.
#' @param seeds_path Character. Path to a TSV file containing all possible seeds to sample. The seeds should be prepared
#'   using the \code{Seed_generator} function.
#' @param C0 Numeric. Probability of no interaction between species (i.e., 0 interaction).
#' @param CN Numeric. Probability of negative interaction (<0) between species.
#' @param Diag_val Numeric. Diagonal values used in the interaction matrix.
#'
#' @return List. A list containing the following elements:
#'   \itemize{
#'     \item \code{Interactions}: A matrix representing the interaction values between species.
#'     \item \code{Growths}: A numeric vector representing the growth rates of the simulated species.
#'     \item \code{Population}: A numeric vector representing the initial abundances of the simulated species.
#'     \item \code{Seeds}: A numeric vector representing the seeds used to generate the interaction matrix, growth rates, and populations.
#'   }
#'
#' @examples
#' # Example usage:
#' wd <- "~/Documents/LAB_ECO"
#' seeds_path <- file.path(wd, "Seeds.tsv")
#' params <- generate(N_species = 2, seeds_path = seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#'
#' @export

generate <- function(N_species, seeds_path, C0, CN, Diag_val) {

  #------------------Read Seeds------------------------------#
  library(data.table)
  seeds <- as.matrix(fread(seeds_path, sep = "\t"))

  #------------------Populations-----------------------------#
  S_p <- sample(seeds, 1)
  set.seed(S_p)
  Pobl <- runif(N_species, min = 0.1, max = 1)

  #--------------------Interactions-------------------------#
  S_i <- sample(seeds, 1)
  set.seed(S_i)

  # Generate random interaction values
  P_neg <- rbinom(N_species * N_species, 1, CN)
  tmp <- rbinom(N_species * N_species, 1, 1 - C0) * ifelse(P_neg != 0, runif(N_species * N_species, min=0, max=1), -runif(N_species * N_species, min=0, max=1))
  inter <- matrix(tmp, nrow = N_species, ncol = N_species)

  # Set diagonal values
  diag(inter) <- Diag_val

  #------------------------Growth Rates---------------------#
  S_g <- sample(seeds, 1)
  set.seed(S_g)
  Grow <- runif(N_species, min = 0.001, max = 1)

  # Collect seeds used
  seed <- c(S_p, S_i, S_g)

  # Return parameters as a list
  params <- list(Interactions = inter,
                 Growths = Grow,
                 Population = Pobl,
                 Seeds = seed)

  return(params)
}



