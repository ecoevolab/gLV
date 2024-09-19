#' Generate Parameters for Simulation
#'
#' This function generates the necessary parameters for a simulation, which can then be used in the \code{Simulate_output} function.
#'
#' @param N_species Numeric. The number of species involved in the simulation.
#' @param seeds_path Character. The file path to a TSV file containing all possible seeds to sample. The seeds should be prepared
#'   using the \code{Seed_generator} function.
#' @param C0 Numeric. The probability of no interaction between species (i.e., 0 interaction).
#' @param CN Numeric. The probability of negative interaction (<0) between species.
#' @param Diag_val Numeric. The diagonal values used in the interaction matrix.
#'
#' @return A list containing the following elements:
#'   \itemize{
#'     \item \code{Interactions}: A matrix representing the interaction values between species.
#'     \item \code{Growths}: A numeric vector representing the growth rates of the simulated species.
#'     \item \code{Population}: A numeric vector representing the initial abundances of the simulated species.
#'     \item \code{Seeds}: A numeric vector representing the seeds used to generate the populations, interaction matrix, and growth rates.
#'   }
#'
#' @importFrom data.table fread
#' @export
#'
#' @examples
#' # Example usage:
#' wd <- "~/Documents/LAB_ECO"
#' seeds_path <- file.path(wd, "Seeds.tsv")
#' params <- init_data(N_species = 2, seeds_path = seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#' print(params)


init_data <- function(N_species, seeds_path, C0, CN, Diag_val) {

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
  params <- list(Population = Pobl,
                 Interactions = inter,
                 Growths = Grow,
                 Seeds = seed)

  return(params)
}



