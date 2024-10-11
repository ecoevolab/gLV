#' Generate Parameters for Simulation
#'
#' This function generates the necessary parameters for a simulation, which can then be used in the \link{run_simulation} function.
#'
#' @param N_species Numeric. The number of species involved in the simulation.
#' @param seeds_path Character. The file path to a TSV file containing all possible seeds to sample. The seeds should be prepared
#'   using the \link{forge_seeds} function.
#' @param C0 Numeric. The probability of no interaction between species (i.e., 0 interaction).
#' @param CN Numeric. The probability of negative interaction (<0) between species.
#' @param Diag_val Numeric. The diagonal values used in the interaction matrix.
#'
#' @details 
#' \itemize{
#'   \item For \strong{generating the populations}, a uniform distribution between 0.1 and 1 is used. This represents the proportion of each species in the medium.
#'   \item For \strong{generating the interaction matrix}, a binomial distribution with a success probability of `CN` is used to determine whether the interaction will be negative.
#'   
#'   The resulting number is then multiplied by another binomial distribution with a success probability of `C0` to determine if the interaction will be null (0). Finally, this product is multiplied by a uniform distribution between 0 and 1 to obtain the final interaction value.
#'   
#'   The diagonal of the interaction matrix is filled with the value of `Diag_val` (user input).
#'   \item For \strong{generating the growth rates}, a uniform distribution between 0.001 and 1 is utilized. This represents the growth rate of each species.
#' }
#' 
#' @return A list containing the following elements:
#'   \itemize{
#'     \item \code{Interactions}: A matrix representing the interaction values between species.
#'     \item \code{Growths}: A numeric vector representing the growth rates of the simulated species.
#'     \item \code{Population}: A numeric vector representing the initial abundances of the simulated species.
#'     \item \code{Seeds}: A numeric vector representing the seeds used to generate the populations, interaction matrix, and growth rates.
#'   }
#'
#' @import stats
#' @importFrom data.table fread
#' @export
#'
#' @examples
#' # Example usage:
#' wd <- "~/Documents/LAB_ECO/Simulations"
#' 
#' # Generate parameters
#' seeds_path <- file.path(wd, "Seeds.tsv")
#' params <- forge_data(N_species = 10, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)


forge_data <- function(N_species, seeds_path, C0, CN, Diag_val) {
  
  # Attach package
  requireNamespace("stats")

  #------------------Read Seeds------------------------------#
  requireNamespace("data.table")
  seeds <- as.matrix(data.table::fread(seeds_path, sep = "\t"))

  #------------------Populations-----------------------------#
  S_p <- sample(seeds, 1)
  set.seed(S_p)
  Pobl <- stats::runif(N_species, min = 0.1, max = 1)

  #--------------------Interactions-------------------------#
  S_i <- sample(seeds, 1)
  set.seed(S_i)

  # Generate random interaction values
  P_neg <- stats::rbinom(N_species * N_species, 1, CN)
  tmp <- stats::rbinom(N_species * N_species, 1, 1 - C0) * ifelse(P_neg != 0, stats::runif(N_species * N_species, min = 0, max = 1), -stats::runif(N_species * N_species, min = 0, max = 1))
  inter <- matrix(tmp, nrow = N_species, ncol = N_species)

  # Set diagonal values
  diag(inter) <- Diag_val

  #------------------------Growth Rates---------------------#
  S_g <- sample(seeds, 1)
  set.seed(S_g)
  Grow <- stats::runif(N_species, min = 0.001, max = 1)

  # Collect seeds used
  seed <- c(S_p, S_i, S_g)

  # Return parameters as a list
  params <- list(Population = Pobl,
                 Interactions = inter,
                 Growths = Grow,
                 Seeds = seed)

  return(params)
}



