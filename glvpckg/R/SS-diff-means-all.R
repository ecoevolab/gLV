#' SS_diff_means_all Function
#'
#' This function searches for steady states in a simulation by analyzing the squared differences in the simulation results over time.
#'
#' @param output Matrix: Simulation results where rows represent species and columns represent time steps.
#' @param tolerance Numeric: Tolerance value for determining the steady state. It is transformed using \code{log(tolerance^2)}.
#'
#' @return A list containing:
#'   \itemize{
#'      \item \code{Diff_means}: Numeric vector. Column means of the squared differences between successive time steps, transformed using \code{log}.
#'      \item \code{Tolerance}: Numeric. Original tolerance value used.
#'      \item \code{Transformed_tolerance}: Numeric. Tolerance value transformed using \code{log(tolerance^2)}, displayed with 5 decimal places.
#'      \item \code{Stable_gen}: Numeric or NA. Generation where the mean of the differences < tolerance transformed.
#'  }
#'
#' @examples
#' # Example usage:
#'
#' wd <- "~/Documents/LAB_ECO"
#'
#' # Initial parameters
#' N_species = 2
#' C0 = 0.45
#' CN = 0.2
#' Diag_val = -0.5
#'
#' # Generate parameters
#' seeds_path <- file.path(wd, "Seeds.tsv" )
#' params <-  init_data(N_species, seeds_path, C0, CN, Diag_val)
#'
#' # Generate simulation
#' times <- 20 # Define the number of generations
#' output <- run_simulation(N_species, params = params, times = times, norm = FALSE)
#'
#' tolerance = 0.05
#' result <- SS_diff_means_all(output, tolerance)
#'
#' @export

SS_diff_means_all <- function(output, tolerance) {

  #------------------------Get differences------------------------#
  gens <- ncol(output) # Number of generations
  specs <- nrow(output) # Number of species
  og_tol <- tolerance
  tolerance <- log(tolerance^2) # Apply transformation to tolerance

  #--------------- Compute differences and square them-----------#
  Diff_output <- matrix(nrow = specs, ncol = gens - 1) # Preallocate matrix
  for (r in 1:specs) {
    v <- output[r, ]
    Diff_output[r, ] <- diff(v)^2
  }

  #---------------------------Create data frame------------------------#
  # Apply log transformation, replace 0 with 10
  ln_mat <- log(ifelse(Diff_output == 0, NA, Diff_output))

  rownames(ln_mat) <- paste("specie", 1:specs, sep = "") # Assign row names
  colnames(ln_mat) <- seq(1, gens - 1) # Assign column names

  #-------------Calculate column means and stable generation-------------------#
  Diff_mean <- colMeans(Diff_output, na.rm = TRUE) # Column means
  Stable_gen <- which(Diff_mean < tolerance)[1] # Generation where difference < tolerance

  if (is.na(Stable_gen)) {
    Stable_gen <- NA # Handle case where no generation meets the tolerance
  } else {
    Stable_gen <- Stable_gen + 1 # Adjust for 1-based index
  }

  return(list("Diff_means" = Diff_mean,
              "Tolerance" = og_tol,
              "Transformed_Tolerance" = format(tolerance, digits = 5, nsmall = 5),
              "Stable_gen" = Stable_gen)
         )
}


