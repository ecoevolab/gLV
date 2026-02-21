#' Simulate Extinctions of Survivor Species
#'
#' Simulates the extinction of each surviving species one by one, 
#' provided that at least 25% of the species have survived. 
#' The function runs a new simulation for each extinction event 
#' and calculates key ecological metrics such as Bray-Curtis 
#' dissimilarity and keystoneness.
#' 
#' @param output A dataframe containing the original simulation results using `ode45`. 
#' @param params A list containing the parameters required for the generalized Lotka-Volterra (gLV) model. 
#'
#' @return A tibble containing the following columns:
#' \itemize{
#'   \item \code{ID} - The simulation ID.
#'   \item \code{Specie} - The index of the species being extinct.
#'   \item \code{new_ext} - The number of additional extinctions caused by removing the species.
#'   \item \code{BC_diss} - The Bray-Curtis dissimilarity between the original and post-extinction communities.
#'   \item \code{K_s} - The keystoneness metric, measuring the impact of the extinction.
#'   \item \code{data} - A list column containing the resulting data frame of the extinction simulation.
#' }
#'
#' If fewer than 25% of species survive, the function returns `NULL` with a message.
#' 
#' @details
#' to access the df run the next dplyr function: ext_tibble[3,] %>% pull(data)
#' @examples
#' # Example usage of the function
#' result <- sim_all_ext(output, params)
#'

generate_extinctions <- function(params, path_core) {
  
  n_species <- params$n
  x_before <- params$x0                              # Non perturbated populations  
  keystoness_props <- x_before / sum(x_before)       # Proportions for keystoneness

  # Pre-allocate data frame
  df <- data.frame(
    spec = integer(n_species),
    new_ext = integer(n_species),
    BC_diss = numeric(n_species),
    K_s = numeric(n_species),
    ext_ts = numeric(n_species)
  )
  for(i in seq_len(n_species)) {
    # Remove species i 
    tmp_params <- list(
      x0 = params$x0[-i],
      mu = params$mu[-i],
      M = params$M[-i, -i, drop = FALSE]
    )
    # Run simulation
    new_out <- solve_gLV(times = 1000, tmp_params)
    x_after <- new_out[, ncol(new_out)]                 # Last column 
    # Section: Extinctions
    extinct_after <- x_after <= 1e-6                    # extinct after perturbation DIED
    extinct_before <- tmp_params$x0 > 1e-6              # extinct before perturbation WERE ALIVE
    new_ext <- sum(extinct_after & extinct_before)      # new extinctions
    # Section: Bray-Curtis dissimilarity 
    sum_x_after <- sum(x_after)
    sum_x_before <- sum(tmp_params$x0)
    bray_curtis <- 1 - (2 * sum(pmin(tmp_params$x0, x_after))) / (sum_x_before + sum_x_after)
    # Section: Keystoneness 
    K_s <- bray_curtis * (1 - keystoness_props[i])
    # Time to stability
    ext_ts <- find_stability(new_out)
    # Fill pre-allocated data frame
    df[i, ] <- list(i, new_ext, bray_curtis, K_s, ext_ts)
    # Save output
    ext_path <- file.path(path_core, sprintf("E_%s-S%d.feather", params$id, i))
    arrow::write_feather(as.data.frame(new_out), ext_path)
  }
  cat(">> Extinctions completed for", params$id, ".\n")
  return(df)
}


