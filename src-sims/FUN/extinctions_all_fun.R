#' Simulate Extinctions of all Species
#'
#' Simulates the extinction of each species one by one.
#' The function runs a new simulation for each extinction event 
#' and calculates key ecological metrics.
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

sim_all_ext <- function(params) {
  # Unperturbed community populations
  x0      <- params$x0
  rel_x0  <- x0 / sum(x0)
  # Perturbate community
  dplyr::bind_rows(lapply(1:params$n, function(i) {
    # Remove species i from the parameters
    tmp_params <- list(
      x0 = x0[-i],
      mu = params$mu[-i],
      M  = params$M[-i, -i, drop = FALSE]
    )

    new_out     <- solve_gLV(times = 1000, tmp_params)  # run simulation
    # population after extinction of species i
    x_after     <- new_out[, ncol(new_out)]            
    # population before extinction without species i
    x_before    <- tmp_params$x0                          
    # New extinctions after removing species i
    n_ext       <- sum(x_after < 1e-6 & x_before > 1e-6)
    bray_curtis <- 1 - (2 * sum(pmin(x_before, x_after))) / (sum(x_before) + sum(x_after))

    data.frame(
      specie           = i,
      pop_initial      = x0[i],
      rel_pop_initial  = rel_x0[i],
      n_extinctions    = n_ext,
      prop_extinctions = n_ext / length(x_before),
      dissimilarity_bc = bray_curtis,
      keystoneness     = bray_curtis * (1 - rel_x0[i]),
      time_stability   = find_ts(new_out)
    )
  }))
}


