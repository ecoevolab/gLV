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
  # Extract parameters
  n = params$n
  x_before = params$x0                        # Non perturbated populations  
  rel_x_before <- x_before / sum(x_before)    # Do them proportions
  exts_df <- do.call(rbind, lapply(1:n, function(i) { 
    # Remove species i from the community
    tmp_params <- list(
      x0 = params$x0[-i],
      mu = params$mu[-i],
      M = params$M[-i, -i, drop = FALSE]
    )
    #------------------------------------
    # Run simulation
    new_out = solve_gLV(times = 1000, tmp_params)
    x_after = new_out[, ncol(new_out)]                 # Last column 
    #------------------------------------
    # Section: Extinctions
    # Species already extincted
    alive_before = tmp_params$x0 > 1e-6                # were alive
    extinct_after = x_after <= 1e-6                      # now died
    # NUmber of extinctions
    # testing lines
    # extinct_after = c(TRUE, TRUE, FALSE)
    # extinct_before = rep(TRUE, 3)
    n_extinctions <- sum(extinct_after & alive_before)    # new extinctions
    props_extinctions = n_extinctions/length(extinct_after) # proportion of extinctions
    #------------------------------------
    # Section: Bray-Curtis dissimilarity
    # Remove species i from the original community.
    bray_curtis <- 1 - (2 * sum(pmin(tmp_params$x0, x_after))) / (sum(tmp_params$x0) + sum(x_after))
    #------------------------------------
    # Section: Keystoneness 
    keystoneness <- bray_curtis * (1 - rel_x_before[i])
    #------------------------------------
    # Section: Time to stability 
    time_stability <- find_ts(new_out)
    #------------------------------------
    # Generate data frame
    data.frame(
        specie           = i,
        n_extinctions    = n_extinctions,
        prop_extinctions = round(props_extinctions, 2),
        dissimilarity_bc = bray_curtis,
        keystoneness     = keystoneness,
        time_stability   = time_stability
        )                                            
    }))
  # Add relative abundance of the extinct species before extinction
  exts_df$pop_initial = x_before
  exts_df$rel_pop_initial = rel_x_before   
  # cat(">> Extinctions completed for", params$id, ".\n")
  return(exts_df)
}


