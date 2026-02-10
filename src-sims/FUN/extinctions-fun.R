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

sim_all_ext <- function(params, path_core) {
  # Extract parameters
  n = params$n
  x_before = params$x0                              # Non perturbated populations  
  rel_pop_initial <- x_before / sum(x_before)       # Proportions
  # Pre-allocate data frame
  exts_df <- data.frame(
    specie = integer(n),
    n_extinctions = integer(n),
    prop_extinctions = numeric(n),
    dissimilarity_bc = numeric(n),
    keystoneness = numeric(n),
    time_stability = numeric(n)
  )
  for(i in 1:n){
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
    extinct_after = x_after <= 1e-6                   # NOW DIED
    extinct_before = tmp_params$x0 > 1e-6                  # WERE ALIVE
    n_extinctions <- sum(extinct_after & extinct_before)    # new extinctions
    props_extinctions = n_extinctions/n                     # proportion of extinctions
    #------------------------------------
    # Section: Bray-Curtis dissimilarity 
    # Remove species i from the original community.
    # Calculate
    x_removed = x_before[-i]
    bray_curtis <- 1 - (2 * sum(pmin(x_removed, x_after))) / (sum(x_removed) + sum(x_after))
    #------------------------------------
    # Section: Keystoneness 
    props <- x_after / sum(x_after) # relative abundance
    keystoneness <- bray_curtis * (1 - props[i])
    #--------------------Time to stability-------------------------#
    time_stability <- find_ts(new_out)
    # 
    # Generate data frame
    exts_df[i, "specie"] <- i                                     # specie-extinct
    exts_df[i, "n_extinctions"] <- n_extinctions                  # new-extinctions
    exts_df[i, "prop_extinctions"] <- round(props_extinctions, 2) # proportion-extinctions
    exts_df[i, "dissimilarity_bc"] <- bray_curtis                 # Bray-Curtis
    exts_df[i, "keystoneness"] <- keystoneness                    # Keystoneness             
    exts_df[i, "time_stability"] <- time_stability                # Time-to-stability          
    # Lines to save the output of each extinction 
    # ext_path <- paste0(path_core, "/E_", params$id, "-S", i, ".feather")         
    # arrow::write_feather(x = new_out, sink = ext_path)                                       
  }
  # Add relative abundance of the extinct species before extinction
  exts_df$rel_pop_initial <- rel_pop_initial   
  cat(">> Extinctions completed for", params$id, ".\n")
  return(exts_df)
}


