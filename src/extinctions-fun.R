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

sim_all_ext <- function(output, params) {
  
  old_x0 <- output[[1000]] # Get population near steady state
  survivors <- which(old_x0 > 0)  # Indices of surviving species
  
  # Function to simulate extinction for a single species
  results <- lapply(survivors, function(spec) {
    
    params$x0[spec] <- 0 # Remove species from initial population

    # Rerun simulation
    ext_out <- solve_gLV(times = 1000, params)
    final_pop <- ext_out[[1000]]  # Extract final population

    # Compute Bray-Curtis dissimilarity
    bray_curtis <- 1 - (2 * sum(pmin(old_x0, final_pop))) / (sum(old_x0) + sum(final_pop))

    # Count secondary extinctions
    new_ext <- sum(final_pop == 0 & old_x0 != 0)

    # Compute keystoneness
    props <- final_pop / sum(final_pop) # final population proportions
    K_s <- bray_curtis * (1 - props[spec])

    # Time to Stability
    ext_ts <- find_ts(ext_out)
    list(new_ext = new_ext, BC_diss = bray_curtis, K_s = K_s, data = ext_out, spec = spec, ext_ts = ext_ts)
  })
  
  # Convert results into a tibble
  tib <- tibble::tibble(
    ID = params$id, # ID from params
    specie = purrr::map_dbl(results, "spec"), # Extract species info as numeric
    sub_ext = purrr::map_dbl(results, "new_ext"), # Subsequent extinctions
    BC_diss = purrr::map_dbl(results, "BC_diss"), # Bray-Curtis dissimilarity
    K_s = purrr::map_dbl(results, "K_s"), # Keystoneness
    ts_ext = purrr::map_dbl(results, "ext_ts"), # Time to stability 
    data = purrr::map(results, "data") # Extract data frames from results
  )
  
  cat("Extinctions done for simulation... ", params$id, "\n")
  return(tib)
}


