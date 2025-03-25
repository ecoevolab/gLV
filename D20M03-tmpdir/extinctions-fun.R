#' Simulate Extinctions of Survivor Species
#'
#' Simulates the extinction of each surviving species one by one, 
#' provided that at least 25% of the species have survived. 
#' The function runs a new simulation for each extinction event 
#' and calculates key ecological metrics such as Bray-Curtis 
#' dissimilarity and keystoneness.
#' 
#' @param output_ode A dataframe containing the original simulation results using `ode45`. 
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
#' result <- simulate_all_extinctions(output_ode, params)
#'

simulate_all_extinctions <- function(output_ode, params) {
  
  # Get population near steady state
  survivors <- output_ode[, 700]
  min_survivors <- params$n * 0.25  # 25% of species must survive
  
  # Check if at least 25% of species survived
  if (sum(survivors > 0) < min_survivors) {
    message("Simulation ", params$id, " failed due to low survivors...\n")
    return(NULL)
  }
  
  
  # Set initial conditions
  params$x0 <- survivors
  surv_specs <- which(params$x0 > 0)  # Indices of surviving species
  
  # Function to simulate extinction for a single species
  simulate_extinction <- function(spec) {
    
    # Remove species from initial population
    new_x0 <- survivors
    new_x0[spec] <- 0
    params$x0 <- new_x0

    # Rerun simulation
    new_out <- solve_gLV(times = 500, params)
    vec2 <- new_out[, ncol(new_out)]  # Extract final population

    # Compute Bray-Curtis dissimilarity
    bray_curtis <- 1 - (2 * sum(pmin(survivors, vec2))) / (sum(survivors) + sum(vec2))

    # Count secondary extinctions
    new_ext <- sum(vec2 == 0) - sum(survivors == 0)

    # Compute keystoneness
    vec2_props <- vec2 / sum(vec2)
    K_s <- bray_curtis * (1 - vec2_props[spec])

    list(new_ext = new_ext, BC_diss = bray_curtis, K_s = K_s, data = new_out)
  }
  
  # Use `map()` to apply `simulate_extinction` to all surviving species
  results <- map(surv_specs, simulate_extinction)
  
  
  # Convert results into a tibble
  ext_tibble <- tibble(
    ID = params$id,
    Specie = surv_specs,
    new_ext = map_dbl(results, "new_ext"),
    BC_diss = map_dbl(results, "BC_diss"),
    K_s = map_dbl(results, "K_s"),
    data = map(results, "data")
  )
  
  cat("Extinctions done for simulation... ", params$id, "\n")
  return(ext_tibble)
}


