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
  
  # Function to simulate extinction for a single species
  x = params$x0
  df = data.frame()
  for(i in seq_along(x)){
    x_new = x
    x_new[i] = 0                # extinct-specie i
    params$x0 = x_new           # update params
    new_out = solve_gLV(times = 1000, params)
    x_end = new_out[,1000]

    # Compute Bray-Curtis dissimilarity
    bray_curtis <- 1 - (2 * sum(pmin(x, x_end))) / (sum(x) + sum(x_end))

    # Count secondary extinctions
    # All that was live before (x > 1e-06)
    # But now is dead (x <= 1e-06)
    new_ext <- sum(x_end <= 1e-6 & x > 1e-6)

    # Compute keystoneness
    props <- x_end / sum(x_end) # final population proportions
    K_s <- bray_curtis * (1 - props[i])

    # Time to Stability
    ext_ts <- find_ts(new_out)

    row <- data.frame(
      spec = i,                     # specie-extinct
      new_ext = new_ext,            # new-extinctions
      BC_diss = bray_curtis,        # Bray-Curtis
      K_s = K_s,                    # Keystoness
      ext_ts = ext_ts               # Time-to-stability
    )
    df = rbind(df, row)
    ext_path <- paste0(path_core, "/E_", params$id, "-S", i, ".feather")          # Extinctions-paths
    arrow::write_feather(new_out, ext_path)                                       # Save-extinctions
  }
  cat(">> Extinctions completed for", params$id, ".\n")
  return(df)
}


