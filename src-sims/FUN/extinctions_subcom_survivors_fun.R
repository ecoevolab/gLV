#' Simulate Extinctions of only survival species
#'
#' Simulates the extinction of each species one by one.
#' The function runs a new simulation for each extinction event 
#' and calculates key ecological metrics.
#' 
#' The impact is calculated with the subcommunity composed of survival species
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
#' If fewer species than threshold survive, the function returns `NULL` with a message.
#' 
#' @examples
#' # Example usage of the function
#' result <- sim_all_ext(output, params)
#'

sim_ext_surv_sub <- function(params, ext_threshold) {
    #------------------------------------
    # Section: Extract parameters
    sim_id = params$id
    x_before_full = params$x0                        # Non perturbated populations  
    rel_x_before_full <- x_before_full / sum(x_before_full)    # Do them proportions
    to_filter <- which(rel_x_before_full > 1e-06)    # Which species to filter
    if (!(length(to_filter) >= ext_threshold)) {
        cat(paste0('>> Skipping sub-community impact extinctions for id ', sim_id, '. Not enough species (', length(to_filter),') passed the filter ', 
        ext_threshold ,'\n'))  
        return(NULL)
    }
    #------------------------------------
    # Section: Generate sub-community composed of survival nodes
    subcom_params = list(x0 = params$x0[to_filter], # for extinctions
        M = params$M[to_filter,to_filter],
        mu = params$mu[to_filter], 
        id = sim_id, n = length(to_filter)
    )
    # Recalculate relative abundance
    rel_subcom = subcom_params$x0/sum(subcom_params$x0)
    #------------------------------------
    # Section: Extract parameters
    exts_df <- do.call(rbind, lapply(seq_along(to_filter), function(i) { 
        # Remove species i from the community
        to_extinct = to_filter[i]
        tmp_params <- list(
            x0 = subcom_params$x0[-i],
            mu = subcom_params$mu[-i],
            M = subcom_params$M[-i, -i, drop = FALSE]
        )
        #------------------------------------
        # Run simulation
        new_out = solve_gLV(times = 1000, tmp_params)
        # Population after perturbation
        x_after = new_out[, ncol(new_out)]                 
        #------------------------------------
        # Section: Extinctions
        # Species already extincted
        alive_before = tmp_params$x0 > 1e-6  # were alive
        extinct_after = x_after < 1e-6      # now died
        # Number of new extinctions
        n_extinctions <- sum(extinct_after & alive_before)    
        props_extinctions = n_extinctions/length(tmp_params$x0) # proportions
        #------------------------------------
        # Section: Bray-Curtis dissimilarity
        x_before = tmp_params$x0    # Population before perturbation
        bray_curtis <- 1 - (2 * sum(pmin(x_before, x_after))) / (sum(x_before) + sum(x_after))
        #------------------------------------
        # Section: Keystoneness 
        keystoneness <- bray_curtis * (1 - rel_subcom[i])
        #------------------------------------
        # Section: Time to stability 
        time_stability <- find_ts(new_out)
        #------------------------------------
        # Generate data frame
        data.frame(
            specie           = to_extinct,
            n_extinctions    = n_extinctions,
            prop_extinctions = round(props_extinctions, 2),
            dissimilarity_bc = bray_curtis,
            keystoneness     = keystoneness,
            time_stability   = time_stability
        )                                            
    }))
    # Add relative abundance of the extinct species before extinction
    exts_df$pop_initial = subcom_params$x0
    exts_df$rel_pop_initial = rel_subcom
    # cat(">> Extinctions completed for", params$id, ".\n")
    return(exts_df)       
}