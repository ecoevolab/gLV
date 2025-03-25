#' Generate gLV parameters with normal distribution
#'
#' This function generates simulation parameters from a TSV file containing the required 
#' parameters and associated IDs. The TSV file must be previously created using the `Grid method`.
#'
#' @param n_species Integer. Number of species to simulate.
#'
#' @param p_noint Numeric. Probability of no interaction, i.e., the probability that a 
#' non-diagonal element in the interaction matrix is zero.
#'
#' @details
#' Each parameter is generated as follows:
#'
#' - `x0`: Initial population of each species at time 0, drawn from a uniform distribution between 0 and 1.
#'
#' - `mu`: Growth rate of each species across all generations, drawn from a uniform distribution between 0.001 and 1.
#'
#' - `M`: Interaction matrix with diagonal elements set to `d` and non-diagonal elements drawn from a normal distribution with mean `Norm_mu` and standard deviation `sigma`. Non-diagonal elements may be zero with probability `p_noint`.
#'
#' @return A lsit containing the parameters required for gLV equation.

regenerate_Dnormal <- function(index) {
  
  n_species <- as.numeric(index[["n_species"]])
  
  #------------------Populations-----------------------------#
  set.seed(as.numeric(index[["x0_seed"]]))
  x0 <- stats::runif(n_species, min = 0.1, max = 1)
  
  #------------------------Growth Rates---------------------#
  set.seed(as.numeric(index[["mu_seed"]]))
  mu <- stats::runif(n_species, min = 0.001, max = 1)
  
  #--------------------Interactions-------------------------#
  
  # Create matrix with NA values and fill the diagonal with -0.5
  M <- matrix(NA, nrow = n_species, ncol = n_species)
  diag(M) <- as.numeric(index[["diagonal"]])
  
  #' Define proportions for null interactions
  p_noint <- as.numeric(index[["p_noint"]])
  
  # Define the number of interactions
  num_off_diag <- n_species * (n_species - 1)  # Total off-diagonal elements
  num_noint <- floor(p_noint * num_off_diag)   # Number of null interactions
  num_int <- num_off_diag - num_noint
  
  # Create the interaction vector
  set.seed(as.numeric(index[["A_seed"]]))
  mu <- as.numeric(index[["Norm_mu"]])
  sd <- as.numeric(index[["sigma"]])
  interaction_values <- c(rep(0, num_noint),
                          rnorm(n = num_int , mean = 0, sd = 0.5)
  )
  
  # Shuffle the interaction vector
  interaction_values <- sample(interaction_values)
  
  # Assign to off-diagonal elements
  M[upper.tri(M, diag = FALSE) | lower.tri(M, diag = FALSE)] <- interaction_values
  
  # Optional: Round if needed
  M <- round(M, digits = 4)
  
  # Extract ID
  id <- index[["id"]]
  
  # Return parameters as a list
  params <- list(x0 = x0,
                 M = M,
                 mu = mu,
                 id = id,
                 n = n_species)
  
  return(params)
}