#' Find Stability Time for a Community
#'
#' This function calculates the **stability time** for a community based on the frequencies of species over time.
#' Stability is defined as the point where the absolute difference between successive frequencies of species 
#' is consistently less than or equal to a specified threshold (0.05) for at least 10 consecutive time steps.
#'
#' The function checks the frequencies of each species over time and identifies the first time when the community 
#' reaches a stable state for the given condition. The output is the time at which stability is reached, defined 
#' as the maximum time where at least 10 successive time points have small changes in frequency.
#'
#' @param output A matrix or data frame where rows correspond to species' frequencies 
#' over time and columns correspond to time points.The values in the table represent the frequencies of the species at each time step.
#'
#' @return A numeric value representing the **stability time** of the community. This is the time at which the species' frequencies 
#'         stabilize, based on the condition of having 10 successive time points with absolute changes less than or equal to 0.05.
#'
find_ts <- function(output) {
  
  # Apply the logic across each row of output
  x <- vapply(1:nrow(output), function(i) {
    row <- as.numeric(output[i, ])
    
    # Identify differences less than or equal to 0.05
    test <- abs(diff(row)) <= 0.05
    rle_res <- rle(test)  # Run-length encoding
    
    # Find the positions where TRUE values have length >= 10
    valid_runs <- which(rle_res$values == TRUE & rle_res$lengths >= 10)
    
    # Identify the starting positions (adjusted for diff)
    if (length(valid_runs) > 0) {
      t <- sum(rle_res$length[-valid_runs]) + 1
    } else {
      t <- NA  # In case no valid runs are found
    }
    return(t)
  }, FUN.VALUE = numeric(1))  # Pre-allocate to numeric vector
  
  # Return the maximum time
  max(x, na.rm = TRUE)  # Ensure that NA values are ignored
}
