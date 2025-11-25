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
find_stability <- function(output, threshold = 0.05, min_length = 10) {
  
  # Apply the logic across each row of output
  x <-  apply(output, 1, function(row) {
    test <- abs(diff(row)) <= threshold
      rle_res <- rle(test)
      valid_idx <- which(rle_res$values & rle_res$lengths >= min_length)[1]
      
      if (is.na(valid_idx)) return(NA_real_)
      if (valid_idx == 1) return(1)  # Edge case: stable from start
      
      sum(rle_res$lengths[seq_len(valid_idx - 1)]) + 1
  })
  
  # Return the maximum time
  max(x, na.rm = TRUE)  # Ensure that NA values are ignored
}
