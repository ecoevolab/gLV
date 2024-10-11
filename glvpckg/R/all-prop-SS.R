#' Steady State Search by Proportional Change in Column Means
#'  
#' @details
#' This function searches for steady states in a simulation by analyzing the changes in the logarithm of the mean simulation results over time. 
#' It leverages the logarithmic property \eqn{\log\left(\frac{x}{y}\right) = \log(x) - \log(y)} to simplify the comparison of successive time steps.
#'
#' @param output Matrix. A matrix of simulation results, where rows represent species and columns represent time steps. Typically generated with \link{run_simulation}.
#' @param tolerance Numeric. The tolerance value used to determine steady states. It compares the differences in consecutive column means over time.
#' @param uniqueID Character. A unique ID string used for saving files. This ID should be generated previously using the \link{forge_ID} function.
#' @param wd Character. The path to the working directory where files will be saved. The file will be saved at: `wd/Differences/means_ld.tsv`.
#'
#' @return Numeric. The first generation where the difference between means is less than the tolerance, or `NA` if no steady state is found.
#' 
#' @examples
#' # Example usage:
#' wd <- "~/Documents/LAB_ECO/Simulations"
#' seeds_path <- file.path(wd, "Seeds.tsv")
#' params <- forge_data(N_species = 2, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#' 
#' # Run simulation for 100 generations
#' output <- run_simulation(N_species = 2, params = params, times = 100)
#'
#' # Search for steady states with a tolerance of 0.05
#' tolerance <- 0.05
#' all_prop_SS(output, tolerance, uniqueID, wd)
#'
#' @export


all_prop_SS <- function(output, tolerance, uniqueID, wd) {
  
  # Compute column means directly
  log_output <- ifelse(output <= 0, NA, log(output))   # Apply logarithm safely
  tmp <- diff(log_output)
  tmp <- as.data.frame(t(colMeans(tmp, na.rm = TRUE)))
  
  # Find the first generation where the difference is below tolerance
  stable_gen <- which(abs(tmp) < tolerance)[1]
  
  # Save Column Means Differences
  mean_diff_path <- file.path(wd, "Differences", paste0("means_ld", ".tsv"))
  means_ld <- cbind(uniqueID, tmp)
  
  # Helper function to unify rows with different lengths
  unify_means <- function(old_meansld, means_ld) {
    mld_col <- ncol(means_ld)
    old_col <- ncol(old_meansld)
    
    # Add NA columns to match sizes
    if (mld_col < old_col) { # New simulation has less generations than older simulations
      means_ld <- cbind(means_ld, matrix(NA, nrow = nrow(means_ld), ncol = old_col - mld_col))
    } else if (mld_col > old_col) { # Simulation has more generations than the older ones
      old_meansld <- cbind(old_meansld, matrix(NA, nrow = nrow(old_meansld), ncol = mld_col - old_col))
    }
    
    # Combine rows
    colnames(old_meansld) <- NULL
    colnames(means_ld) <- NULL
    rownames(old_meansld) <- NULL
    rownames(means_ld) <- NULL
    
    
    means_save <- rbind(as.matrix(old_meansld), as.matrix(means_ld))
    means_save <- as.data.frame(means_save)
    return(means_save)
  }
  
  # Check if the file exists. If so, verify that the columns are the same size. If they are not, add NAs to the table or the 
  # column means to ensure the columns have the same size.
  if (file.exists(mean_diff_path)) {
    old_meansld <- read.delim(mean_diff_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    means_save <- unify_means(old_meansld, means_ld)
  } else {
    means_save <- means_ld
  }
  
  # Save the combined data
  utils::write.table(means_save, file = mean_diff_path, sep = "\t", row.names = FALSE, col.names = FALSE)
  message("file created")
  
  return(stable_gen)
}


