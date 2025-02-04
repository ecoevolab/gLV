#' Individual Steady State Search by Proportion of Change
#'
#' This function searches for steady states in a simulation by analyzing the \eqn{\log(t / (t+1))} of species abundance across generations.
#'
#' @param uniqueID Character: A unique identifier for the simulation run, typically generated using \code{forge_id}.
#' @param output Matrix: Simulation results where rows represent species and columns represent time steps. Each element in the matrix represents the abundance of a species at a given time step.
#' @param tolerance Numeric: Tolerance value for determining the steady state. The tolerance represents the proportional change between successive generations, calculated using \eqn{\log(t / (t+1))}.
#' @param wd Character: The working directory where the results will be saved. The directory must exist prior to running the function.
#'
#' @details
#' The function computes the logarithmic differences of species abundance to identify stable generations. 
#' It determines a steady state when the absolute difference in logarithmic values between consecutive generations falls 
#' below the specified tolerance. This calculation leverages the logarithmic property \eqn{\log\left(\frac{x}{y}\right) = \log(x) - \log(y)} to facilitate the analysis.
#' 
#' Values of zero in the \code{output} matrix are treated as \code{NA} to avoid 
#' taking the logarithm of zero. If no stable generations are found for a species, the output will reflect this 
#' with an \code{NA} entry for that species.
#'
#' @return A numeric vector indicating the generation where each species reached the steady state. 
#' If a species did not reach a steady state, the corresponding value will be \code{NA}.
#'
#' @examples
#' wd = "~/Documents/LAB_ECO/Simulations"
#' seeds_path <- file.path(wd, "Seeds.tsv" )
#' 
#' # Generate parameters
#' params <- init_data(N_species = 2, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#' 
#' # Run simulation
#' times <- 100  # Define the number of generations
#' output <- run_simulation(N_species = 2, params = params, times = times)
#' 
#' # Generate unique ID
#' uniqueID <- forge_id(wd)
#' 
#' # Example usage:
#' tolerance=0.05
#' individual_prop_SS(uniqueID, output, tolerance, wd)
#'
#' @export

individual_prop_SS <- function(uniqueID, output, tolerance, wd) {
  
  # Ensure required package is loaded
  requireNamespace("utils")
  
  #-------------------------- Declare Data -------------------------------#
  specs <- nrow(output)  # Number of species
  gens <- ncol(output)
  
  # Preallocate stable generations and log_diff matrix
  stable_gen <- numeric(specs)  
  log_diff_df <- matrix(NA, nrow = specs, ncol = gens - 1)  
  
  #---------------------------- Precompute Log Safely ------------------#
  # Replace 0 values in output with NA to avoid log(0) = -Inf
  log_output <- ifelse(output <= 0, NA, log(output) )   # Apply logarithm safely
  
  #---------------------------- Search Stability ------------------------#
  for (s in seq_len(specs)) {
    log_diff <- diff(log_output[s, ])
    log_diff_df[s, ] <- log_diff
    
    # Get stable generations
    stable_points <- which(abs(log_diff) < tolerance)  # Ignore NAs automatically
    stable_gen[s] <- ifelse(length(stable_points) > 0, stable_points[1], NA)
  }
  
  # Naming stable generations
  names(stable_gen) <- paste0("Specie", seq_len(specs))
  
  #------------------------------- Save Differences -----------------------------#
  diff_path <- file.path(wd, "Differences", paste0("ld", uniqueID, ".tsv"))
  
  # Convert log_diff_df to data frame for saving
  log_diff_df <- as.data.frame(log_diff_df)
  colnames(log_diff_df) <- paste0("Diff", seq_len(gens - 1))
  
  # Save log differences
  utils::write.table(log_diff_df, file = diff_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  #------------------------------- Save Column Means Differences -----------------------------#
  mdiff_path <- file.path(wd, "Differences", paste0("means_ld", ".tsv"))
  
  # Calculate column means and format for saving
  means_ld <- as.data.frame(t(colMeans(log_diff_df, na.rm = TRUE)))
  IDs <- uniqueID 
  means_ld <- cbind(IDs, means_ld)
  colnames(means_ld) <- NULL
  #colnames(means_ld) <- paste0("Diff", seq_len(gens - 1)) 
    
  # Function to read and unify rows with different lengths
  read_unified_tsv <- function(file_path, means_ld) {
    
    # Read existing data and find column information
    old_meansld <- read.delim(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE) 
    mld_col <- ncol(means_ld)
    old_col <- ncol(old_meansld)
    old_row <- nrow(old_meansld)
    
    # Determine the number of columns to add
    num_cols_to_add <- max(0, old_col - mld_col)  # Columns to add if means_ld has fewer columns
    na_columns <- data.frame(matrix(NA, nrow = nrow(means_ld), ncol = num_cols_to_add))
    
    # Combine old and new data
    colnames(old_meansld) <- NULL
    colnames(means_ld) <- NULL
    rownames(old_meansld) <- NULL
    rownames(means_ld) <- NULL
    
    # If simulation has more generations than the table add NAs to the table
    if (mld_col == old_col) {
      
      means_save <- data.frame(matrix(NA, nrow = old_row + 1, ncol = old_col))
      
      for (r in 1:(old_row + 1)) {
        means_save[r, ] <- old_meansld[r, ]
        
        if (r == old_row + 1) {
          means_save[r, ] <- means_ld[1, ]
        }
      }
      
    }
    
    # If simulation has less generations than the table add NAs to the output
    if (mld_col < old_col) {
      
      means_save <- data.frame(matrix(NA, nrow = old_row + 1, ncol = old_col))
      
      means_ld <- cbind(means_ld, na_columns)
      
      for (r in 1:(old_row + 1)) {
        means_save[r, ] <- old_meansld[r, ]
        
        if (r == old_row + 1) {
          means_save[r, ] <- means_ld[1, ]
        }
      }
    }
    
    # If simulation has more generations than the table add NAs to the table
    if (mld_col > old_col) {
      
      means_save <- data.frame(matrix(NA, nrow = old_row + 1, ncol = mld_col))
      
      # Create a new data frame with NA values to match mld_col
      padded_rows <- apply(old_meansld, 1, function(row) {
        c(as.character(row), rep(NA, mld_col - length(row) ) )
      })
      old_meansld <- as.data.frame(t(padded_rows), stringsAsFactors = FALSE)
      #old_meansld[, -1] <- as.numeric(old_meansld[, -1])
      
      # Convert non-first columns to numeric while keeping NA values intact
      old_meansld[, -1] <- lapply(old_meansld[, -1], function(col) {
        # Convert to character first to handle possible factors, and then to numeric
        as.numeric(as.character(col))
        # NAs will automatically be preserved in this process
      })
      
      # Fill means_save with old_meansld values
      means_save[1:old_row, ] <- old_meansld
      
      # Assign the last row from means_ld
      means_save[old_row + 1, ] <- means_ld[1, ]
      
      # for (r in 1:old_row + 1) {
      #   if (r == old_row + 1) {
      #     means_save[r, ] <- means_ld[1, ]
      #   }
      #   means_save[r, ] <- old_meansld[r, ]
      #   print(r)
      #   
      # }
    }
    
    return(means_save)
  }
  
  if (file.exists(mdiff_path)) {
    means_save <- read_unified_tsv(file_path = mdiff_path, means_ld)
    
    # Save the combined data (new + old)
    utils::write.table(means_save, file = mdiff_path, sep = "\t", row.names = FALSE, col.names = FALSE)  
    
  } else {
    # Save the data
    utils::write.table(means_ld, file = mdiff_path, sep = "\t", row.names = FALSE, col.names = FALSE)  
    message("file created")
  }
  
  
  #--------------------- Print Messages -----------------------------------------#
  # cat("Log differences steady state search done and saved\n",
  #     "\tWith ID:", uniqueID, "\n",
  #     "\tLog differences saved at path:", diff_path, "\n")
  
  return(stable_gen)
}
