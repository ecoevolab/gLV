#' Load individual SS RDS file
#'
#' This function reads the individual RDS file and returns it as a table.
#' 
#' @param wd Character. The working directory where the \code{Scan/SS_individual.rds} file is located.
#' @param uniqueID A unique identifier for the simulation.
#' 
#' @return A data frame where each row represents a simulation and each column contains different information.
#'
#' @examples
#' # Example usage:
#' wd <- "~/Documents/LAB_ECO/testing"
#' uniqueID <- "1f52dc"
#' table <- load_individual_ss_rds(wd, uniqueID)
#' 
#' @export

load_individual_ss_rds <- function(wd, uniqueID) {
  
  # Create path
  rds_path <- file.path(wd, "Scan/SS_Individual.rds")
  
  # Initialize the join_table
  join_table <- data.frame()
  
  # Read RDS file
  data <- readRDS(rds_path)
  
  # Use length(data) directly
  for (e in seq_along(data)) {
    
    element <- data[[e]]  # Access the e-th element of the list
    
    if (element$ID==uniqueID) {
      Differences_vector <- element$Differences_method
      Rolling_vector <- element$Rolling_variance_method
      
      # Combine all parts into a final data frame
      result_df <- data.frame(
        ID = element$ID,
        Tolerance = element$Tolerance,
        Diff = t(Differences_vector),
        Roll_var = t(Rolling_vector)
      )
      
      join_table <- rbind(join_table, result_df)
    }
  }
  
  return(join_table)
  
}

