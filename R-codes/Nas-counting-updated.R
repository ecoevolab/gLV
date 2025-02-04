
# Define the root path
# root_path <- "/home/rivera/Cluster/testing/Tolerances"
root_path <- "/mnt/atgc-d3/sur/users/mrivera/testing/Tolerances"

# Get all combinations of tolerances
tol_r <- list.files(root_path, pattern = "^tol_r", full.names = TRUE)
tol_a <- list.files(tol_r[1], pattern = "^tol_a", full.names = FALSE)

# Create a function to calculate NAs
calculate_nas <- function(file_path) {
  data <- data.table::fread(file_path)
  tmp <- sum(is.na(data))
  return(tmp)
}

# Initialize an empty list to store results
na_counts_list <- list()

# Loop over tol_r and tol_a to calculate NA counts
# The next code iterates over each tolerance value and checks all IDs whitin the path
for (r in tol_r) {
  for (a in tol_a) {
    # Get the directory path for the current tolerance combination
    dir_path <- file.path(r, a)
    
    # List .tsv files 
    tsv_files <- list.files(dir_path, pattern = "\\.tsv$", full.names = TRUE)
    # print(tsv_files)
    
    # Calculate NAs for each file
    # cat("The file is: ", tsv_files[1], "\n")
    na_counts <- sapply(X = tsv_files, FUN = calculate_nas)
    # cat("The number of NAs are: ", na_counts, "\n")
    
    # Create a data frame for this combination of tolerances
    combination_name <- paste(basename(r), basename(a), sep = "-")
    cat("Combination name: ", combination_name, "completed...\n")
    
    # Extract IDs
    ids <- sub(".*id_(.*)\\.tsv$", "\\1", tsv_files)
    # print(ids)
    
    # Create list of the combinations
    na_counts_list[[combination_name]] <- data.frame(
      TSV_ID = ids,
      NA_Counts = na_counts,
      stringsAsFactors = FALSE  # Avoid converting strings to factors
    )
    names(na_counts_list[[combination_name]])[2] <- combination_name
  }
}


# Combine all the results into a single data frame
library(dplyr)

# Combine all data frames
final_df <- Reduce(function(x, y) merge(x, y, by = "TSV_ID", all = TRUE), na_counts_list)

# Optionally, save the final data frame to a CSV
save_path <- "/mnt/atgc-d3/sur/users/mrivera/testing/NA_counts_test.tsv"
data.table::fwrite(final_df, file = save_path, sep = "\t")

