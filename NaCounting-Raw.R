cat(
  paste0(rep("=", 20), collapse = ""), "  Running code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)

# Define the root path
root_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Unified/RawUnified-D10M02Y24"


rtol_values <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1) # Relative tolerance
atol_values <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1) # Absolute tolerance
grid <- expand.grid(rtol = rtol_values, atol = atol_values)

# Create a function to calculate NAs
calculate_nas <- function(file_path) {
  data <- data.table::fread(file_path)
  tmp <- sum(is.na(data))
  return(tmp)
}

na_counts_list <- lapply(1:nrow(grid), function(i) {
  
  tol <- grid[i, ]
  a <- format(tol["atol"], scientific = TRUE, digits = 3)
  r <- format(tol["rtol"], scientific = TRUE, digits = 3)
  
  dir_path <- file.path(root_path, paste0("tol_r", r), paste0("tol_a", a))
  
  tsv_files <- list.files(dir_path, pattern = "\\.tsv$", full.names = TRUE)
  
  ids <- sub(".*_(.*?)\\.tsv$", "\\1", basename(tsv_files))  # Extract IDs directly
  
  na_counts <- sapply(X = tsv_files, FUN = calculate_nas)
  
  combination_name <- paste(r, a, sep = "-")
  
  df <- data.frame(
    TSV_ID = ids,
    NA_Counts = na_counts,
    stringsAsFactors = FALSE
  )
  names(df)[2] <- combination_name
  
  df  # Return data frame
})

# Assign names correctly
names(na_counts_list) <- apply(grid, 1, function(tol) paste(tol["rtol"], tol["atol"], sep = "-"))


# Combine all the results into a single data frame
library(dplyr)

# Combine all data frames
final_df <- Reduce(function(x, y) merge(x, y, by = "TSV_ID", all = TRUE), na_counts_list)

# Optionally, save the final data frame to a CSV
save_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Unified/NA-RawUnified-D10M02Y24.tsv"
data.table::fwrite(final_df, file = save_path, sep = "\t")

cat(
  paste0(rep("=", 20), collapse = ""), "  Ending code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)


