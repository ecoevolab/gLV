# This code is stored on the cluster

# Set the path to your directory
#wd <- "/home/rivera/Cluster/testing"
wd <- "/mnt/atgc-d3/sur/users/mrivera/testing"

#---------------------------------Get Outputs IDs-------------------------------#
outs_path <- file.path(wd, "Outputs")
outs_files <- list.files(outs_path, pattern = "\\.tsv$", full.names = TRUE) # List all TSV files

# The next line is only for testing purposes:
#outs_files <- outs_files[1:20]

# Extract unique IDs directly
nids <- sub("O_(.*)\\.tsv", "\\1", basename(outs_files))

#-------------Get IDs that have already been calculated the differences---------#

# Load means_ld.tsv
means_ld_path <- file.path(wd, "Differences", "means_ld.tsv")
means_ld_ids <- if (file.exists(means_ld_path)) {
  read.delim(means_ld_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 1]
} else {
  character(0)
}

# Items in nids but not in means_ld_ids
result <- setdiff(nids, means_ld_ids)

# Vectorize operations by pre-computing file paths
file_paths <- file.path(wd, "Outputs", paste0("O_", result, ".tsv"))

#-----------------------------------------Calculate differences---------------------#
#
## Load required libraries
library(parallel)
library(dplyr)
library(data.table)

# Function to split file paths into batches
split_into_batches <- function(paths, batch_size) {
  split(paths, ceiling(seq_along(paths) / batch_size))
}

# Number of cores to use
num_cores <- detectCores() - 1


# Split file_paths into batches
file_batches <- split_into_batches(file_paths, batch_size = 1000)

# Process each batch in parallel
mclapply(seq_along(file_batches), function(batch_idx) {
  
  # Get the current batch
  current_batch <- file_batches[[batch_idx]]

  # Initialize a data.table to store results for this batch
  batch_results <- data.table()

  # Process each file in the current batch
  for (res_path in current_batch) {

    cat("Getting differences of ", res_path, "\n")

    # Read the TSV file
    output <- fread(res_path, sep = "\t")

    # Log transformation and calculation in one step
    tmp_diff <- diff(abs(colMeans(log(output))))

    # Add ID to the vector of the differences
    id <- sub(".*/O_(.*)\\.tsv$", "\\1", res_path)
    tmp_diff <- as.data.frame(t(c(ID = id, round(tmp_diff, 8), row.names = NULL)))

    # Convert to data.table and append to batch results
    batch_results <- rbindlist(list(batch_results, as.data.table(tmp_diff)), use.names = FALSE, fill = FALSE)
  }
  
  # Save the current batch to a file
  batch_file_name <- file.path(wd, "Differences", sprintf("Par_diff_%02d.tsv", batch_idx))   # e.g., Differences_01.tsv
  cat(batch_file_name)
  fwrite(batch_results, file = batch_file_name, sep = "\t", col.names = TRUE, na = "NA")
  message(paste("Batch file", batch_file_name, "created"))
  
}, mc.cores = num_cores)  # Number of parallel processes

#------------------------------Merge data---------------------------#

# List all the batch files from the 'Differences' folder
batch_files <- list.files(file.path(wd, "Differences"), pattern = "Par_diff_\\d{2}\\.tsv", full.names = TRUE)

# Read all batch files into a list and merge them at once
final_results <- rbindlist(lapply(batch_files, function(file) {
  fread(file, sep = "\t", header = TRUE, colClasses = c(ID = "character"))
}), use.names = FALSE, fill = FALSE)

# Optionally, remove all batch files after merging if no longer needed
unlink(batch_files)

# Check if the file exists and merge with old data if it does
if (file.exists(means_ld_path)) {
  old_meansld <- fread(means_ld_path, sep = "\t", header = TRUE, colClasses = c(ID = "character"))
  means_save <- rbindlist(list(old_meansld, final_results), use.names = TRUE, fill = TRUE)
} else {
  means_save <- final_results
}

# Save the merged data to the final file
fwrite(means_save, file = means_ld_path, sep = "\t", col.names = TRUE, na = "NA")
message(paste("Final merged file", means_ld_path, "created"))

