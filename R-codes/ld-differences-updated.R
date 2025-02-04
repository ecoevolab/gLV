
#---------------------------------Get Outputs IDs-------------------------------#
## Set the path to your directory
# wd <- "/home/rivera/Cluster/testing"
wd <- "/mnt/atgc-d3/sur/users/mrivera/testing"

outs_path <- file.path(wd, "Outputs")
outs_files <- list.files(outs_path, pattern = "\\.tsv$", full.names = TRUE) # List all TSV files

# The next line is only for testing purposes:
# outs_files <- outs_files[1:20]

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

#---------------------------------------------------------------------------------#
# Initialize an empty data.table to store batch results
library(data.table)
batch_results <- data.frame()

# Process each file in the current batch
# file_paths <- file_paths[20]
for (file in file_paths) {
  
  cat("Getting differences of ", file, "\n")
  
  # Read the TSV file
  output <- fread(file, sep = "\t")
  
  # Log transformation and calculation in one step
  tmp_diff <- diff(abs(colMeans(log(output))))
  
  # Add ID to the vector of the differences
  id <- sub(".*/O_(.*)\\.tsv$", "\\1", file)
  tmp_diff_df <- as.data.frame(t(c(ID = id, round(tmp_diff, 8))))
  
  # Append to batch_results using rbind
  batch_results <- rbind(batch_results, tmp_diff_df, stringsAsFactors = FALSE)
}

#------------------Save the data frame of the differences-----------------------#

# Check if the file exists and merge with old data if it does
if (file.exists(means_ld_path)) {
  old_meansld <- fread(means_ld_path, sep = "\t", header = TRUE, colClasses = c(ID = "character"))
  means_save <- rbindlist(list(old_meansld, batch_results), use.names = TRUE, fill = TRUE)
} else {
  means_save <- batch_results
}

# Save the merged data to the final file
fwrite(means_save, file = means_ld_path, sep = "\t", col.names = TRUE, na = "NA")
cat("Differences completed. The number of differences added were: ", nrow(batch_results), "created\n")




