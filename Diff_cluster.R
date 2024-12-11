# This code is stored on the cluster

# Set the path to your directory
functions_path <- "/mnt/atgc-d3/sur/users/mrivera/Functions"
# functions_path <- "/home/rivera/Cluster/Functions/"

# List all R files in the directory and source them
r_files <- list.files(functions_path, pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_files, source))

#---------------------------------Get IDs-------------------------------#
outs_path <- "/mnt/atgc-d3/sur/users/mrivera/Simulations/Outputs/"
# outs_path <- "/home/rivera/Cluster/Simulations/Outputs/"
outs_files <- list.files(outs_path, pattern = "\\.tsv$", full.names = TRUE) # List all TSV files

# Extract unique IDs directly
nids <- sub("O_(.*)\\.tsv", "\\1", basename(outs_files))

#---------------------------------Get IDs that have already been calculated the differences-------------------------------#
wd <- "/mnt/atgc-d3/sur/users/mrivera/Simulations"
# wd <-  "/home/rivera/Cluster/Simulations"

# Load means_ld.tsv
means_ld_path <- file.path(wd, "Differences/means_ld.tsv")
means_ld_ids <- if (file.exists(means_ld_path)) {
  read.delim(means_ld_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 1]
} else {
  character(0)
}

# Items in nids but not in means_ld_ids
result <- setdiff(nids, means_ld_ids)

# Process each file in result
for (res in result) {
  tmp <- file.path(wd, "Outputs", paste0("O_", res, ".tsv"))
  cat("The ID is:", res, "\n")
  
  # Read the TSV file
  output <- read.delim(tmp, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Call the all_prop_SS function if the second column isn't all NA
  if (!all(is.na(output[, 2]))) {
    all_prop_SS(output, tolerance = 0.005, uniqueID = res, wd)
    cat("Differences of the output", res, "completed\n")
  }
}



uniqueID <- "50b1cb"
tmp <- file.path(wd, "Outputs", paste0("O_", uniqueID, ".tsv"))
output <- read.delim(tmp, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
