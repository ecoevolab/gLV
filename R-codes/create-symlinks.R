#--------------------------------Specify directories--------------------------------#

# Define source and target directories
source_parent <- "/mnt/atgc-d3/sur/users/mrivera/JAN-test/workers"

target_parent <- "/mnt/atgc-d3/sur/users/mrivera/JAN-test/unified"

cat("The source directory is:", source_parent, "\n")

cat("The target directory is:", target_parent, "\n")

# Get all files in the worker directory recursively
all_files <- list.files(source_parent, recursive = TRUE, full.names = TRUE)

# Filter files matching the pattern Outputs/O_*.tsv
outs_files <- grep("Outputs/O_.*\\.tsv$", all_files, value = TRUE)
params_files <- grep("Parameters/Seeds_save.tsv$", all_files, value = TRUE)

#--------------------------------Create outputs symbolic links--------------------#
# Extract relative paths for all files at once
relative_paths <- sub(paste0("^", source_parent, "/core_\\d+/"), "", outs_files)

# Construct target paths
target_paths <- file.path(target_parent, relative_paths)

# Get unique target directories
target_dirs <- unique(dirname(target_paths))

# Create target directories only if they don't exist
for (dir in target_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Create symbolic links in a vectorized way
mapply(file.symlink, outs_files, target_paths)

cat(sprintf("The symbolic links for %d output files were completed successfully.\n", length(outs_files)))

#--------------------------------Create parameters symbolic links-------------------#
library(data.table)

# Combine all TSVs into a single data.table
master_table <- rbindlist(lapply(params_files, function(file) {
  tryCatch({
    fread(file)  # Read the TSV file
  }, error = function(e) {
    cat(sprintf("Error reading file %s: %s\n", file, e$message))
    return(NULL)  # Skip problematic files
  })
}), use.names = TRUE, fill = TRUE)

# Specify the output path for the master table
master_table_path <- file.path(target_parent, "Parameters", "master_table.tsv")
master_table_dir <- dirname(master_table_path)
if (!dir.exists(master_table_dir)) dir.create(master_table_dir, recursive = TRUE)

# Save the combined master table to a TSV file
fwrite(master_table, master_table_path, sep = "\t")

cat(sprintf("Master table has been successfully saved to %s.\n", master_table_path))



