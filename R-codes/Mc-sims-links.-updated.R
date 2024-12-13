#--------------------------------Specify directories--------------------------------#

# Define source and target directories
source_parent <- "/mnt/atgc-d3/sur/users/mrivera/tmp-times/C40"
# source_parent <- "/home/rivera/Cluster/tmp-times/C40"

target_parent <- "/mnt/atgc-d3/sur/users/mrivera/tmp-times/unify"
# target_parent <- "/home/rivera/Cluster/tmp-times/unify"

cat("The source directory is:", source_parent, "\n")

cat("The target directory is:", target_parent, "\n")

# Get all files in the worker directory recursively
all_files <- list.files(source_parent, recursive = TRUE, full.names = TRUE)

# Filter files matching the pattern Outputs/O_*.tsv
outs_files <- grep("Outputs/O_.*\\.tsv$", all_files, value = TRUE)
params_files <- grep("Parameters/Seeds_save.tsv$", all_files, value = TRUE)

#--------------------------------Create outputs symbolic links--------------------------------#
# Generate outputs symbolic links
for (file in outs_files) {
    
    # Extract relative path
    relative_path <- sub(paste0("^", source_parent, "/core_\\d+/"), "", file)
    
    # Construct the target path
    target_path <- file.path(target_parent, relative_path)
    
    # Ensure the target directory exists
    target_dir <- dirname(target_path)
    if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
    
    # Create symbolic link
    file.symlink(file, target_path)
}
  
cat("The symbolic links of output was completed...\n")

#--------------------------------Create parameters symbolic links--------------------------------#
# Generate Paramaeters symbolic links
library(data.table)  # Ensure data.table package is loaded for fread and fwrite

for (file in params_files) {
  
  # Extract relative path
  relative_path <- sub(paste0("^", source_parent, "/core_\\d+/"), "", file)
  
  # Construct the target path
  target_path <- file.path(target_parent, relative_path)
  
  # Load the data from the source file
  data <- fread(file)  # fread is more efficient than read.csv for large files
  
  # Ensure the target directory exists
  target_dir <- dirname(target_path)
  if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
  
  # Check if the target file exists
  if (file.exists(target_path)) {
    # Append data to the existing file (without writing header)
    fwrite(data, target_path, append = TRUE, col.names = FALSE, sep = "\t")
  } else {
    # Write the data to a new file (with header)
    fwrite(data, target_path, sep = "\t")
  }
}

cat("The master table of parameters was completed...\n")
