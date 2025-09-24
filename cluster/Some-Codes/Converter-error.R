library(data.table)
library(arrow)
library(parallel)
library(tools)

# Function to convert a single TSV to feather
converter <- function(file) {
  tryCatch({
    df <- fread(file, sep = "\t") # Read the TSV file
    feather_file <- paste0(file_path_sans_ext(file), ".feather")# Generate feather name
    write_feather(df, feather_file) # Write to feather
    file.remove(file) # Remove file
    # Return success message
    cat(paste("Converted: ", file, "/n"))
  }, error = function(e) {
    cat(paste("Error converting", file, ":", e$message))
  })
}

# Main function with parallel processing
convert_all_tsvs <- function(base_path, pattern = "^O_.*\\.tsv$", num_cores = parallel::detectCores() - 1) {
  # Find all work directories
  work_dirs <- list.files(base_path, recursive = FALSE, full.names = TRUE)
  
  # Find all TSV files across all work directories
  lapply(work_dirs, function(dir) {
    all_files <- list.files(dir, recursive = FALSE, full.names = TRUE, pattern = pattern)
    message("Found ", length(all_files), " TSV files to convert...\n")
    mclapply(all_files, converter, mc.cores = num_cores)
  })
  
 message("Conversion complete!\n")
}

# Example usage:
path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Exp06-D29-Apr/mc-apply"
convert_all_tsvs(path)

#=====================  Testing conversion ============== 
work_dirs <- list.files(path, recursive = FALSE, full.names = TRUE)

x = lapply(work_dirs, function(dir){
    message("Starting directory ", dir, " ...\n")
    files <- list.files(dir, recursive = FALSE, full.names = TRUE, pattern = "^O_.*\\.tsv$")
    cat("Missing conversion ", length(files), " ...\n")
    return(length(files))
})


#===================== Rename Info files ============== 
lapply(work_dirs, function(dir){
    message("Starting directory ", dir, " ...\n")
    files <- list.files(dir, recursive = FALSE, full.names = TRUE, pattern =  "^E_.*-Info\\.tsv$")
    cat("Missing correct naming ", length(files), " ...\n")

    # name <- tools::file_path_sans_ext(files)
    # new_files <- paste0(name, ".feather")  
    # file.rename(files, new_files) # Fix naming
})

#=====================  Generate new symlinks ============== 
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/src/forge-symls.R") # generate symbolic links
exp_dir <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Exp06-D29-Apr"
purrr::walk(work_dirs, ~gen_syml(.x, exp_dir))

