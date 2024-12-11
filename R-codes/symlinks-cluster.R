#-----------------------------------Get valid simulations------------------------------#
# Set working directory
# wd <- "/home/rivera/Cluster/testing"
wd <- "/mnt/atgc-d3/sur/users/mrivera/testing"

library(data.table)

# Filter simulations
diff_path <- file.path(wd, "Differences", "means_ld.tsv" )
table <- fread(diff_path, sep = "\t", header = TRUE, colClasses = c(ID = "character"))

# Use `complete.cases()` to filter rows directly and get valid IDs
controls <- table[complete.cases(table), ]

# Get the simulations we will repeat
raw_reps_ids <- setdiff(table$ID, controls$ID)
control_ids <- controls$ID

# Get IDs of the simulations that have been already done 
# rds_path <- file.path(wd, "Rp_tol", "ids-completed.rds")
# loaded_vector <- readRDS("my_vector.rds")

#--------------------------------Create symbolic links--------------------------------#
# 
# # Define source and target directories
# source_parent <- "/mnt/atgc-d3/sur/users/mrivera/testing/test02"
# # source_parent <- "/home/rivera/Cluster/testing/test02"
# 
# target_parent <- "/mnt/atgc-d3/sur/users/mrivera/testing/Tolerances"
# # target_parent <- "/home/rivera/Cluster/testing/Tolerances"
# 
# library(parallel)
# 
# # Get all worker directories
# worker_dirs <- list.dirs(source_parent, recursive = FALSE, full.names = TRUE)
# 
# # Define a function to process each worker directory
# process_worker <- function(worker) {
#   # Get all files in the worker directory recursively
#   files <- list.files(worker, recursive = TRUE, full.names = TRUE)
#   
#   # Iterate over each file
#   for (file_path in files) {
#     # Extract relative path
#     relative_path <- sub(paste0("^", source_parent, "/worker_\\d+/"), "", file_path)
#     
#     # Construct the target path
#     target_path <- file.path(target_parent, relative_path)
#     
#     # Ensure the target directory exists
#     target_dir <- dirname(target_path)
#     if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
#     
#     # Create symbolic link
#     file.symlink(file_path, target_path)
#   }
#   
#   msg <- cat("The symbolic links of worker ", worker, " were completed...\n")
#   return(msg)
# }
# 
# num_cores <- detectCores() - 1  # Use one less than the total number of cores
# cat("The number of cores that will be used are: ", num_cores, "\n")
# 
# # Use mclapply for parallel processing (for Unix-like systems)
# mclapply(worker_dirs, process_worker, mc.cores = num_cores)


#-----------------------------------Testing---------------------------------------------#

path <- "/mnt/atgc-d3/sur/users/mrivera/testing/Tolerances/tol_r1e-06/tol_a1e-06"
# path <- "/home/rivera/Cluster/testing/Tolerances/tol_r1e-06/tol_a1e-06"

files <- list.files(path, recursive = TRUE, full.names = FALSE)
clean_files <- sub("id_(\\w+)\\.tsv", "\\1", files)

print(clean_files[1])


tmp <- setdiff(table$ID, clean_files)

cat("--------------------")
cat("\n The number of IDs on the table is of: ", length(table$ID))
cat("\n The number of IDs that should be left are: ", length(table$ID) - length(raw_reps_ids))
cat("\n\n The number of IDs on the table after subset is of: ", length(tmp), "\n")


