# 11-November-2024
# Manuel Rivera
# The next code is to enable parallelization of the simulations of the cluster, it works as follows:
#   1) It reads the seeds table generated with random.org
#   2) It selects the list of seeds that are available per parameter. This is possible by compairing the file `Seeds_save` and reading the seeds that have been already used for that parameter.
#   3) We then specify the number of simulations per batches and the number of batches to generate the simulations.
#---------------------------------------Declare functions----------------------------------------#
functions_path <- "/mnt/atgc-d3/sur/users/mrivera/Functions"
r_files <- list.files(functions_path, pattern = "\\.R$", full.names = TRUE) # List all R files in the directory
# Source each file without printing output
invisible( lapply(r_files, function(file) {
  capture.output(source(file)) })
)

#-----------------------------------------Forge directories for testing purposes----------------#
wd <- "/mnt/atgc-d3/sur/users/mrivera/testing"
forge_directories(wd)

#---------------------------------------Read seeds table----------------------------------------#
requireNamespace("data.table")
seeds_path <- "/mnt/atgc-d3/sur/users/mrivera/combined_seeds.tsv"
seeds <- as.matrix(data.table::fread(seeds_path, sep = "\t"))

#-----------------------------------Select available seeds--------------------------------------#

# Select seeds that have not been used
params_path <- file.path(wd, "Parameters", "Seeds_save.tsv")

# Check if the params file exists
if (file.exists(params_path)) {
  
  # Read the parameters file and extract relevant columns
  Pms_seeds <- data.table::fread(params_path, sep = "\t", select = c("Population_seed", "Interactions_seed", "Growth_seed"))
  
  # Get possible seeds by subtracting from the loaded seeds
  seeds_list <- list(
    avail_pop_seed = setdiff(seeds, Pms_seeds$Population_seed),
    avail_interacs_seeds = setdiff(seeds, Pms_seeds$Interactions_seed),
    avail_grw_seeds = setdiff(seeds, Pms_seeds$Growth_seed)
  )
  
  cat("The subtraction of seeds was performed to avoid using repetitive seeds.\n")
  
} else {
  
  # Default to using the original seeds if no params file exists
  seeds_list <- list(
    avail_pop_seed = seeds,
    avail_interacs_seeds = seeds,
    avail_grw_seeds = seeds
  )
}


#-----------------------------------Do Parallel simulations--------------------------------#
sims_per_batch <- 400
batches <- 10

# Load required libraries
library(parallel)
library(dplyr)
library(data.table)

# Set the number of cores to use
cores <- detectCores() - 1

for (b in seq_len(batches)) {
  
  # Pre-calculate seeds
  sim_pop_seeds <- sample(seeds_list$avail_pop_seed, sims_per_batch, replace = FALSE)
  sim_int_seeds <- sample(seeds_list$avail_interacs_seeds, sims_per_batch, replace = FALSE)
  sim_grw_seeds <- sample(seeds_list$avail_grw_seeds, sims_per_batch, replace = FALSE)
  
  # Update available seeds
  seeds_list$avail_pop_seed <- setdiff(seeds_list$avail_pop_seed, sim_pop_seeds)
  seeds_list$avail_interacs_seeds <- setdiff(seeds_list$avail_interacs_seeds, sim_int_seeds)
  seeds_list$avail_grw_seed <- setdiff(seeds_list$avail_grw_seed, sim_grw_seeds)
  
  # Parameters and IDs for the batch
  C0_values <- round(runif(sims_per_batch, min = 0, max = 1), 3)
  CN_values <- round(runif(sims_per_batch, min = 0, max = 1), 3)
  N_species_values <- sample(5:100, sims_per_batch, replace = TRUE)
  Diag_val <- -0.5
  uniqueIDs <- sapply(1:sims_per_batch, function(i) forge_id(wd))
  
  # Using mclapply for parallel simulations
  batch_results <- do.call(rbind, mclapply(seq_len(sims_per_batch), function(i) {
    
    # Params for each simulation
    params <- forge_data(N_species_values[i], C0_values[i], CN_values[i], Diag_val, seeds_list)
    output <- run_simulation(N_species_values[i], params = params, times = 700)
    capture.output( output_saver(output, uniqueIDs[i], wd) )
    
    # Create result row
    data.frame(
      ID_simulation = uniqueIDs[i],
      N_specs = N_species_values[i],
      Prob_0 = C0_values[i],
      Prob_neg = CN_values[i],
      Diagonal = Diag_val,
      Population_seed = params$Seeds[1],
      Interactions_seed = params$Seeds[2],
      Growth_seed = params$Seeds[3]
    )
  }, mc.cores = cores))
  
  # Save batch results
  table_path <- file.path(wd, "Parameters", paste("Pbatch_0", b, ".tsv", sep = ""))
  fwrite(batch_results, file = table_path, sep = "\t")
}

#-----------------------------Combine tables into one-----------------------------#

# Define the path where your files are stored
path <- file.path(wd, "Parameters") # Replace with your directory path

# Get a list of all files matching the pattern in the specified directory
file_list <- list.files(path = path, pattern = "^Pbatch_\\d+\\.tsv$", full.names = TRUE)

# Read and combine all files into one data frame
combined_table <- rbindlist(lapply(file_list, fread))

#------------------------------Combine table with the old one---------------------#
Ps_path <- file.path(wd, "Parameters", "Master_table.tsv")

if (file.exists(Ps_path)) {
  # Read the existing table
  old_params <- data.table::fread(Ps_path)
  
  # Append new data using rbindlist (more efficient than rbind)
  final_table <- rbindlist(list(old_params, combined_table))
  
  # Directly save the all_batches table
  fwrite(final_table, file = Ps_path, sep = "\t")
  message("Seeds_save file UPDATED")
  
} else {
  # Directly save the all_batches table
  fwrite(combined_table, file = Ps_path, sep = "\t")
  message("No Seeds_save file, creating a new one...")
}

#----------------------------Delete old tables--------------------------------------#

for (file in file_list) {
  file.remove(file)
  message(file, " removed... \n")
}
