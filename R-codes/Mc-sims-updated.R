#---------------------------------------Declare functions----------------------------------------#
functions_path <- "/mnt/atgc-d3/sur/users/mrivera/Functions"
r_files <- list.files(functions_path, pattern = "\\.R$", full.names = TRUE) # List all R files in the directory
# Source each file without printing output
invisible( lapply(r_files, function(file) {
  capture.output(source(file)) })
)

#-----------------------------------------Forge directories for testing purposes----------------#
wd <- "/mnt/atgc-d3/sur/users/mrivera/tmp-times/C40"

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


#-------------------------------------------Separate IDs by core------------------------------#

library(parallel)

n_cores <- 40  # Use one less than the total number of cores
cat("The number of cores that will be used are: ", n_cores, "\n")

#----------------------------Generate parameters-------------------------------#
# Path to the main directory
main_dir <- "/mnt/atgc-d3/sur/users/mrivera/tmp-times/C40"

n_sims <- 1000
sim_batches <- split(1:n_sims, rep(1:n_cores, length.out = n_sims)) # Split simulations into batches

# Pre-calculate seeds
sim_pop_seeds <- sample(seeds_list$avail_pop_seed, n_sims, replace = FALSE)
sim_int_seeds <- sample(seeds_list$avail_interacs_seeds, n_sims, replace = FALSE)
sim_grw_seeds <- sample(seeds_list$avail_grw_seeds, n_sims, replace = FALSE)

# Update available seeds
seeds_list$avail_pop_seed <- setdiff(seeds_list$avail_pop_seed, sim_pop_seeds)
seeds_list$avail_interacs_seeds <- setdiff(seeds_list$avail_interacs_seeds, sim_int_seeds)
seeds_list$avail_grw_seed <- setdiff(seeds_list$avail_grw_seed, sim_grw_seeds)

# Parameters and IDs for the batch
C0_values <- round(runif(n_sims, min = 0, max = 1), 3)
CN_values <- round(runif(n_sims, min = 0, max = 1), 3)
N_species_values <- sample(5:100, n_sims, replace = TRUE)
Diag_val <- -0.5
uniqueIDs <- sapply(1:n_sims, function(i) forge_id(wd))

#----------------------------Parallel Processing-------------------------------#

# Create core-specific directories
core_dirs <- sapply(1:n_cores, function(core_id) {
  dir_name <- file.path(wd, paste0("core_", core_id))
  if (!dir.exists(dir_name)) dir.create(dir_name)
  return(dir_name)
})

# Function for running simulations on each core
simulate_batch <- function(batch, core_id) {
  core_dir <- core_dirs[core_id]
  for (i in batch) {
    # Check if i matches 1, 10, 100, or 1000
    if (i %in% c(1, 10, 100, 1000)) {
      readable_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      cat("Core", core_id, "- Time when i =", i, "is:", readable_time, "\n")
    }
    
    forge_directories(core_dir)
    
    # Generate simulations
    params <- forge_data(N_species_values[i], C0_values[i], CN_values[i], Diag_val, seeds_list)
    output <- run_simulation(N_species_values[i], params, 700)
    
    # Save output in the core's directory
    output_saver(output, uniqueIDs[i], core_dir)
    
    # Save parameters in the core's directory
    params_seed_saver(N_species_values[i], C0_values[i], CN_values[i], Diag_val, params, uniqueIDs[i], core_dir)
  }
}

# Run simulations in parallel
mclapply(1:n_cores, function(core_id) {
  simulate_batch(sim_batches[[core_id]], core_id)
}, mc.cores = n_cores)

