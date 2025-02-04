
cat(
  paste0(rep("=", 20), collapse = ""), "  Running code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)

#---------------------------------------Declare functions----------------------------------------#
functions_path <- "/mnt/atgc-d3/sur/users/mrivera/Functions"
r_files <- list.files(functions_path, pattern = "\\.R$", full.names = TRUE) # List all R files in the directory
# Source each file without printing output
invisible( lapply(r_files, function(file) {
  capture.output(source(file)) })
)

#-----------------------------------------Forge directories for testing purposes----------------#
wd <- "/mnt/atgc-d3/sur/users/mrivera/JAN-test"

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
  
  cat("No seeds to substract. \n")
}


#-------------------------------------------Separate IDs by core------------------------------#

library(parallel)

n_cores <- 20  # Use one less than the total number of cores
cat("The number of cores that will be used are: ", n_cores, "\n")

#----------------------------Generate parameters-------------------------------#
n_sims <- 5000
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


#-----------------------------Generalized Lotka-Volterra equation---------------#

requireNamespace("deSolve")

ode_function <- function (times, params) {
  
  # Define the equation
  glv_model <- function(t, x, params) {
    r <- params$Growths         # Growth rate vector
    A <- params$Interactions          # Interaction matrix
    
    # Compute dx/dt for each species
    dx <- x * (r + A %*% x)
    list(dx)
  }
  
  time_seq <- seq(0, times, by = 1)  # Define the time sequence
  
  # Get solution
  results <- deSolve::ode(y = params$Population, times = time_seq, func = glv_model, parms = params, method = "ode45",
                          rtol = 1e-06, 
                          atol = 1e-06)
  
  # Remove the column of times and get the transversal, where rows are species and column generations
  ode_df <- as.matrix(t(results)) [-1, ] 
  return(ode_df)
}

#----------------------------Parallel Processing-------------------------------#
cores_path <- file.path(wd, "workers")
dir.create(cores_path)

# Create core-specific directories
core_dirs <- sapply(1:n_cores, function(core_id) {
  dir_name <- file.path(cores_path, paste0("core_", core_id))
  if (!dir.exists(dir_name)) dir.create(dir_name)
  return(dir_name)
})
 
lapply(core_dirs, forge_directories)
 
# Function for running simulations on each core
simulate_batch <- function(batch, core_id) {
  core_dir <- core_dirs[core_id]
  for (i in batch) {
  
    # Generate parameters
    params <- forge_data(N_species_values[i], C0_values[i], CN_values[i], Diag_val, seeds_list)
    
    # Run simulation
    output <- ode_function(times = 700, params)
    
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


cat(
  paste0(rep("=", 20), collapse = ""), "  Ending code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)
