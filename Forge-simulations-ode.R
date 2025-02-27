
# ==== Load data and declare functions ====

# Load Parameters table
params_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Exp03-D25M02.tsv")

# Source function to regenerate parameters
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/Forge-gLV-Parameters.R")

# ODE solver
ode.simulate <- function(times, params) {
  
  # Define the equation
  glv_model <- function(t, x0, params) {
    r <- params$mu         # Growth rate vector
    A <- params$M          # Interaction matrix
    
    # Compute dx/dt for each species
    dx <- x0 * (r + A %*% x0)
    list(dx)
  }
  
  time_seq <- seq(1, times, by = 1)  # Define the time sequence
  
  # Get solution
  results <- tryCatch(
    R.utils::withTimeout(deSolve::ode(y = params$x0, times = time_seq, func = glv_model, 
                                      parms = params,
                                      method = "ode45",
                                      rtol = 1e-06, 
                                      atol = 1e-06),
                         timeout = 600), 
    error = function(e) {
      message(">> Simulation failed... skipping")
      return(NULL) # Return NULL instead of an NA matrix
    })
  
  
  # Remove the first column (`time`) and transpose it so columns represent generations
  # Check for valid output and return transposed results
  if (!is.null(results) && ncol(results) > 1) {
    return(t(results[, -1]))
  } else {
    return(matrix(NA, nrow = nrow(params$M), ncol = times)) # Return properly shaped NA matrix
  }
}

# ==== Divide data into chunks====

library(parallel)

num_cores <- detectCores() - 1  # Use one less than the total number of cores
cat("The number of cores that will be used are: ", num_cores, "\n")

split_table <- function(df, n_chunks) {
  split(df, cut(seq_len(nrow(df)), breaks = n_chunks, labels = FALSE))
}

chunks <- split_table(params_table, num_cores)

# ==== Generate directories for each core ====
# 
# Generate workers directories
# mias_dir <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Exp03-D25M02/Simulate_miaSim"
ode_dir <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Exp03-D25M02/Simulate_ODE"

# Function to create main and worker directories
create_dirs <- function(main_dir, num_cores) {
  if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)
  
  # Create worker directories
  worker_dirs <- file.path(main_dir, paste0("worker_", seq_len(num_cores)))
  invisible(lapply(worker_dirs, dir.create, showWarnings = FALSE))
  
  return(worker_dirs)
}

# Create directories for both simulations
# workers_mia <- create_dirs(mias_dir, num_cores)
workers_ODE <- create_dirs(ode_dir, num_cores)

# ==== Data Import ====
parallel.sims <- function(index, path_ODE) {
  
  # Generate parameters
  params <- regenerate(index)
  
  # Run simulation
  output_ode <- ode.simulate(times = 700, params)
  # output_mia <- sim_glv(params = params, n_t = 700)
  
  # Define paths
  id <- index$id
  save_ode <- file.path(path_ODE, paste0("O_", id, ".tsv"))
  # save_mia <- file.path(path_miaSim, paste0("O_", id, ".tsv"))
  
  # Calculate NAs
  NA.ODE <- sum(is.na(output_ode))
  #NA.mia <- sum(is.na(output_mia))
  
  # Save simulation
  utils::write.table(output_ode, file = save_ode, sep = "\t", row.names = FALSE, col.names = TRUE)
  #utils::write.table(output_mia, file = save_mia, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  return(c(id = id,
           NA.ODE = NA.ODE
           #NA.mia = NA.mia,
  )
  )
}
