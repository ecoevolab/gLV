
library("tictoc", lib = "/mnt/atgc-d3/sur/modules/pkgs/tidyverse_mrc")
tictoc::tic("Section 0: Total running time:")

# ==== Load data and declare functions ====

# Load Parameters table
params_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Exp03-D25M02.tsv")

# Source function to regenerate parameters
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/Forge-gLV-Parameters.R")

# ==== ODE solver ====
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

tictoc::tic("Section 1: Divide data into chunks")
num_cores <- parallel::detectCores() - 1  # Use one less than the total number of cores
cat("The number of cores that will be used are: ", num_cores, "\n")

split_table <- function(df, n_chunks) {
  split(df, cut(seq_len(nrow(df)), breaks = n_chunks, labels = FALSE))
}

chunks <- split_table(params_table, num_cores)
toc() # For section 1

# ==== Generate directories for each core ====
tictoc::tic("Section 2: Generate directories for each core")

# Generate workers directories
ode_dir <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Exp03-D25M02/Simulate_ODE"

# Function to create main and worker directories
create_dirs <- function(main_dir, num_cores) {
  if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)
  
  # Create worker directories
  worker_dirs <- file.path(main_dir, paste0("worker_", seq_len(num_cores)))
  invisible(lapply(worker_dirs, dir.create, showWarnings = FALSE))
  
  return(worker_dirs)
}

# Create directories 
workers_ODE <- create_dirs(ode_dir, num_cores)
toc() # For section 2

# ==== Wrapper for running all required steps ====
parallel.sims <- function(index, path_ODE) {
  
  # Generate parameters
  params <- regenerate(index)
  
  # Run simulation
  output_ode <- ode.simulate(times = 700, params)
  
  # Define paths
  id <- index$id
  save_ode <- file.path(path_ODE, paste0("O_", id, ".tsv"))
  
  # Calculate NAs
  Total.NAs <- sum(is.na(output_ode))
  
  # Save simulation
  utils::write.table(output_ode, file = save_ode, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  return(c(id = id,
           Total.NAs = Total.NAs)
  )
}

# ==== Parallelize it ====
tictoc::tic("Section 3: Run simulations using the parallel package")

NAs_vecs <- parallel::mclapply(1:num_cores, function(core_id) {
  
  message("Starting worker ", core_id, "....\n")
  
  core_chunk <- chunks[[core_id]]  # rows assigned to this core
  
  # cat("\nODE path", workers_ODE[core_id], "\n")
  
  na.vec <- lapply(1:nrow(core_chunk), function(i) {
    parallel.sims(core_chunk[i, ], 
            path_ODE = workers_ODE[core_id])
  })
  
  message("Ending worker ", core_id, "....\n")
  
  return(na.vec)
  
}, mc.cores = num_cores)
toc() # For section 3

# ==== Get NAs number on simulations====
tictoc::tic("Section 4: Count total number of NAs")
NAs.counts <- unlist(NAs_vecs)
counts_df <- as.data.frame(matrix(NAs.counts, ncol = 2, byrow = TRUE))
colnames(counts_df) <- unique(names(NAs.counts))

# Save Parameters as TSV
save.path = "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Exp03-D25M02/Nas-counting.tsv"
data.table::fwrite(x = counts_df, file = save.path, sep = "\t")
toc() # For section 4

# ==== Create symbolic links====
tictoc::tic("Section 5: Generate symbolic links")

source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/forge_symlinks.R")
# Define source and target directories
source_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Exp03-D25M02/Simulate_ODE"
target_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Exp03-D25M02/Simulate_ODE/Unified"
generate_symlinks(source_path = source_path, target_path = target_path)

toc() # For section 5


toc() # For Total running time