cat(
  paste0(rep("=", 20), collapse = ""), "  Running code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)

#' Load Parameters table
params_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Data-D25M02.tsv")

#' Add function to generate the parameters
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/Forge-gLV-Parameters.R")

#' The next function is for solvide gLV model with ode45
requireNamespace("deSolve")

ode_function <- function (times, params) {
  
  # Define the equation
  glv_model <- function(t, x, params) {
    r <- params$mu         # Growth rate vector
    A <- params$M          # Interaction matrix
    
    # Compute dx/dt for each species
    dx <- x * (r + A %*% x)
    list(dx)
  }
  
  time_seq <- seq(0, times, by = 1)  # Define the time sequence
  
  # Get solution
  results <- deSolve::ode(y = params$x0, times = time_seq, func = glv_model, parms = params,
                          method = "ode45",
                          rtol = 1e-06, 
                          atol = 1e-06)
}

#' The next function is for solving gLV equation with miaSim
sim_glv <- function(params = params, n_t = n_t){
  
  # Check that matrrix is square
  if (nrow(params$M) != ncol(params$M)) {
    stop("Matrix is not n*n", call. = TRUE)
  }
  
  
  # Use miaSim to simulate standard gLV
  # Added timeout for dealing with rare instance where simulation
  # keeps going forever. An issue is that we cannto distinguish failure
  # by timeout from other types of failure.
  msim <- tryCatch(
    R.utils::withTimeout(miaSim::simulateGLV(n_species = nrow(params$M), 
                                             names_species = names(params$x0),
                                             A = params$M,
                                             x0 = params$x0,
                                             growth_rates = params$mu,
                                             sigma_migration = 0,
                                             epoch_p = 0,
                                             t_external_events = NULL,
                                             t_external_durations = NULL,
                                             stochastic = FALSE,
                                             migration_p = 0,
                                             error_variance = 0,
                                             norm = FALSE,
                                             t_end = n_t),
                         timeout = 600), 
    error = function(e) {
      message(">> Simulation failed... skipping")
      return(NULL) # Return NULL instead of an NA matrix
    })
  
  # Ensure msim is valid before extracting counts
  if (inherits(msim, "TreeSummarizedExperiment")) {
    counts <- SummarizedExperiment::assay(msim, "counts")
    colnames(counts) <- round(msim$time, 1)
    return(counts)
  } else {
    return(matrix(NA, nrow = nrow(params$M), ncol = n_t)) # Return properly shaped NA matrix
  }
}

#'-----------------------Separate table by chunks-------------#

library(parallel)

num_cores <- detectCores() - 3  # Use one less than the total number of cores
cat("The number of cores that will be used are: ", num_cores, "\n")

split_table <- function(df, n_chunks) {
  split(df, cut(seq_len(nrow(df)), breaks = n_chunks, labels = FALSE))
}

chunks <- split_table(params_table, num_cores)

#'-------------------------Generate workers directory-------------#

mias_dir <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/D25M02/Simulate_miaSim"
ode_dir <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/D25M02/Simulate_ODE"

# Function to create main and worker directories
create_dirs <- function(main_dir, num_cores) {
  if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)
  
  # Create worker directories
  worker_dirs <- file.path(main_dir, paste0("worker_", seq_len(num_cores)))
  invisible(lapply(worker_dirs, dir.create, showWarnings = FALSE))
  
  return(worker_dirs)
}

# Create directories for both simulations
workers_mia <- create_dirs(mias_dir, num_cores)
workers_ODE <- create_dirs(ode_dir, num_cores)

#'-------------------------Function for repeating simulations--------#

parsims <- function(index, workers_ODE, workers_miaSim) {
    
    # Generate parameters
    params <- regenerate(index)
    
    # Run simulation
    output_ode <- ode_function(times = 700, params)
    output_mia <- sim_glv(params = params, n_t = 700)
    
    # Define paths
    id <- params$ID
    save_ode <- file.path(workers_ODE, paste0("O_", id, ".tsv"))
    save_mia <- file.path(workers_miaSim, paste0("O_", id, ".tsv"))
    
    # Save simulation
    utils::write.table(output_ode, file = save_ode, sep = "\t", row.names = FALSE, col.names = TRUE)
    utils::write.table(output_mia, file = save_mia, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    return(id)
}


#'-----------------------------------------------------------------#

completed_ids <- mclapply(1:num_cores, function(core_id) {
  
  cat("Starting worker ", core_id, "....\n")
  
  core_chunk <- chunks[[core_id]]  # rows assigned to this core
  ids_vector <- lapply(1:nrow(core_chunk), function(i) parsims(core_chunk[i, ], workers_ODE, workers_mia))
  
  cat("Ending worker ", core_id, "....\n")
  
  return(as.vector(unlist(ids_vector)) )
  
}, mc.cores = num_cores)


cat("The number of simulations performed were: ", 
    length(as.vector(unlist(completed_ids))), "\n")

cat(
  paste0(rep("=", 20), collapse = ""), "  Ending code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)

