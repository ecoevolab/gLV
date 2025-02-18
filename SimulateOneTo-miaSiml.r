cat(
  paste0(rep("=", 20), collapse = ""), "  Running code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)

#`-----------------------------Load master table------------#'
params_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/D13M02Y25.tsv")

#'-----------------------function is for generating the parameters-------------#'
regenerate <- function(index) {
  
  N_species <- as.numeric(index[["n_species"]])
  
  #------------------Populations-----------------------------#
  set.seed(as.numeric(index[["Pop_seed"]]))
  Pobl <- stats::runif(N_species, min = 0.1, max = 1)
  
  #------------------------Growth Rates---------------------#
  set.seed(as.numeric(index[["Growth_seed"]]))
  Grow <- stats::runif(N_species, min = 0.001, max = 1)
  
  #--------------------Interactions-------------------------#c
  set.seed(as.numeric(index[["A_seed"]]))
  
  # Probability of negative interaction and vector of 0's and 1's
  p_neg <- as.numeric(index[["p_neg"]])
  V_neg <- stats::rbinom(N_species * N_species, 1, p_neg)
  
  # Probability of null interaction and vector of 0's and 1's
  p_noint <- as.numeric(index[["p_noint"]])
  V_noint <- stats::rbinom(N_species * N_species, 1, 1 - p_noint)
  
  tmp <- V_noint * ifelse(V_neg != 0,
                          stats::runif(N_species * N_species, min = 0, max = 1),
                          -stats::runif(N_species * N_species, min = 0, max = 1)
  )
  
  inter <- matrix(tmp, nrow = N_species, ncol = N_species)
  
  # Set diagonal values
  diag(inter) <- -0.5
  
  # Extract ID
  id <- index[["id"]]
  
  # Return parameters as a list
  params <- list(x0 = Pobl,
                 M = inter,
                 mu = Grow,
                 ID = id)
  
  return(params)
}

#' Simulate gLV
#' 
#' Wrapper for miaSim. Takes output from enerate_params and runs
#' gLV simulation. No stochasticity, measurement error, perturbation or
#' immigration is considered.
#' 
#' Checks for numerical issues after the simulation
#'
#' @param params A list with starting conditions (x0), growth rates (mu), and
#' interaction matrix (M) for gLV simulation
#' @param n_t Number of timepoints to simulate
#' @param timeout Number of seconds to wait for simulation to end before
#' stopping it
#'
#' @return A tibble where column sim has the results of the simulation
#' 
#' @export
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

main_dir <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/D13M02Y25-miaSim"

# Create the main directory if it doesn't exist
if (!dir.exists(main_dir)) {
  dir.create(main_dir, recursive = TRUE)
}

# Ensure each worker has its own directory
worker_dirs <- file.path(main_dir, paste0("worker_", seq_len(num_cores)))
invisible( lapply(worker_dirs, dir.create, showWarnings = FALSE) )

#'-------------------------Function for repeating simulations--------#

parsims <- function(index, worker_path) {
  
  # Generate parameters
  params <- regenerate(index)
  
  # Run simulation
  output <- sim_glv(params = params, n_t = 700)
  
  # Define paths
  id <- params$ID
  save_path <- file.path(worker_path, paste0("O_", id, ".tsv"))
  
  # Create the main directory if it doesn't exist
  dir.create(dirname(save_path), recursive = TRUE,  showWarnings = FALSE)
  
  # Save simulation
  utils::write.table(output, file = save_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  return(id)
}


#'-----------------------------------------------------------------#

completed_ids <- mclapply(1:num_cores, function(core_id) {
  
  cat("Starting worker ", core_id, "....\n")
  
  core_chunk <- chunks[[core_id]]  # rows assigned to this core
  worker_path <- worker_dirs[[core_id]]  # Worker directory
  ids_vector <- lapply(1:nrow(core_chunk), function(i) parsims(core_chunk[i, ], worker_path))
  
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

