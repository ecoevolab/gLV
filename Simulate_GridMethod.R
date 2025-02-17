cat(
  paste0(rep("=", 20), collapse = ""), "  Running code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)

#------------Load master table------------#
params_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/D13M02Y25.tsv")

#'-----------------------function is for generating the parameters-------------#
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
  params <- list(Population = Pobl,
                 Interactions = inter,
                 Growths = Grow,
                 ID = id)

  return(params)
}

#'-----------------------function is for solving the gLV equation-------------#
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
  results <- deSolve::ode(y = params$Population, times = time_seq, func = glv_model, parms = params,
                          method = "ode45",
                          rtol = 1e-06, 
                          atol = 1e-06)
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

main_dir <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/D13M02Y25"

# Create the main directory if it doesn't exist
if (!dir.exists(main_dir)) {
  dir.create(main_dir, recursive = TRUE)
}

# Ensure each worker has its own directory
worker_dirs <- file.path(main_dir, paste0("worker_", seq_len(num_cores)))
lapply(worker_dirs, dir.create, showWarnings = FALSE)

#'-------------------------Function for repeating simulations--------#

parsims <- function(index, worker_path) {
    
    # Generate parameters
    params <- regenerate(index)
    
    # Run simulation
    output <- ode_function(times = 700, params)
    
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

