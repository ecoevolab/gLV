cat(
  paste0(rep("=", 20), collapse = ""), "  Running code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)

#------------Load master table#------------#
params_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Params-exp01.tsv")

#'-----------------------function is for generating the parameters-------------#
regenerate <- function(index) {
  
  N_species <- as.numeric(index["N_specs"])
  
  #------------------Populations-----------------------------#
  set.seed(as.numeric(index["Population_seed"]))
  Pobl <- stats::runif(N_species, min = 0.1, max = 1)
  
  #------------------------Growth Rates---------------------#
  set.seed(as.numeric(index["Growth_seed"]))
  Grow <- stats::runif(N_species, min = 0.001, max = 1)
  
  #--------------------Interactions-------------------------#c
  set.seed(as.numeric(index["Interactions_seed"]))
  
  # Generate random interaction values
  P_neg <- stats::rbinom(N_species * N_species, 1, as.numeric(index["Prob_neg"]))
  tmp <- stats::rbinom(N_species * N_species, 1, 1 - as.numeric(index["Prob_0"])) * ifelse(P_neg != 0, stats::runif(N_species * N_species, min = 0, max = 1), -stats::runif(N_species * N_species, min = 0, max = 1))
  inter <- matrix(tmp, nrow = N_species, ncol = N_species)
  
  # Set diagonal values
  diag(inter) <- as.numeric(index["Diagonal"])
  
  # Return parameters as a list
  params <- list(Population = Pobl,
                 Interactions = inter,
                 Growths = Grow)
  
  return(params)
}

#------------------------------------Regenerate simulations----------------------------#

requireNamespace("deSolve")

ode_function <- function (times, params, atol, rtol) {
  
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
                          rtol = rtol, 
                          atol = atol)
  
  # Remove the column of times and get the transversal, where rows are species and column generations
  ode_df <- as.matrix(t(results))[-1, -ncol(results)]
  
  # Get the proportion of each specie by column value/columnsum
  for (c in 1:ncol(ode_df)) {
    csum <- sum(ode_df[,c])
    ode_df[ , c] =  ode_df[ , c] / csum
    
  }
  
  return(ode_df)
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

main_dir <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Experiment-01/Mc_run01"

# Create the main directory if it doesn't exist
if (!dir.exists(main_dir)) {
  dir.create(main_dir, recursive = TRUE)
}

# Ensure each worker has its own directory
worker_dirs <- file.path(main_dir, paste0("worker_", seq_len(num_cores)))
lapply(worker_dirs, dir.create, showWarnings = FALSE)

#'-------------------------Function for repeating simulations--------#

rtol_values <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1) # Relative tolerance
atol_values <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1) # Absolute tolerance
grid <- expand.grid(rtol = rtol_values, atol = atol_values)

reps_fun <- function(grid, index, params, worker_path) {
  
  # Loop in reverse order
  apply(grid, 1, function(tol) {
    
    # Get tolerances
    a <- tol["atol"]
    r <- tol["rtol"]
    
    # Run simulation
    new_out <- ode_function(times = 700, params, a, r)
    
    # Define paths
    tmp_path <- file.path(worker_path, paste0("tol_r", format(r, scientific = TRUE, digits = 3) ), 
                          paste0("tol_a", format(a, scientific = TRUE, digits = 3) ) )
    
    # Create the main directory if it doesn't exist
    if (!dir.exists(tmp_path)) {
      dir.create(tmp_path, recursive = TRUE)
    }
    
    # Save simulation
    save_path <- file.path(tmp_path, paste0("id_", index["ID_simulation"], ".tsv")  )
    utils::write.table(new_out, file = save_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  })
}


#'-----------------------------------------------------------------#

completed_ids <- mclapply(1:num_cores, function(core_id) {
  
  core_chunk <- chunks[[core_id]]  # rows assigned to this core
  worker_path <- worker_dirs[[core_id]]  # Worker directory
  
  ids_vector <- lapply(seq_len(nrow(core_chunk)), function(i) {
    index <- as.list(core_chunk[i, ])
    params <- regenerate(index)
    
    reps_fun(grid, index, params, worker_path)
    return(index["ID_simulation"])
  })
  
  return(as.vector(unlist(ids_vector)) )
  
}, mc.cores = num_cores)

cat("The number of simulations repeated with the combination of tolerances was: ", 
    length(as.vector(unlist(completed_ids))), "\n")

cat(
  paste0(rep("=", 20), collapse = ""), "  Ending code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)



