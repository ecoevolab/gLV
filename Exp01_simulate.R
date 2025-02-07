cat(
  paste0(rep("=", 20), collapse = ""), "  Running code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)

#------------Load master table#------------#
params_table <- read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Params-exp01.tsv", sep = "\t", header = TRUE)

#'-----------------------function is for generating the parameters-------------#
regenerate <- function(index) {
  
  # Convert all required index values once
  N_species <- as.numeric(index["N_specs"])
  seeds <- as.numeric(index[c("Population_seed", "Growth_seed", "Interactions_seed")])
  probs <- as.numeric(index[c("Prob_neg", "Prob_0", "Diagonal")])
  
  # Generate Population parameter
  set.seed(seeds[1])
  Pobl <- runif(n = N_species, min = 0.1, max = 1)
  
  # Generate Growth rates parameter
  set.seed(seeds[2])
  Grow <- runif(n = N_species, min = 0.001, max = 1)
  
  # Generate Interactions parameter
  set.seed(seeds[3])
  P_neg <- rbinom(n = N_species^2, size = 1, prob = probs[1])
  interaction_values <- runif(n = N_species^2, min = 0, max = 1) * ifelse(P_neg == 1, -1, 1)
  V_noint <- rbinom(n = N_species^2, size = 1, prob = 1 - probs[2])
  inter <- matrix(V_noint * interaction_values, N_species)
  
  # Set diagonal values
  diag(inter) <- probs[3]
  
  # Return as a list
  list(Population = Pobl, Interactions = inter, Growths = Grow)
}


#------------------------------------Regenerate simulations----------------------------#

requireNamespace("deSolve")

ode_function <- function (times, params, atol, rtol) {
  
  # Define the equation
  glv_model <- function(t, x, params) {
    dx <- x * (params$Growths + params$Interactions %*% x)
    list(dx)
  }
  
  time_seq <- seq_len(times) # Faster sequence generation
  
  # Solve the ODE
  results <- deSolve::ode(y = params$Population, times = time_seq, func = glv_model, 
                          parms = params, method = "ode45", rtol = rtol, atol = atol)
  
  # Transform results: remove time column, transpose
  ode_df <- t(results[-1, -1, drop = FALSE])  # Efficient subsetting
  
  # Normalize each column (species proportions)
  col_sums <- colSums(ode_df)
  ode_df <- sweep(ode_df, 2, col_sums, "/")  # Faster than apply()
  
  return(ode_df)
}

#'-----------------------Separate table by chunks-------------#

library(parallel)

num_cores <- 20  # Use one less than the total number of cores
cat("The number of cores that will be used are: ", num_cores, "\n")

split_table <- function(df, n_chunks) {
  split(df, cut(seq_len(nrow(df)), breaks = n_chunks, labels = FALSE))
}

chunks <- split_table(params_table, num_cores)

#'-------------------------Generate workers directory-------------#

main_dir <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Experiment-01/Mc-run_02"

# Create the main directory if it doesn't exist
if (!dir.exists(main_dir)) {
  invisible(dir.create(main_dir, recursive = TRUE,  showWarnings = FALSE))
}

# Ensure each worker has its own directory
worker_dirs <- file.path(main_dir, paste0("worker_", seq_len(num_cores)))
invisible( lapply(worker_dirs, dir.create, showWarnings = FALSE) )
cat("Directories for the cores created, in total ", num_cores, " cores were created \n")

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

library("doParallel", lib = "/mnt/atgc-d3/sur/modules/pkgs/tidyverse_mrc")
cl <- makeCluster(num_cores)
registerDoParallel(cl)

completed_ids <- foreach(core_id = 1:num_cores, .combine = c, .packages = c("deSolve")) %dopar% {
  cat("Starting worker:", core_id, "\n")
  core_chunk <- chunks[[core_id]]
  worker_path <- worker_dirs[[core_id]]
  
  ids_vector <- sapply(seq_len(nrow(core_chunk)), function(i) {
    index <- as.list(core_chunk[i, ])
    params <- regenerate(index)
    
    reps_fun(grid, index, params, worker_path)
    return(index["ID_simulation"])
  })
  
  cat("Worker #", core_id, " completed.\n")
  return(ids_vector)
}

stopCluster(cl)

cat("The number of simulations repeated with the combination of tolerances was: ", 
    length(as.vector(unlist(completed_ids))), "\n\n")

#------------------------------------------------------------------#
cat(
  paste0(rep("=", 20), collapse = ""), "  Ending code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)



