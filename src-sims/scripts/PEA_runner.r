# This script is done to generate the results for the simulation study.

# ===================================
# Generate directories
# ===================================
tictoc::tic("Section 0: Total running time")

tictoc::tic("Section 1: Time for Parameter Generation")
#' Indicate directories paths
pdir <- "/mnt/data/sur/users/mrivera/Controls"    # parent directory
exp_id <- 'PEA'                                   # experiment name                                     
exp_dir <- file.path(pdir, exp_id)                # experiment directory

params_path <- file.path(exp_dir, "sim-params.tsv")     # Parameters-TSV
mc_dir <- file.path(exp_dir, "mc-apply")                # Workers-dir
info_path <- file.path(exp_dir, "pea-summary.feather")  # Information-TSV
cat(">> The experiment path is:", exp_dir,"\n", sep="")

# ===================================
# Generate parameters for simulation
# ===================================
generate_params <- function (){
  p_neg = seq(0, 1, by = 0.1)       # negative-interactions
  p_noint = seq(0, .9, by = 0.1)     # null-interactions                    
  n_species = rep(seq(20, 100, by = 20), 10) # number of species               
  dt <- data.table::CJ(n_species, p_neg, p_noint )

  # Add Columns with data table operator `:=`
  n_total <- nrow(dt)
  all_seeds <- sample.int(3e6L, 3L *n_total, replace = FALSE)
  
  dt[, `:=`(
    id = ids::random_id(n = n_total, bytes = 3),
    x0_seed = all_seeds[1:n_total],
    mu_seed = all_seeds[(n_total+1):(n_total*2)],
    A_seed = all_seeds[(2 * n_total + 1):(3 * n_total)]
  )]
  
  return(dt)
}

params_df <- generate_params()
# Verify if ids are unique and in case they are, save the parameters.
while (nrow(params_df) != length(unique(params_df$id))) {
    params_df <- generate_params() # Repeat function
}

data.table::fwrite(x = params_df, file = params_path, sep = "\t", quote = FALSE, row.names = FALSE) # Save parameters
message("\nParameteres generated and saved at path:\n", params_path, "\n")
cat(">> The number of simulations are:", nrow(params_df),"\n", sep=" ")
tictoc::toc() # For section 1

# ===================================
# Divide data into chunks
# ===================================
tictoc::tic("Section 2: Divide data into chunks")

library(parallel)
num_cores <- parallel::detectCores() - 1  # Use one less than the total number of cores
cat("The number of cores that will be used are: ", num_cores, "\n")

split_table <- function(df, n_chunks) {
  split(df, cut(seq_len(nrow(df)), breaks = n_chunks, labels = FALSE))
}

# Divide data into chunks
chunks <- split_table(params_df, num_cores)
message("\nData split completed...\n")
tictoc::toc() # For section 2

# ===================================
# Create-worker-directories
# ===================================
tictoc::tic("Section 3: Generate directories for each core")

# Generate workers directories
create_dirs <- function(mc_dir, num_cores) {
  if (!dir.exists(mc_dir)) dir.create(mc_dir, recursive = TRUE)
  
  # Create worker directories
  worker_dirs <- character(num_cores) # Preallocate vector
  worker_dirs <- file.path(mc_dir, paste0("worker_", seq_len(num_cores)))
  sapply(worker_dirs, dir.create, showWarnings = FALSE, recursive = FALSE)
  invisible(lapply(worker_dirs, dir.create, showWarnings = FALSE))
  
  return(worker_dirs)
}
# Create directories 
workers_ODE <- create_dirs(mc_dir, num_cores)
message("\nWorking directories created at path:\n", mc_dir,"\n")
tictoc::toc() # For section 3

# ===================================
# Source functions and simulate
# ===================================
codes = list.files('/mnt/data/sur/users/mrivera/gLV/src-sims/FUN', full.names=TRUE)

lapply(codes, function(file){
  cat(">> Sourcing function: ", file, "\n")
  capture.output(source(file))
  return()
})

# Wrap functions
wrapper <- function(index, path_core) {
  #=================== Output ===================
  params <- regenerate(index)                                               # Generate-parameters
  output <- solve_gLV(times = 1000, params)                                 # Run-simulation
  out_path <- file.path(path_core, paste0("O_", params$id, ".feather"))     # Simulation-path
  A_path <- file.path(path_core, paste0("A_", params$id, ".feather"))       # Mat-path
  arrow::write_feather(x = output, sink = out_path)                         # Saving-ODE
  arrow::write_feather(as.data.frame(params$M), A_path)                     # Save-MAT
  #=================== Information ===================
  NA_count <- sum(is.na(output))                            # simulation-NAs
  tts_out <- find_ts(output)                                # time-to-stability OUTPUT
  cat(">> Simulation ", params$id, " completed.\n")
  return(list(id = params$id, na_ct = NA_count, tts_out = tts_out))
}

# Parallelize code
tictoc::tic("Section 4: Run simulations and extinctions using the parallel package")

sims_info <- parallel::mclapply(1:num_cores, function(core_id) {
  message("Starting worker ", core_id, "....\n")
  core_chunk <- chunks[[core_id]]  # rows assigned to this core
  result <- lapply(1:nrow(core_chunk), function(i) {
    wrapper(index = core_chunk[i, ], path_core = workers_ODE[core_id])
  })
  result <- data.table::rbindlist(result, use.names = TRUE) # Convert list to df
  message("Ending worker ", core_id, "....\n")
  return(result)
}, mc.cores = num_cores)

# ===================================
# Save simulation information
# ===================================
# Generate TSV file
sims_info_df <- data.table::rbindlist(sims_info) # Convert list (of df) to df
arrow::write_feather(sims_info_df, info_path)
tictoc::toc() # For section 4

# Create symbolic links
tictoc::tic("Section 5: Generate symbolic links")
gen_syml(mc_dir)
files=list.files(mc_dir, recursive=TRUE)
unlink(mc_dir, recursive = TRUE)
tictoc::toc() # For section 5
tictoc::toc() # For Total running time
