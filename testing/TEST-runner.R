#' ---
#' title: "Testing notebook"
#' author: "Manuel Rivera"
#' date: "2025-05-20"
#' output: 
#'   html_document:
#'   toc: true           # Enable Table of Contents
#' theme: sandstone
#' toc_float: true
#' collapsed: false
#' smooth_scroll: true
#' ---


#' Indicate directories paths
BASE_DIR <- "/mnt/atgc-d3/sur/users/mrivera/glv-research"
DATA_DIR <- file.path(BASE_DIR, "Data")
RESULTS_DIR <- file.path(BASE_DIR, "Results")
exp_id <- paste0("TEST01-D", format(Sys.Date(), "%d-%B"))

#========================================================
exp_dir <- file.path(RESULTS_DIR, exp_id) # Experiment directory
params_path <- file.path(DATA_DIR, paste0(exp_id, ".tsv")) # Parameters TSV
mc_dir <- file.path(RESULTS_DIR, exp_id, "mc-apply") # Workers directory
outs_path <- file.path(RESULTS_DIR, exp_id, "Outputs") # Outputs directory
info_path <- file.path(RESULTS_DIR, exp_id, "sims-summary.tsv") # Information TSV

#' First we generate the parameters for simulation:
#+ eval=FALSE

tictoc::tic("Section 0: Total running time")

# ==== Generate parameters ====
tictoc::tic("Section 1: Time for Parameter Generation")

#' Generate grid of parameters:
#+ eval=FALSE
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

generate_params <- function (){
 dt <- tidyr::expand_grid(
    n_species = seq(20, 100, by = 20),
    p_neg = seq(0, 1, by = 0.1), 
    p_noint = seq(0, 1, by = 0.1)
  ) %>%
  dplyr::filter(!(p_noint == 0 & p_neg == 0)) %>% # Exclude where p_noint=p_neg = 0
  tidyr::uncount(weights = 10) %>%
  dplyr::mutate(
    id = ids::random_id(n = dplyr::n(), bytes = 5),
    x0_seed = sample(1:1e6, dplyr::n(), replace = FALSE),
    mu_seed = sample(1:1e6, dplyr::n(), replace = FALSE),
    A_seed = sample(1:1e6, dplyr::n(), replace = FALSE)
  )
  
  return(dt)
}

#========================================================

# Generate parameters
params_df <- generate_params()

# Verify if ids are unique and in case they are, save the parameters.
while (nrow(params_df) != length(unique(params_df$id))) {
  params_df <- generate_params() # Repeat function
}

# Save them
data.table::fwrite(x = params_df, file = params_path, 
  sep = "\t",
  quote = FALSE,      # Disable quoting for faster writing
  row.names = FALSE,  # Don't write row names
)
message("\nParameteres generated and saved at path:\n", params_path, "\n")

tictoc::toc() # For section 1

#' We split the data into chunks of `n_cores`:
#+ eval=FALSE
# ==== Split data into chunks====
tictoc::tic("Section 2: Divide data into chunks")
library(parallel)
num_cores <- 60  # Use one less than the total number of cores

cat(rep("=", 40), "\n")
cat("The number of cores that will be used are: ", num_cores, "\n")
cat("Each core will do: ", nrow(params_df)/num_cores, " simulations\n")
cat(rep("=", 40), "\n")

split_table <- function(df, n_chunks) {
  split(df, cut(seq_len(nrow(df)), breaks = n_chunks, labels = FALSE))
}

chunks <- split_table(params_df, num_cores)
message("\nData split completed...\n")

tictoc::toc() # For section 2
#====================================================================
# Testing line
# num_cores <- 5
# chunks <- split_table(params_df[1:10], num_cores)
#====================================================================


#' We create a separate directory for each core to prevent race conditions (when two cores access the same directory simultaneously).
#+ eval=FALSE
# ==== Generate directories for each core ====
tictoc::tic("Section 3: Generate directories for each core")

# Generate workers directories
create_dirs <- function(main_dir, num_cores) {
  if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)
  
  # Create worker directories
  worker_dirs <- character(num_cores) # Preallocate vector
  worker_dirs <- file.path(main_dir, paste0("worker_", seq_len(num_cores)))
  sapply(worker_dirs, dir.create, showWarnings = FALSE, recursive = FALSE)
  invisible(lapply(worker_dirs, dir.create, showWarnings = FALSE))
  
  return(worker_dirs)
}
# Create directories 
workers_ODE <- create_dirs(mc_dir, num_cores)
message("\nWorking directories created at path:\n", mc_dir,"\n")
tictoc::toc() # For section 3

#' We source the function to generate the gLV parameters using the initial parameters and for solve it:


#===================================================
# Testing lines
# index = params_df[1,]
# times=500
# params = regenerate(index)
#===========================

#+ eval=FALSE
# Source functions to:
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/testing/generate-params.R") # regenerate parameters
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/testing/solver-gLV.R") # solve gLV equations

#' We wrap the code for parallelizing the simulations.
#+ eval=FALSE
# ==== Wrapper for running all required steps ====
sims_wrapper <- function(index, path_core) {
  params <- regenerate(index) # Generate parameters
  output <- solve_gLV(times = 500, params) # Run simulation

  # Define paths
  path_ode <- file.path(path_core, paste0("O_", params$id, ".feather")) # Output ODE
  NA_count <- as.integer(any(is.na(output))) # Calculate NAs
  arrow::write_feather(output, path_ode) # Save simulation

  cat(">> Simulation", params$id, "completed...\n", sep = " ")
  return(list(id = params$id, na_ct = NA_count))
}

#' We Parallelize the code and get the summary of simulations
#+ eval=FALSE
# ==== Parallelize it ====
tictoc::tic("Section 4: Run simulations using the parallel package")

core_id = 1
sims_info <- parallel::mclapply(1:num_cores, function(core_id) {
  
  message("Starting worker ", core_id, "....\n")
  
  core_chunk <- chunks[[core_id]]  # rows assigned to this core
  
  start_time <- Sys.time()  # Start timing
  result <- lapply(1:nrow(core_chunk), function(i) {
    sims_wrapper(index = core_chunk[i, ], path_core = workers_ODE[core_id])
  })

  end_time <- Sys.time()  # End timing
  elapsed <- difftime(end_time, start_time, units = "secs")
  message("Worker ", core_id, " finished in ", round(elapsed, 2), " seconds.\n")

  result <- data.table::rbindlist(result, use.names = TRUE) # Convert list to df
  message("Ending worker ", core_id, "....\n")
  
  return(result)
  
}, mc.cores = num_cores)

# Generate TSV file
sims_info_df <- data.table::rbindlist(sims_info) # Convert list (of df) to df
data.table::fwrite(x = sims_info_df, file = info_path, sep = "\t")
tictoc::toc() # For section 4

tictoc::toc() # For Total running time