#' ---
#' title: "Experiment05 code raw"
#' output:
#'   pdf_document:
#'     keep_tex: true
#' date: 2025-03-20
#' ---

#' First we generate the parameters for simulation:
#+ eval=FALSE

# Load libraries
library(tictoc)
library(tidyverse)
library(tidyr)

tictoc::tic("Section 0: Total running time")
tictoc::tic("Section 1: Time for Parameter Generation:")

# ==== Generate parameters ====
# Generate experiment name
generate_id <- function() {
  day <- format(Sys.Date(), "%d")
  month <- format(Sys.Date(), "%m")
  
  # char_string <- paste0(sample(c(LETTERS, 0:9), 4, replace = TRUE), collapse = "")
  paste0("Exp05", "-D", day, "M", month)
}

wd <- "/mnt/atgc-d3/sur/users/mrivera/glv-research"
exp_id <- generate_id()
params_path <- file.path(wd, "Data", paste0(exp_id, ".tsv"))

#' Generate grid of parameters:
#+ eval=FALSE
# params_to_sim <- expand_grid(n_species = rep(c(20, 100), times = 300), p_neg = 1, p_noint = seq(0, 1, by = 0.05)) %>%
#   mutate(id = ids::random_id(n = length(n_species), bytes = 5)) %>%
#   mutate(x0_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE))) %>%
#   mutate(mu_seed = as.vector(sample(1:1e6, length(n_species),replace = FALSE))) %>%
#   mutate(A_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE)))
# 
# # Verify if ids are unique and in case they are, save the parameters.
# if (nrow(params_to_sim) == length(unique(params_to_sim$id))) {
#   # Save Parameters as TSV
#   data.table::fwrite(x = params_to_sim, file = params_path, sep = "\t")
#   message("\nParameteres generated and saved with ID: ",  paste0(exp_id, ".tsv"), "\n")
# } else {
#   message("\nAn error ocurred. IDs are not unique\n")
# }

# Read TSV file
params_to_sim <- data.table::fread(params_path, sep = "\t")
tictoc::toc()


#' We split the data into chunks of `n_cores`:
#+ eval=FALSE
# ==== Split data into chunks====
tictoc::tic("Section 2: Divide data into chunks")
library(parallel)
num_cores <- parallel::detectCores() - 1  # Use one less than the total number of cores
cat("The number of cores that will be used are: ", num_cores, "\n")

split_table <- function(df, n_chunks) {
  split(df, cut(seq_len(nrow(df)), breaks = n_chunks, labels = FALSE))
}

chunks <- split_table(params_to_sim, num_cores)
message("\nData split completed...\n")
tictoc::toc() # For section 2

#' We create a separate directory for each core to prevent race conditions (when two cores access the same directory simultaneously).
#+ eval=FALSE
# ==== Generate directories for each core ====
tictoc::tic("Section 3: Generate directories for each core")

# Generate workers directories
ode_dir <- file.path(wd, "Results", exp_id, "Parallel-02")

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
message("\nWorking directories created at path:\n", ode_dir,"\n")
tictoc::toc() # For section 3

#' We source the function to generate the gLV parameters using the initial parameters and for solve it:

#+ eval=FALSE
# Source function to regenerate parameters
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/D20M03-tmpdir/gLV-Params.R")
# print(regenerate)

#+ eval=FALSE
# Source function to solve gLV equations
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/D20M03-tmpdir/solve-gLV.R")
# print(solve_gLV)

#+ eval=FALSE
# Source function to generate extinctions
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/D20M03-tmpdir/extinctions-fun.R")
# print(simulate_all_extinctions)

#' We wrap the code for parallelizing the simulations.
#+ eval=FALSE
# ==== Wrapper for running all required steps ====
parallel.sims <- function(index, path_core) {
  
  # Generate parameters
  params <- regenerate(index)
  
  # Run simulation
  output_ode <- solve_gLV(times = 1000, params)
  
  # Run extinctions
  ext_tb <- simulate_all_extinctions(output_ode, params)
  
  # Define paths
  path_ode <- file.path(path_core, paste0("O_", params$id, ".tsv"))
  path_ext <- file.path(path_core, paste0("Ext_", params$id, ".rds"))
  
  # Calculate NAs
  NA_count <- sum(is.na(output_ode))
  
  # Save simulation
  utils::write.table(output_ode, file = path_ode, sep = "\t", row.names = FALSE, col.names = TRUE)
  saveRDS(ext_tb, file = path_ext)
  
  return(c(id = params$id, NA_count = NA_count))
}

#' We Parallelize the code and get the `NA counting`. 
#+ eval=FALSE
# ==== Parallelize it ====
tictoc::tic("Section 4: Run simulations and extinctions using the parallel package")

nas_counts <- parallel::mclapply(1:num_cores, function(core_id) {
  
  message("Starting worker ", core_id, "....\n")
  
  core_chunk <- chunks[[core_id]]  # rows assigned to this core
  
  na.vec <- lapply(1:nrow(core_chunk), function(i) {
    parallel.sims(core_chunk[i, ], path_core = workers_ODE[core_id])
  })
  
  message("Ending worker ", core_id, "....\n")
  
  return(na.vec)
  
}, mc.cores = num_cores)
tictoc::toc() # For section 4


# ==== Get NAs number on simulations====
tictoc::tic("Section 5: Count total number of NAs")
counts_df <- as.data.frame(matrix(unlist(nas_counts), ncol = 2, byrow = TRUE))
colnames(counts_df) <- unique(names(nas_counts))

# Save Parameters as TSV
na_path <- file.path(wd, "Results", exp_id, "Na-count-02.tsv")
data.table::fwrite(x = counts_df, file = na_path, sep = "\t")
tictoc::toc() # For section 5


#' We create symbolic links of the simulation...
#+ eval=FALSE
# ==== Create symbolic links====
tictoc::tic("Section 6: Generate symbolic links")

source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/D20M03-tmpdir/Forge_symlinks.R")

# Define source and target directories
workers_path <- file.path(wd, "Results", exp_id, "Parallel-02")
outs_path <- file.path(wd, "Results", exp_id, "Outputs-02")
exts_path <- file.path(wd, "Results", exp_id, "Extinctions-02")
generate_symlinks(workers_path, outs_path, exts_path)
tictoc::toc() # For section 6
tictoc::toc() # For Total running time