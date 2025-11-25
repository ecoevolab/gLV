#' ---
#' title: "ExpM279 notebook"
#' author: "Manuel Rivera"
#' date: "2025-10-30"
#' output: 
#'   html_document:
#'   toc: true           # Enable Table of Contents
#' theme: sandstone
#' toc_float: true
#' collapsed: false
#' smooth_scroll: true
#' ---

#============================================================================
# SECTION: Generate-ID-and-paths
tictoc::tic("Section 0: Total running time")

#' Indicate directories paths
pdir <- '/mnt/data/sur/users/mrivera/Controls'                              # Parent-dir
# exp_id <- substr(ids::uuid(1, use_time = TRUE), 1, 13)                    # exp-id
exp_id <- paste0('exp_', format(Sys.Date(), "%Y%m%d")) 
exp_dir = file.path(pdir, exp_id)                                      

# Create directory experiment directory
dir.create(exp_dir)   
mc_dir <- file.path(exp_dir, "mc-apply")                                        # Workers-dir
info_path <- file.path(exp_dir, "Summary-simulations.feather")                  # Information-TSV
cat(">> The experiment path is:", exp_dir,"\n", sep="")

#============================================================================
#' First we generate the parameters for simulation:
#+ eval=FALSE
# SECTION: generate-parameters
tictoc::tic("Section 1: Time for Parameter Generation")

#' Generate grid of parameters:
#+ eval=FALSE
library(data.table)

generate_params <- function (){
  reps = 1000
  dt <-data.frame(p_noint = rep(0.7, reps), n_species = rep(30, reps), keys = sample(1:30, reps, replace = TRUE))
  dt <- as.data.table(dt)  # Convert to data.table
  all_seeds <- sample.int(3e6L, 3L * reps, replace = FALSE)
  
  dt[, `:=`(
    id = ids::random_id(n = reps, bytes = 5),
    x0_seed = all_seeds[1:reps],
    mu_seed = all_seeds[(reps+1):(reps*2)],
    A_seed = all_seeds[(2 * reps + 1):(3 * reps)]
  )]
  
  return(dt)
}
# Generate unique parameters 
params_df <- generate_params()
n_unique <- uniqueN(params_df$id)  

while (nrow(params_df) != n_unique) {
  params_df <- generate_params()
  n_unique <- uniqueN(params_df$id)
}

# Save parameters
params_path <- file.path(exp_dir, "Simulation-parameters.tsv")  
data.table::fwrite(params_df, params_path, sep = "\t")  
message("\nParameters generated and saved at path:\n", params_path)

# Print summary
n_sims <- nrow(params_df)
cat(">> The number of simulations is:", n_sims, "\n")
cat(">> The number of extinctions to do is:", n_sims * 30, "\n")  
tictoc::toc() # For section 1

#============================================================================
# SECTION: Divide-data-into-chunks

#' We split the data into chunks of `n_cores`:
#+ eval=FALSE
tictoc::tic("Section 2: Divide data into chunks")

library(parallel)
num_cores <- parallel::detectCores() - 1  # Use one less than the total number of cores
if (num_cores > nrow(params_df)) {
  num_cores <- nrow(params_df)
}
cat("The number of cores that will be used are: ", num_cores, "\n")

split_table <- function(df, n_chunks) {
  split(df, cut(seq_len(nrow(df)), breaks = n_chunks, labels = FALSE))
}

chunks <- split_table(params_df, num_cores)
message("\nData split completed...\n")
tictoc::toc() # For section 2

#============================================================================
# SECTION: Create-worker-directories
#' We create a separate directory for each core to prevent race conditions (when two cores access the same directory simultaneously).
#+ eval=FALSE
# ==== Generate directories for each core ====
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

#============================================================================
# SECTION: Source-codes
#' We source the function to generate the gLV parameters using the initial parameters and for solve it:

#+ eval=FALSE
# Source functions to:
codes = list.files('/mnt/data/sur/users/mrivera/gLV/src-sims/CTRL_FUN', full.names=TRUE)

lapply(codes, function(file){
  cat(">> Sourcing function: ", file, "\n")
  capture.output(source(file))
  return()
})

#============================================================================
# SECTION: Wrap-functions
#' We wrap the code for parallelizing the simulations.
#+ eval=FALSE
wrapper <- function(index, path_core) {
  #============== Generate parameters ===============
  params <- ctrl_regenerate(index)
  id_str <- params$id
  #============== Build paths ===============
  out_path <- file.path(path_core, paste0("O_", id_str, ".feather"))
  A_path <- file.path(path_core, paste0("A_", id_str, ".feather"))
  preds_path <- file.path(path_core, paste0("tgt_", id_str, ".feather"))
  #============== Run simulation ===============
  output <- solve_gLV(times = 1000, params)
  #============== Write files ===============
  arrow::write_feather(output, out_path)
  arrow::write_feather(as.data.frame(params$M), A_path)
  #============== Information ===============
  NA_count <- sum(is.na(output))  # Keep for counting
  tts_out <- find_stability(output)
  cat(">> Simulation", id_str, "completed.\n")
  #============== Extinctions ==============
  params$x0 <- output[, 1000]                         # Stable-population 
  preds_df <- generate_extinctions(params, path_core)          # Generate-extinctions
  arrow::write_feather(preds_df, preds_path)          # Save predictions
  list(
    id = id_str,
    na_ct = NA_count,
    tts_out = tts_out,
    tts_ext = max(preds_df$ext_ts)
  )
}

#============================================================================
# SECTION: Parallelize-code
#' We Parallelize the code and get the summary of simulations
#+ eval=FALSE
# ==== Parallelize it ====
tictoc::tic("Section 4: Run simulations and extinctions using the parallel package")

sims_info <- parallel::mclapply(seq_len(num_cores), function(core_id) {

  message("Starting worker ", core_id, "....\n")
  # Get chunk for this core
  core_chunk <- chunks[[core_id]]
  n_rows <- nrow(core_chunk)

  # Pre-allocate list for results
  result <- vector("list", n_rows)
  path_core <- workers_ODE[core_id]  

  for (i in seq_len(n_rows)) {
    result[[i]] <- wrapper(index = core_chunk[i, ], path_core = path_core)
  }
  message("Ending worker ", core_id, "....\n")
  data.table::rbindlist(result, use.names = TRUE) # Return list as df
}, mc.cores = num_cores)

# Generate TSV file
sims_info_df <- data.table::rbindlist(sims_info) # Convert list (of df) to df
arrow::write_feather(sims_info_df, info_path)
# arrow::read_feather( info_path)
tictoc::toc() # For section 4

#============================================================================
# SECTION: Create-symbolic-links
#' We create symbolic links of the simulation...
#+ eval=FALSE
tictoc::tic("Section 5: Generate symbolic links")
gen_syml(mc_dir)
files=list.files(mc_dir, recursive=TRUE)
unlink(mc_dir, recursive = TRUE)
tictoc::toc() # For section 5

tictoc::toc() # For Total running time
