#' ---
#' title: "ExpM279 notebook"
#' author: "Manuel Rivera"
#' date: "2025-03-04"
#' output: 
#'   html_document:
#'   toc: true           # Enable Table of Contents
#' theme: sandstone
#' toc_float: true
#' collapsed: false
#' smooth_scroll: true
#' ---


#' Indicate directories paths
pdir <- "/home/mrivera"                                                                 # Parent-dir
dat_dir <- file.path(pdir, "Data")                                                      # data-dir
res_dir <- file.path(pdir, "Experiments")                                               # experiment-dir
exp_id <- paste0(substr(ids::proquint(), 1, 5),"-D", format(Sys.Date(), "%d%b"))        # exp-id
exp_dir <- file.path(res_dir, exp_id)                                                   # Experiment directory
params_path <- file.path(dat_dir, paste0(exp_id, ".tsv"))               # Parameters TSV
mc_dir <- file.path(res_dir, exp_id, "mc-apply")                        # Workers directory
info_path <- file.path(res_dir, exp_id, "sims-summary.tsv")             # Information TSV

#' First we generate the parameters for simulation:
#+ eval=FALSE

tictoc::tic("Section 0: Total running time")

# ==== Generate parameters ====
tictoc::tic("Section 1: Time for Parameter Generation")

#' Generate grid of parameters:
#+ eval=FALSE
generate_params <- function (){
  dt <- data.table::CJ(n_species = rep(c(20, 100), times = 300), p_neg = 1, p_noint = seq(0, 1, by = 0.05))

  # Add Columns with data table operator `:=`
  dt[, `:=`(
    id = ids::random_id(n = .N, bytes = 5),
    x0_seed = sample(1:1e6, .N, replace = FALSE),
    mu_seed = sample(1:1e6, .N, replace = FALSE),
    A_seed = sample(1:1e6, .N, replace = FALSE)
  )]
  
  return(dt)
}
#====================================================================
# Testing line
params_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Exp06-D29-Apr.tsv"
#========================================================
# Save Parameters as TSV or load them
if (!(file.exists(params_path))) {

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
} else {
  params_df <- data.table::fread(params_path)
  message("\nParameters loaded..\n", params_path)
}
tictoc::toc() # For section 1

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

#====================================================================
# Testing line
# num_cores <- 5
# chunks <- split_table(params_df[1:10], num_cores)
#====================================================================
chunks <- split_table(params_df, num_cores)
message("\nData split completed...\n")
tictoc::toc() # For section 2

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

#+ eval=FALSE
# Source functions to:
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/src/generate-params.R") # regenerate parameters
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/src/solver-gLV.R") # solve gLV equations
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/src/forge-symls.R") # generate symbolic links
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/src/find-time-stability.R") # source function for ts
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/src/extinctions-fun.R") # Extinction function

#' We wrap the code for parallelizing the simulations.
#+ eval=FALSE
# ==== Wrapper for running all required steps ====
sims_wrapper <- function(index, path_core) {
  params <- regenerate(index)                       # Generate parameters
  output <- solve_gLV(times = 1000, params)         # Run simulation

  # Define paths
  path_ode <- file.path(path_core, paste0("O_", params$id, ".tsv"))       # Output
  ext_path <-  file.path(path_core, paste0("E_", params$id, "-Info.tsv")) # Extinctions

  NA_count <- sum(is.na(output))            # Calculate NAs
  ts_out <- find_ts(output)                 # time-to-stability (ts) OUTPUT

  # Save simulation ODE
  # data.table::fwrite(x = output, file = path_ode, sep = "\t", row.names = FALSE, col.names = FALSE) 
  arrow::write_feather(x = output, sink = path_ode)
  
  #=================== Extinctions ===================
  ext_list <- sim_all_ext(output, params)                                     # Generate extinctions
  ts_ext <- max(ext_list$info[["ext_ts"]])                                    # Find tts on EXTINCTIONS
  tmp <- names(ext_list$data)                                                 # #-Spec (S1, S2, etc)
  ext_paths <- paste0(path_core, "/E_", params$id, "-", tmp, ".feather")      # Extinctions paths

  # Save extinctions-OUTPUT
  lapply(seq_along(ext_list$data), function(e) {
    arrow::write_feather(ext_list$data[[e]], ext_paths[e])
  })
  # Save extinctions-INFO
  arrow::write_feather(ext_list$info, ext_path) 

  cat("Simulation ", params$id, " completed...\n")
  return(list(id = params$id, na_ct = NA_count, tts_out = ts_out, tts_ext = ts_ext))
}

#' We Parallelize the code and get the summary of simulations
#+ eval=FALSE
# ==== Parallelize it ====
tictoc::tic("Section 4: Run simulations and extinctions using the parallel package")

sims_info <- parallel::mclapply(1:num_cores, function(core_id) {
  
  message("Starting worker ", core_id, "....\n")
  
  core_chunk <- chunks[[core_id]]  # rows assigned to this core
  
  result <- lapply(1:nrow(core_chunk), function(i) {
    sims_wrapper(index = core_chunk[i, ], path_core = workers_ODE[core_id])
    
  })

  result <- data.table::rbindlist(result, use.names = TRUE) # Convert list to df
  message("Ending worker ", core_id, "....\n")
  
  return(result)
  
}, mc.cores = num_cores)

# Generate TSV file
sims_info_df <- data.table::rbindlist(sims_info) # Convert list (of df) to df
data.table::fwrite(x = sims_info_df, file = info_path, sep = "\t")
tictoc::toc() # For section 4

#' We create symbolic links of the simulation...
#+ eval=FALSE
# ==== Create symbolic links====
# tictoc::tic("Section 5: Generate symbolic links")
# purrr::walk(workers_ODE, ~gen_syml(.x, exp_dir))
# tictoc::toc() # For section 5
tictoc::toc() # For Total running time