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
BASE_DIR <- "/mnt/atgc-d3/sur/users/mrivera/glv-research"
DATA_DIR <- file.path(BASE_DIR, "Data")
RESULTS_DIR <- file.path(BASE_DIR, "Results")
exp_id <- paste0("Exp06-D", format(Sys.Date(), "%d-%b"))

params_path <- file.path(DATA_DIR, paste0(exp_id, ".tsv")) # Parameters TSV
mc_dir <- file.path(RESULTS_DIR, exp_id, "mc-apply") # Workers directory
outs_path <- file.path(RESULTS_DIR, exp_id, "Outputs") # Outputs directory
info_path <- file.path(RESULTS_DIR, exp_id, "Info.tsv") # Information TSV

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

params_df <- generate_params()

# Verify if ids are unique and in case they are, save the parameters.
while (nrow(params_df) != length(unique(params_df$id))) {
  params_df <- generate_params() # Repeat function
}

# Save Parameters as TSV
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
num_cores <- parallel::detectCores() - 1  # Use one less than the total number of cores
cat("The number of cores that will be used are: ", num_cores, "\n")

split_table <- function(df, n_chunks) {
  split(df, cut(seq_len(nrow(df)), breaks = n_chunks, labels = FALSE))
}

# Testing line:
num_cores <- 5
chunks <- split_table(params_df[1:10], num_cores)

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

#' We wrap the code for parallelizing the simulations.
#+ eval=FALSE
# ==== Wrapper for running all required steps ====
par_sims <- function(index, path_core) {
  params <- regenerate(index) # Generate parameters
  output <- solve_gLV(times = 1000, params) # Run simulation

  # Define paths
  path_ode <- file.path(path_core, paste0("O_", params$id, ".tsv")) # Output
  ext_path <-  file.path(path_core, paste0("E_", params$id, "-Info.tsv")) # Extinctions

  NA_count <- sum(is.na(output)) # Calculate NAs
  ts_out <- find_ts(output) # Find time-to-stability (ts) OUTPUT
  data.table::fwrite(x = output, file = path_ode, sep = "\t", row.names = FALSE, col.names = FALSE) # Save simulation
  
  # Simulate extinctions
  ext_list <- sim_all_ext(output, params) # Generate extinctions
  ts_ext <- max(ext_list$info[["ext_ts"]]) # Find ts on EXTINCTIONS

  tmp <- names(ext_list$data)
  ext_paths <- paste0(path_core, "/E_", params$id, "-", tmp, ".feather") # Extinctions paths

  # Save EXTINCTIONS output
  lapply(seq_along(ext_list$data), function(e) {
    arrow::write_feather(ext_list$data[[e]], ext_paths[e])
  })
  arrow::write_feather(ext_list$info, ext_path) # Save EXTINCTIONS info
  
  return(list(id = params$id, na_ct = NA_count, ts_out = ts_out, ts_ext = ts_ext ))
}

#' We Parallelize the code and get the `NA counting`. 
#+ eval=FALSE
# ==== Parallelize it ====
tictoc::tic("Section 4: Run simulations and extinctions using the parallel package")

sims_info <- parallel::mclapply(1:num_cores, function(core_id) {
  
  message("Starting worker ", core_id, "....\n")
  
  core_chunk <- chunks[[core_id]]  # rows assigned to this core
  
  result <- lapply(1:nrow(core_chunk), function(i) {
    par_sims(index = core_chunk[i, ], path_core = workers_ODE[core_id])
  })

  result <- rbindlist(result, use.names = TRUE) # Convert list to df
  message("Ending worker ", core_id, "....\n")
  
  return(result)
  
}, mc.cores = num_cores)


sims_info_df <- data.table::rbindlist(sims_info) # Convert list (of df) to df
data.table::fwrite(x = sims_info_df, file = info_path, sep = "\t")
tictoc::toc() # For section 4



#' We create symbolic links of the simulation...
#+ eval=FALSE
# ==== Create symbolic links====
tictoc::tic("Section 6: Generate symbolic links")

# Define source and target directories
generate_symlinks(mc_dir, outs_path)
tictoc::toc() # For section 6



#' We look for stability time...
#+ eval=FALSE
# ==== Create symbolic links====
tictoc::tic("Section 7: Finding stable time")

source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/src/find-time-stability.R") # source function for ts

# List all outputs files
files <- list.files(outs_path, full.names = TRUE)
chunks <- split(files, cut(x = seq_along(files), breaks = num_cores, labels = FALSE))

tsa_list <- parallel::mclapply(seq_len(num_cores), function(core) {

  message("Starting worker ", core, "....\n")

  files <- chunks[[core]]  # rows assigned to this core

  res <- lapply(files, FUN = function(file) {
    id <- sub(".*/O_(\\w+)\\.tsv$", "\\1", file) # Sub id
    output <- data.table::fread(file, sep = "\t") # Load output
    ts <- find_ts(output) # Find time-to-stability (ts)
    return(list(id=id,ts=ts))
  })
  message("Ending worker ", core, "....\n")
  return(res)
}, mc.cores = num_cores)


# Create data frame with time to stability for each community
ts_df <- do.call(rbind, lapply(seq_along(tsa_list), function(i) {
  core_list <- tsa_list[[i]]
  do.call(rbind, lapply(core_list, function(e) {
    data.frame(simulation = e$id, t_stable = e$ts)
  }))
}))

# Save data frame
wd <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results"
path <- file.path(wd, exp_id, "tsa-table.tsv")
data.table::fwrite(x = ts_df, file = path, sep = "\t")
tsa <- max(ts_df$t_stable)
print(paste0("Time required for reaching stable state: ", tsa))
tictoc::toc() # For section 7
tictoc::toc() # For Total running time