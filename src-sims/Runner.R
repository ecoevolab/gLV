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

#============================================================================
# SECTION: Generate-ID-and-paths
tictoc::tic("Section 0: Total running time")

#' Indicate directories paths
pdir <- "/mnt/data/sur/users/mrivera"                                     # Parent-dir
exp_id <- substr(ids::uuid(1, use_time = TRUE), 1, 13)                    # exp-id


# Verify-experiment-unique-id
res_dir <- file.path(pdir, "Train-sims")                                               
direx = basename(list.dirs(res_dir, recursive = FALSE))

while (exp_id %in% direx)  {
  exp_id <- substr(ids::uuid(1, use_time = TRUE), 1, 13) 
}                                           

exp_dir <- file.path(res_dir, exp_id)                                                   # Experiment-dir
params_path <- file.path(pdir, "Data", paste0(exp_id, ".tsv"))                          # Parameters-TSV
mc_dir <- file.path(res_dir, exp_id, "mc-apply")                                        # Workers-dir
info_path <- file.path(res_dir, exp_id, "raw-sims.feather")                             # Information-TSV
cat(">> The experiment path is:", exp_dir,"\n", sep="")

#============================================================================
#' First we generate the parameters for simulation:
#+ eval=FALSE
# SECTION: generate-paths
tictoc::tic("Section 1: Time for Parameter Generation")

#' Generate grid of parameters:
#+ eval=FALSE
generate_params <- function (){
  reps = 300                                               # Replicas-per-combination
  p_noint = seq(0, 1, by = 0.05)                          # interaction=0
  n_species = rep(c(20, 100), times = reps)               
  dt <- data.table::CJ(n_species, p_neg = 1, p_noint )

  # Add Columns with data table operator `:=`
  n_total <- nrow(dt)
  all_seeds <- sample.int(3e6L, 3L *n_total, replace = FALSE)
  
  dt[, `:=`(
    id = ids::random_id(n = n_total, bytes = 5),
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
cat(">> The number of simulations is:", nrow(params_df),"\n", sep=" ")
unique_species <- unique(params_df$n_species)
rows_per_species <- nrow(params_df) / length(unique_species)
result <- sum(unique_species * rows_per_species)
cat(">> The number of extinctions to do is:", result,"\n", sep=" ")
tictoc::toc() # For section 1

#============================================================================
# SECTION: Divide-data-into-chunks

#' We split the data into chunks of `n_cores`:
#+ eval=FALSE
tictoc::tic("Section 2: Divide data into chunks")

library(parallel)
num_cores <- parallel::detectCores() - 1  # Use one less than the total number of cores
cat("The number of cores that will be used are: ", num_cores, "\n")

split_table <- function(df, n_chunks) {
  split(df, cut(seq_len(nrow(df)), breaks = n_chunks, labels = FALSE))
}

# TEST
chunks <- split_table(params_df[1:27,], num_cores)
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
codes = list.files(file.path(pdir, "gLV/src/FUN" ), full.names=TRUE)

lapply(codes, function(file){
  cat(">> Sourcing function: ", file, "\n")
  capture.output(source(file))
  return()
})

#============================================================================
# SECTION: Wrap-functions
#' We wrap the code for parallelizing the simulations.
#+ eval=FALSE
# ==== Wrapper for running all required steps ====
wrapper <- function(index, path_core) {

  # Review: testing
  # path_core = workers_ODE[1]
  # index=params_df[1,]
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

  #=================== Extinctions ===================
  params$x0 = output[,1000]                                                  # Stable-population
  preds_df = sim_all_ext(params, path_core)                                  # Generate-extinctions
  tts_ext = max(preds_df$ext_ts)                                             # time-to-stability EXTINCTIONS
  preds_path <- paste0(path_core, "/tgt_", params$id, ".feather")          # Extinctions-paths
  arrow::write_feather(preds_df, preds_path) 
  
  return(list(id = params$id, na_ct = NA_count, tts_out = tts_out, tts_ext = tts_ext))
}

#============================================================================
# SECTION: Parallelize-code
#' We Parallelize the code and get the summary of simulations
#+ eval=FALSE
# ==== Parallelize it ====
tictoc::tic("Section 4: Run simulations and extinctions using the parallel package")

# review testing
# chunks <- split_table(params_df[1:27,], num_cores)

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



# REVIEW fixing
# params_df
# dat = readr::read_tsv("/mnt/data/sur/users/mrivera/Data/2263e52c-8384.tsv")
# sum(unique(dat$n_species) * (nrow(dat)/length(unique(dat$n_species))))                                            # number-extinctions-files TOTAL
# tgt_dir = "/mnt/data/sur/users/mrivera/Experiments/4018ff92-8377/"
# source_dir = "/mnt/data/sur/users/mrivera/Experiments/4018ff92-8377/mc-apply"