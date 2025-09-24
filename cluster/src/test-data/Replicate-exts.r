
# This code is for redoing the extinctions whit the extinction coefficient fixed. 

# Section: Load-data
tictoc::tic("Section 1: Divide data into chunks")

# Load parameters table
params_df = read.csv("/mnt/data/sur/users/mrivera/Data/c748247a-8dc2.tsv", sep = "\t")

library(parallel)
num_cores <- parallel::detectCores() - 2  # Use one less than the total number of cores
cat("The number of cores that will be used are: ", num_cores, "\n")

split_table <- function(df, n_chunks) {
  split(df, cut(seq_len(nrow(df)), breaks = n_chunks, labels = FALSE))
}

chunks <- split_table(params_df, num_cores)
message("\n Section: Data split completed.\n")
tictoc::toc() # For section 1

# Section: Generate-new-workers
tictoc::tic("Section 2: Generate directories for each core")

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
# Create new-workers-directory
wd = "/mnt/data/sur/users/mrivera/Experiments/c748247a-8dc2/Replica2"
mc_dir = file.path(wd, "mc-exts")
workers_ODE <- create_dirs(mc_dir, num_cores)
message("\nWorking directories created at path:\n", mc_dir,"\n")
tictoc::toc() # For section 3

#============================================================================
# SECTION: Source-codes
#' We source the function to generate the gLV parameters using the initial parameters and for solve it:

#+ eval=FALSE
# Source functions to:
codes = list.files("/mnt/data/sur/users/mrivera/gLV/src/FUN", full.names=TRUE)

lapply(codes, function(file){
  cat(">> Sourcing function: ", file, "\n")
  capture.output(source(file))
  return()
})

#============================================================================
# SECTION: Wrap-functions
# We wrap the code for parallelizing the simulations.
wrapper <- function(index, path_core) {

    # Review: testing
    # path_core = workers_ODE[1]
    # index=params_df[1,]
    #=================== Output ===================
    params <- regenerate(index)                                               # Generate-parameters
    Opath = file.path("/mnt/data/sur/users/mrivera/Experiments/c748247a-8dc2/raw-ODEs", paste0("O_", params$id, ".feather"))
    x = arrow::read_feather(Opath, col_select = 1000)
    params$x0 = unlist(x, use.names = FALSE)
    #=================== Extinctions ===================
    preds_df = sim_all_ext(params, path_core)                                       # Generate-extinctions
    tts_ext = max(preds_df$ext_ts)                                                  # time-to-stability EXTINCTIONS
    preds_path <- paste0(path_core, "/tgt_", params$id, ".feather")            # Extinctions-targets-for-GCN
    arrow::write_feather(preds_df, preds_path) 
    return(list(id = params$id, tts_ext = tts_ext))
}


#============================================================================
# SECTION: Parallelize-code
#' We Parallelize the code and get the summary of simulations
tictoc::tic("Section 3: Run simulations and extinctions using the parallel package")

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
info_path = file.path(wd, "RepExts-info.feather")
info_df <- data.table::rbindlist(sims_info) # Convert list (of df) to df
arrow::write_feather(info_df, info_path)
tictoc::toc() # For section 4


#============================================================================
# SECTION: Create-symbolic-links
#' We create symbolic links of the simulation...
#+ eval=FALSE
tictoc::tic("Section 4: Generate symbolic links")
mc_dir = file.path(wd, "mc-exts")
gen_syml(mc_dir)
tictoc::toc() # For section 5

tictoc::toc() # For Total running time
