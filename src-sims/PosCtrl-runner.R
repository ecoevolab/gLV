# This script is for running the simulations of the gLV model in parallel using multiple cores.
# It generates a unique experiment ID, creates necessary directories, generates parameters,

#============================================================================
# SECTION: Generate-ID-and-paths
tictoc::tic("Section 0: Total running time")

#' Indicate directories paths
pdir <- "/mnt/data/sur/users/mrivera/Controls"      # Parent-dir                                                
exp_id = "PosCtrl-V1"                               # Experiment-ID      
exp_dir <- file.path(pdir, exp_id)                  # Experiment-dir

params_path <- file.path(exp_dir,"simulation-params.tsv")    # Parameters-TSV
mc_dir <- file.path(exp_dir, "mc-apply")                     # Workers-dir
info_path <- file.path(exp_dir, "summary.feather")           # Information-TSV
cat(">> The experiment path is:", exp_dir,"\n", sep=" ")
 

#============================================================================
#' First we generate the parameters for simulation:
#+ eval=FALSE
# SECTION: generate-paths
tictoc::tic("Section 1: Time for Parameter Generation")

#' Generate grid of parameters:
#+ eval=FALSE
generate_params <- function (){
  dt <- data.table::CJ(n_species = 30, p_neg = 1, p_noint = seq(0, 1, by = 0.1))
  dt <- dt[rep(1:.N, each = 100)]     # Replicate each row 'reps' times
  n_total <- nrow(dt)
  all_seeds <- sample.int(3e6L, 3L *n_total, replace = FALSE)
  
  dt[, `:=`(
    key = sample(x = 1:30, size = n_total, replace = TRUE),  # key species to not go extinct
    id = ids::random_id(n = n_total, bytes = 3),
    x0_seed = all_seeds[1:n_total],
    mu_seed = all_seeds[(n_total+1):(n_total*2)],
    A_seed = all_seeds[(2 * n_total + 1):(3 * n_total)]
  )]
  
  return(dt)
}

df <- generate_params()
# Verify if ids are unique and in case they are, save the parameters.
while (nrow(df) != length(unique(df$id))) {
    df <- generate_params() # Repeat function
}
      
# Save parameters
dir.create(exp_dir, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(x = df, file = params_path, sep = "\t", quote = FALSE, row.names = FALSE) 
message("\nParameteres generated and saved at path:\n", params_path, "\n")
cat(">> The number of extinctions to do is:", 30 * nrow(df),"\n", sep=" ")
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
chunks <- split_table(df, num_cores)
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
codes = list.files("/mnt/data/sur/users/mrivera/gLV/src-sims/FUN", full.names=TRUE)

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
library(igraph)
# install.packages("igraph", lib = "/mnt/data/sur/modules/pkgs/R_libs")

topology <- function(A) {
  # Create graph from adjacency matrix
  g <- graph_from_adjacency_matrix(abs(A), mode="directed", weighted=TRUE)
  # Generate topology
  df <- cbind(
    in_degree = degree(g, mode="in"),
    out_degree = degree(g, mode="out"),
    total_degree = degree(g, mode="all"),
    strength_in = strength(g, mode="in"),         # Weighted IN degree
    strength_out = strength(g, mode="out"),  # Weighted OUT degree
    betweenness = betweenness(g, directed=TRUE),
    closeness = closeness(g, mode="out")
  )
  # Convert NaN to 0
  df[is.nan(df)] <- 0
  df = round(df,3)
  # Add species
  # result = cbind(species = paste0("Sp", 1:nrow(A)), df)
  return(df)
}
A = params$M
wrapper <- function(index, path_core) {

  # Review: testing
  # path_core = workers_ODE[1]
  # index=df[1,]
  #=================== Output ===================
  params <- posctrl_regenerate(index)                                       # Generate-parameters
  output <- solve_gLV(times = 1000, params)                                 # Run-simulation
  
  # Generate filenames
  sim_id <- params$id
  out_path <- file.path(path_core, paste0("O_", sim_id, ".feather"))     # Simulation-path
  A_path <- file.path(path_core, paste0("A_", sim_id, ".feather"))       # Mat-path
  preds_path <- file.path(path_core, paste0("tgt_", sim_id, ".feather")) # Extinctions-paths

  arrow::write_feather(x = output, sink = out_path)                         # Saving-ODE
  arrow::write_feather(x = as.data.frame(params$M), sink = A_path)          # Save-MAT
  
  #=================== Information ===================
  NA_count <- sum(is.na(output))                            # simulation-NAs
  tts_out <- find_ts(output)                                # time-to-stability OUTPUT
  cat(">> Simulation ", params$id, " completed.\n")

  # Generate topology
  topo_df <- topology(params$M)
 

  #=================== Extinctions ===================
  params$x0 = output[[ncol(output)]]                    # Stable-population
  summary_exts = sim_all_ext(params, path_core)         # Generate-extinctions
  tts_ext = max(summary_exts$ext_ts)                    # time-to-stability EXTINCTIONS

  arrow::write_feather(x = summary_exts, sink = preds_path) 
  arrow::read_feather(preds_path)
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