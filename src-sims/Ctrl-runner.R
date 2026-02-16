# This script is for running the simulations of the gLV model in parallel using multiple cores.
# It generates a unique experiment ID, creates necessary directories, generates parameters,
#
# It is purpose is to generate positive controls simulations.

#============================================================================
# SECTION: Generate-ID-and-paths
tictoc::tic("Section 0: Total running time")

#' Indicate directories paths
pdir <- "/mnt/data/sur/users/mrivera/Controls"      # Parent-dir                                                
exp_id = "NegCtrl-V1"                               # Experiment-ID      
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
  dt <- data.table::CJ(n_species = 30, p_neg = 1, p_noint = seq(0, 0.9, by = 0.1))
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

#------------------------------------------------------
# Function to generate positive controls: build_posctrl
# testing lines
# index = df[900,]
# params = build_params(index)
# path_core = workers_ODE[1]
# num_cores = 5
# chunks = split_table(df[1:5,], num_cores)
# workers_ODE <- create_dirs(mc_dir, num_cores)
#------------------------------------------------------

#============================================================================
# SECTION: Topolgy function
# Next line is for testing
# A = params$M

library(igraph)
build_topology <- function(A) {
  # Create graph from adjacency matrix
  g <- graph_from_adjacency_matrix(abs(A), mode="directed", weighted=TRUE)
  # Generate topology
  top_df <- cbind(
    # In-degree
    in_pos = apply(A, 2, function(x) sum(x>0)),  # 2 means columns
    in_neg = apply(A, 2, function(x) sum(x<0)),
    # Out-degree
    out_pos = apply(A, 1, function(x) sum(x>0)),  
    out_neg = apply(A, 1, function(x) sum(x<0)),  # 1 means rows
    # Total degree
    total_degree = degree(g, mode="all"),
    strength_in = strength(g, mode="in"),    # Weighted IN degree
    strength_out = strength(g, mode="out"),  # Weighted OUT degree
    betweenness = betweenness(g, directed=TRUE),
    closeness = closeness(g, mode="out"),
    pagerank = page_rank(g)$vector,               # Google page-rank
    transitivity = transitivity(g, type="local"), # Transitivity
    in_coreness = coreness(g, mode = 'in'),  # In coreness
    out_coreness = coreness(g, mode = 'out')  # out coreness
  )
  # Next line is for testing if total_degree is correct
  # rowSums(top_df[,1:4]) == top_df[,"total_degree"]
  #
  # Convert NaN to 0
  top_df[is.nan(top_df)] <- 0
  top_df = round(top_df,3)
  # Add species
  # result = cbind(species = paste0("Sp", 1:nrow(A)), df)
  return(top_df)
}

#----------------------------------------------
# Section: Wrapper function

wrapper <- function(index, path_core) {
  #-----------------------------
  # Section: Generate parameters and run simulation
  params <- build_params(index)              # Generate-parameters
  output <- solve_gLV(times = 1000, params)   # Run-simulation
  #-----------------------------
  # Section: Generate filenames and save files
  sim_id <- params$id
  out_path <- file.path(path_core, paste0("RawOutput_", sim_id, ".feather"))        # simulation output
  A_path <- file.path(path_core, paste0("A_", sim_id, ".feather"))                  # interactions matrix
  topology_path <- file.path(path_core, paste0("Topology_", sim_id, ".feather"))    # network topology
  preds_path <- file.path(path_core, paste0("ExtSummary_", sim_id, ".feather"))     # extinctions summary
  # Save files
  # df_subset <- df[, c(1, seq(10, ncol(df), by = 10))]       # Create df with every 10 cols
  arrow::write_feather(x = output, sink = out_path)                      # Save output
  arrow::write_feather(x = as.data.frame(params$M), sink = A_path)       # Save interactions matrix
  #-----------------------------
  # Section: Generate topology and summary of the simulation
  na_count <- sum(is.na(output))                              # simulation-NAs
  out_stability_time <- find_ts(output)                       # time-to-stability OUTPUT
  topology_df <- build_topology(params$M)                     # network topology
  arrow::write_feather(x = as.data.frame(topology_df), sink = topology_path) # Save topology
  cat(">> Simulation ", params$id, " completed.\n")   
  #-----------------------------
  # Section: Simulate extinctions and summary of extinctions
  params$x0 = output[[ncol(output)]]                            # Stable-population
  summary_exts = sim_all_ext(params, path_core)                 # Generate-extinctions
  extinction_stability_time = max(summary_exts$time_stability)  # time-to-stability 
  # Save files
  arrow::write_feather(x = summary_exts, sink = preds_path)   
  #-----------------------------
  return(list(id = params$id, na_ct = na_count, tts_out = out_stability_time, tts_ext = extinction_stability_time))
}

#----------------------------------------------
# Section: Parallelize-code and get the summary of simulations
tictoc::tic("Section 4: Run simulations and extinctions using the parallel package")

simulation_summary = parallel::mclapply(1:num_cores, function(core_id) {
  message("Starting worker ", core_id, "....\n")
  core_chunk <- chunks[[core_id]]  # rows assigned to this core
  #----------------------------
  result <- lapply(1:nrow(core_chunk), function(i) {
    wrapper(index = core_chunk[i, ], path_core = workers_ODE[core_id])
  })
  # Convert list to df
  result <- data.table::rbindlist(result, use.names = TRUE) 
  message("Ending worker ", core_id, "....\n")
  return(result)
}, mc.cores = num_cores)

# Generate TSV file
simulation_summary_df <- data.table::rbindlist(simulation_summary) # Convert list (of df) to df
arrow::write_feather(simulation_summary_df, info_path)
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
