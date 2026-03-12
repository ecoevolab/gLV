# This script is for running the simulations of the gLV model in parallel using multiple cores.
# It generates a unique experiment ID, creates necessary directories, generates parameters,
#
# Its purpose is to generate positive controls simulations.

#--------------------------------------------------------------------------
# Section: Generate-ID-and-paths
tictoc::tic("Section 0: Total running time")

# Parent directory
parent_directory <- "/mnt/data/sur/users/mrivera/Training-Data"  

# Generate experiment ID
experiment_id = paste0('Batch_', substr(ids::random_id(),1, 12))
exist = list.dirs(path = parent_directory, recursive = FALSE)
while (experiment_id %in% basename(exist)) {
  experiment_id = paste0('Batch_', substr(ids::random_id(),1, 12)) # Regenerate if ID already exists
}     

# Generate directory paths
experiment_dir <- file.path(parent_directory, experiment_id)                    # Experiment-dir
dir.create(experiment_dir)
cat(">> The experiment path is:", experiment_dir,"\n", sep=" ")
 
# Section: Generate directories
dirs <- c('RawOutputs', 'Interactions', 'Topologies', 'ExtSummaries')
sapply(file.path(experiment_dir, dirs), dir.create)

#--------------------------------------------------------------------------
# Section: generate parameters
generate_params <- function (n_reps = 100){
  dt <- data.table::CJ(n_species = 30, p_neg = 1, p_noint = seq(0, 0.9, by = 0.1))
  dt <- dt[rep(1:.N, each = n_reps)]     # Replicate each row 'reps' times
  n_total <- nrow(dt)
  all_seeds <- sample.int(3e6L, 3L *n_total, replace = FALSE)
  
  dt[, `:=`(
    key = sample(x = 1:30, size = n_total, replace = TRUE),  # key species to not go extinct
    id = ids::random_id(n = n_total, bytes = 4),
    x0_seed = all_seeds[1:n_total],
    mu_seed = all_seeds[(n_total+1):(n_total*2)],
    A_seed = all_seeds[(2 * n_total + 1):(3 * n_total)]
  )]
  
  return(dt)
}

n_reps = 1
df_params <- generate_params(n_reps)
# Verify if ids are unique and in case they are, save the parameters.
while (nrow(df_params) != length(unique(df_params$id))) {
  df_params <- generate_params(n_reps) # Repeat function
}
      
# Save parameters
params_path <- file.path(experiment_dir,"simulation-params.tsv")    # Parameters-TSV
data.table::fwrite(x = df_params, file = params_path, sep = "\t", quote = FALSE, row.names = FALSE) 
message("\nParameters generated and saved at path:\n", params_path, "\n")
cat(">> The number of extinctions to do is:", 30 * nrow(df_params),"\n", sep=" ")
tictoc::toc() # For section 1

#--------------------------------------------------------------------------
# Section: Source-codes
#' We source the function to generate the gLV parameters using the initial parameters and for solve it:

#+ eval=FALSE
# Source functions to:
codes = list.files("/mnt/data/sur/users/mrivera/gLV/src-sims/FUN", full.names=TRUE)

lapply(codes, function(file){
  cat(">> Sourcing function: ", file, "\n")
  capture.output(source(file))
  return()
})

#--------------------------------------------------------------------------
# Section: Node statistics function

library(igraph)
build_topology <- function(A) {
  # Create graph from adjacency matrix
  g <- graph_from_adjacency_matrix(abs(A), mode="directed", weighted=TRUE)
  #----------------------------
  # Section: Calculate degrees
  in_degrees <- apply(A, 2, function(x) c(in_pos = sum(x > 0), in_neg = sum(x < 0)))
  out_degrees <- apply(A, 1, function(x) c(out_pos = sum(x > 0), out_neg = sum(x < 0)))
  # Generate topology
  top_df <- cbind(
    in_pos  = in_degrees["in_pos", ],
    in_neg  = in_degrees["in_neg", ],
    out_pos = out_degrees["out_pos", ],
    out_neg = out_degrees["out_neg", ],
    # Total degree
    # total_degree = degree(g, mode="all"),
    strength_in = strength(g, mode="in"),    # Weighted IN degree
    strength_out = strength(g, mode="out"),  # Weighted OUT degree
    betweenness = betweenness(g, directed=TRUE),
    closeness = closeness(g, mode="out"),
    pagerank = page_rank(g)$vector,               # Google page-rank
    transitivity = transitivity(g, type="local"), # Transitivity
    in_coreness = coreness(g, mode = 'in'),  # In coreness
    out_coreness = coreness(g, mode = 'out')  # out coreness
  )
  # Convert NaN to 0
  top_df[is.nan(top_df)] <- 0
  top_df = round(top_df,3)
  return(top_df)
}

#----------------------------------------------
# Section: Wrapper function
wrapper <- function(index, df_params) {
  #-----------------------------
  # Section: Generate parameters and run simulation
  row = df_params[index, ] 
  sim_id = row$id
  # Generate parameters
  params <- gen_training_params(row)                              # Generate-parameters
  output <- solve_gLV(times = 1000, params)                       # Run-simulation
  output_subset <- output[, c(1, seq(50, ncol(output), by = 50))]  # Save output every 50 cols
  #-----------------------------
  # Section: Generate filenames and save files
  out_path <- file.path(experiment_dir, 'RawOutputs', paste0("RawOutput_", sim_id, ".feather"))        # simulation output
  A_path <- file.path(experiment_dir, 'Interactions', paste0("A_", sim_id, ".feather"))                  # interactions matrix
  topology_path <- file.path(experiment_dir, 'Topologies', paste0("Topology_", sim_id, ".feather"))    # network topology
  preds_path <- file.path(experiment_dir, 'ExtSummaries', paste0("ExtSummary_", sim_id, ".feather"))     # extinctions summary
  # Save files
  arrow::write_feather(x = output_subset, sink = out_path)             # Save output
  arrow::write_feather(x = as.data.frame(params$M), sink = A_path)     # Save interactions matrix
  #-----------------------------
  # Section: Generate topology and summary of the simulation
  na_count <- sum(is.na(output))                              # simulation-NAs
  out_stability_time <- find_ts(output)                       # time-to-stability OUTPUT
  topology_df <- build_topology(params$M)                     # network topology
  arrow::write_feather(x = as.data.frame(topology_df), sink = topology_path) # Save topology
  cat(">> Simulation ", sim_id, " completed.\n")   
  #-----------------------------
  # Section: Simulate extinctions and summary of extinctions
  params$x0 = output[[ncol(output)]]                            # Stable-population
  summary_exts = sim_all_ext(params, path_core=NULL)            # Generate-extinctions
  extinction_stability_time = max(summary_exts$time_stability)  # time-to-stability 
  # Save files
  arrow::write_feather(x = summary_exts, sink = preds_path)   
  #-----------------------------
  return(list(id = sim_id, na_ct = na_count, tts_out = out_stability_time, tts_ext = extinction_stability_time))
}

#----------------------------------------------
# Section: Parallelize-code and get the summary of simulations
tictoc::tic("Section 4: Run simulations and extinctions using the parallel package")

library(parallel)
results_summary <- mclapply(
    seq_len(nrow(df_params)),           # iterate over row indices
    wrapper,
    df_params,
    mc.cores = max(1, detectCores() - 1, na.rm = TRUE)
)

# Generate information file
result_df <- data.table::rbindlist(results_summary, use.names = TRUE) 
info_path <- file.path(experiment_dir, "summary.feather")           # Information-TSV
arrow::write_feather(x = result_df, sink = info_path)
tictoc::toc() # For section 4

# Sanity check
all_dirs = list.dirs(experiment_dir, recursive = TRUE)[-1]
for (dir in all_dirs) {
  x = list.files(dir)
  p = paste0('>> The number of files at ', dir, ' are: ', length(x))
  print(p)
}