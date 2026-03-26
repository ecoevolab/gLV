# This script is for running the simulations of the gLV model in parallel using multiple cores.
# It generates a unique experiment ID, creates necessary directories, generates parameters,
#
# Its purpose is to generate positive controls simulations.

#--------------------------------------------------------------------------
# Section: Generate-ID-and-paths
tictoc::tic("Section 1: Parameter generation")
start = Sys.time()

# Indicate directories paths
pdir <- "/mnt/data/sur/users/mrivera/clean_controls"  # Parent-dir                                                
experiment_id <- substr(ids::uuid(1, drop_hyphens = TRUE, use_time = TRUE), start=1, stop=12)            

# Generate directory paths
experiment_dir <- file.path(pdir, experiment_id)                    # Experiment-dir
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

n_reps = 100
df_params <- generate_params(n_reps)
# Verify if ids are unique and in case they are, save the parameters.
while (nrow(df_params) != length(unique(df_params$id))) {
  df_params <- generate_params(n_reps) # Repeat function
}
      
# Save parameters
params_path <- file.path(experiment_dir,"simulation_params.tsv")    # Parameters-TSV
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

invisible(lapply(codes, function(file){
  cat(">> Sourcing function: ", file, "\n")
  capture.output(source(file))
  return()
}))

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
make_path <- function(experiment_dir, folder, prefix, sim_id) {
  file.path(experiment_dir, folder, paste0(prefix, sim_id, ".feather"))
}

wrapper <- function(index, df_params, ext_threshold, experiment_dir) {
  #-----------------------------
  # Section: Generate parameters and run simulation
  row = df_params[index, ] 
  sim_id = row$id
  # Generate parameters
  params <- gen_Kboost_params(row)                                # Generate-parameters
  output <- solve_gLV(times = 1000, params)                       # Run-simulation
  output_subset <- output[, c(1, seq(50, ncol(output), by = 50))]  # Save output every 50 cols
  # Population at quasi stable state
  final <- output[, ncol(output)]  
  #-----------------------------
  # Section: Generate filenames and save files
  out_path <- make_path(experiment_dir, 'RawOutputs', 'RawOutput_', sim_id)     # output
  A_path   <- make_path(experiment_dir, 'Interactions', 'A_', sim_id)           # interactions
  topology_path <- make_path(experiment_dir, 'Topologies', 'Topology_', sim_id) # node statistics
  # Summary extinctions full community impact
  summary_full_path <- make_path(experiment_dir, 'Full_ExtSummaries','ExtSummary_', sim_id)   
  # Save files
  arrow::write_feather(x = as.data.frame(output_subset), sink = out_path)             # Save output
  arrow::write_feather(x = as.data.frame(params$M), sink = A_path)     # Save interactions matrix
  #-----------------------------
  # Section: Generate topology and summary of the simulation
  na_count <- sum(is.na(output))                              # simulation-NAs
  out_stability_time <- find_ts(output)                       # time-to-stability OUTPUT
  topology_df <- build_topology(params$M)                     # node statistics
  arrow::write_feather(x = as.data.frame(topology_df), sink = topology_path) 
  #-----------------------------
  # Section: Simulate survivor nodes extinctions
  params$x0 = final                                             # Stable-population
  # Extinctions of survival nodes with impact in full community
  summary_full = sim_ext_surv_full(params, ext_threshold)    # Full community impact
  # If no extinctions were performed
  if (is.null(summary_full)){
    cat(">> Simulation ", sim_id, " completed.\n")
    return(list(id = sim_id, na_ct = na_count, tts_output = out_stability_time, ext_threshold = ext_threshold, ext_performed = FALSE, tts_ext = NA))
  }
  # Save extinctions at full community
  extinction_stability_time = max(summary_full$time_stability)  # time-to-stability 
  arrow::write_feather(x = as.data.frame(summary_full), sink = summary_full_path)   
  #-----------------------------
  cat(">> Simulation ", sim_id, " completed.\n")  
  return(list(id = sim_id, na_ct = na_count, tts_output = out_stability_time, ext_threshold = ext_threshold, ext_performed = TRUE, tts_ext = extinction_stability_time))
}

# Test line:
# wrapper(index=1, df_params, ext_threshold)
#----------------------------------------------
# Section: Parallelize-code and get the summary of simulations
tictoc::tic("Section 2: Run simulations and extinctions using the parallel package")

# Set extinction threshold 
ext_threshold = 5
cat('>> The extinctions threshold is of:', ext_threshold, '\n')

library(parallel)
ncore = max(1, detectCores() - 1, na.rm = TRUE)
cat('>> The number of cores to use are: ', ncore, '\n')
results_summary <- mclapply(
  seq_len(nrow(df_params)),
  wrapper,
  df_params       = df_params,
  ext_threshold   = ext_threshold,
  experiment_dir  = experiment_dir,
  mc.cores        = ncore
)

# Generate information file
result_df <- data.table::rbindlist(results_summary, use.names = TRUE) 
info_path <- file.path(experiment_dir, "simulation_summary.feather")           # Information-TSV
arrow::write_feather(x = as.data.frame(result_df,) sink = info_path)
tictoc::toc() 

#--------------------------------------------
# Section: Write TXT file
elapsed = Sys.time() - start
summary_text <- paste0(
  "Experiment summary\n",
  "=========================\n",
  "Method used: Boost keystone column by K \n",
  "Time for solver: 1000 \n",
  "Number of cores: ", ncore, "\n",
  "Extinction threshold: ", ext_threshold, "\n",
  "Communities simulated: ", nrow(df_params), "\n",
  "Communities with extinctions: ", sum(result_df$ext_performed), "\n",
  'Elapsed time: ', round(elapsed,5), ' minutes \n'
)

writeLines(summary_text, file.path(experiment_dir, "simulation_notes.txt"))