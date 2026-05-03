# =============================================================================
# Impact quantification analysis full community vs sub-community.
# =============================================================================
#
# Interaction matrix setup:
#   - Diagonal:          fixed at -0.5
#   - Other off-diagonal: sparse mix of zeros and negatives from U(0,1)
#
# =============================================================================


#--------------------------------------------------------------------------
# Section: Generate-ID-and-paths
tictoc::tic("Section 1: Parameter generation")
start = Sys.time()

# Indicate directories paths
pdir <- "/mnt/data/sur/users/mrivera/Data/clean_controls"  # Parent-dir                                                
id <- substr(ids::uuid(1, drop_hyphens = TRUE, use_time = TRUE), start=1, stop=6)            
experiment_id <- paste0("Kboost_proof_", id) # Experiment-ID

# Generate directory paths
experiment_dir <- file.path(pdir, experiment_id)                    
cat(">> The experiment path is:", experiment_dir,"\n", sep=" ")
dir.create(experiment_dir)

 
# Section: Generate directories
dirs <- c('RawOutputs', 'Interactions', 'Topologies', 'Full_ExtSummaries', 'Sub_ExtSummaries')
sapply(file.path(experiment_dir, dirs), dir.create)


generate_params <- function(samples = 100) {
  dt <- data.table::CJ(n_species = 30, p_neg = 1, p_noint = seq(0, 0.9, by = 0.1))
  # Calculate how many times to replicate each row
  reps <- samples / nrow(dt)        
  dt <- dt[rep(1:.N, each = reps)]     # Replicate each row 'reps' times
  n_total   <- nrow(dt)
  all_seeds <- matrix(sample.int(3e6L, 3L * n_total, replace = FALSE), ncol = 3) 
  dt[, `:=`(
      key     = sample(x = 1:30, size = n_total, replace = TRUE),  # keystone specie
      id       = ids::random_id(n = n_total, bytes = 4),
      x0_seed  = all_seeds[, 1],                                                  
      mu_seed  = all_seeds[, 2],                                                  
      A_seed   = all_seeds[, 3]                                                   
  )]
}

df_params <- generate_params(samples = 1000)
# gen_Kboost_params(row = df_params[1, ])
# Verify if ids are unique and in case they are, save the parameters.
while (nrow(df_params) != length(unique(df_params$id))) {
  df_params <- generate_params(samples = 1000) # Repeat function
}
      
# Save parameters
params_path <- file.path(experiment_dir,"simulation_params.tsv")    # Parameters-TSV
data.table::fwrite(x = df_params, file = params_path, sep = "\t", quote = FALSE, row.names = FALSE) 
message("\nParameters generated and saved at path:\n", params_path, "\n")
cat(">> The number of extinctions to do is:", 30 * nrow(df_params),"\n", sep=" ")
tictoc::toc() # For section 1

#----------------------------------------------
# Section: Source functions
codes = list.files("/mnt/data/sur/users/mrivera/gLV/src-sims/FUN", full.names=TRUE)

invisible(lapply(codes, function(file){
  cat(">> Sourcing function: ", file, "\n")
  capture.output(source(file))
  return()
}))

#----------------------------------------------
# Section: Node statistics function
#----------------------------------------------
library(igraph)
build_topology <- function(A) {
    g <- graph_from_adjacency_matrix(abs(A), mode = "directed", weighted = TRUE)
    in_deg  <- apply(A, 2, function(x) c(in_pos  = sum(x > 0), in_neg  = sum(x < 0)))
    out_deg <- apply(A, 1, function(x) c(out_pos = sum(x > 0), out_neg = sum(x < 0)))
    top_df <- cbind(
        in_pos       = in_deg["in_pos", ],
        in_neg       = in_deg["in_neg", ],
        out_pos      = out_deg["out_pos", ],
        out_neg      = out_deg["out_neg", ],
        strength_in  = strength(g, mode = "in"),
        strength_out = strength(g, mode = "out"),
        betweenness  = betweenness(g, directed = TRUE),
        closeness    = closeness(g, mode = "out"),
        pagerank     = page_rank(g)$vector,
        transitivity = transitivity(g, type = "local"),
        in_coreness  = coreness(g, mode = "in"),
        out_coreness = coreness(g, mode = "out")
    )
    top_df[is.nan(top_df)] <- 0
    round(top_df, 3)
}

#----------------------------------------------
# Section: All extinctions function
#----------------------------------------------
sim_all_ext <- function(params) {
  # Unperturbed community populations
  x0      <- params$x0
  rel_x0  <- x0 / sum(x0)
  # Perturbate community
  dplyr::bind_rows(lapply(1:params$n, function(i) {
    # Remove species i from the parameters
    tmp_params <- list(
      x0 = x0[-i],
      mu = params$mu[-i],
      M  = params$M[-i, -i, drop = FALSE]
    )
    # Run simulation
    new_out     <- solve_gLV(times = 1000, tmp_params)  
    # population after extinction of species i
    x_after     <- new_out[, ncol(new_out)]            
    # population before extinction without species i
    x_before    <- tmp_params$x0                          
    # New extinctions after removing species i
    n_ext       <- sum(x_after < 1e-6 & x_before > 1e-6)
    bray_curtis <- 1 - (2 * sum(pmin(x_before, x_after))) / (sum(x_before) + sum(x_after))
    # return df
    data.frame(
      specie           = i,
      pop_initial      = x0[i],
      rel_pop_initial  = rel_x0[i],
      n_extinctions    = n_ext,
      prop_extinctions = n_ext / length(x_before),
      dissimilarity_bc = bray_curtis,
      keystoneness     = bray_curtis * (1 - rel_x0[i]),
      time_stability   = find_ts(new_out)
    )
  }))
}

#----------------------------------------------
# Section: Sub-community impact function
#----------------------------------------------
library(dplyr)
sim_sub_ext <- function(params, ext_threshold = 1e-06) {
  # Filter to surviving species
  rel_x0      <- params$x0 / sum(params$x0)
  to_filter   <- which(rel_x0 > ext_threshold)
  # Sub-community
  sub_x0  <- params$x0[to_filter]
  sub_mu  <- params$mu[to_filter]
  sub_M   <- params$M[to_filter, to_filter]
  rel_sub <- sub_x0 / sum(sub_x0)
  dplyr::bind_rows(lapply(seq_along(to_filter), function(i) {
      tmp_params <- list(
          x0 = sub_x0[-i],
          mu = sub_mu[-i],
          M  = sub_M[-i, -i, drop = FALSE]
      )
      new_out <- solve_gLV(times = 1000, tmp_params)
      x_after <- new_out[, ncol(new_out)] # final population after extinction of species i
      x_before <- tmp_params$x0           # population before extinction
      n_ext      <- sum(x_after < 1e-6)   # new extinctions after removing species i
      bray_curtis <- 1 - (2 * sum(pmin(x_before, x_after))) / (sum(x_before) + sum(x_after))

      data.frame(
          specie           = to_filter[i],
          pop_initial      = sub_x0[i],
          rel_pop_initial  = rel_sub[i],
          n_extinctions    = n_ext,
          prop_extinctions = n_ext / length(x_before),
          dissimilarity_bc = bray_curtis,
          keystoneness     = bray_curtis * (1 - rel_sub[i]),
          time_stability   = find_ts(new_out)
      )
  }))
}

#----------------------------------------------
# Section: Wrapper function
#----------------------------------------------
make_path <- function(experiment_dir, folder, prefix, sim_id) {
    file.path(experiment_dir, folder, paste0(prefix, sim_id, ".feather"))
}

save_feather <- function(x, path) arrow::write_feather(as.data.frame(x), sink = path)  # helper

wrapper <- function(index, df_params, ext_threshold, experiment_dir) {
    row    <- df_params[index, ]
    sim_id <- row$id
    # Simulate
    params         <- gen_Kboost_params(row)
    output         <- solve_gLV(times = 1000, params)
    output_subset  <- output[, c(1, seq(50, ncol(output), by = 50))]
    params$x0      <- output[, ncol(output)]                          # update to stable state
    # Paths list
    paths <- list(
        out      = make_path(experiment_dir, "RawOutputs",       "RawOutput_",  sim_id),
        A        = make_path(experiment_dir, "Interactions",      "A_",          sim_id),
        topology = make_path(experiment_dir, "Topologies",        "Topology_",   sim_id),
        full_ext = make_path(experiment_dir, "Full_ExtSummaries", "ExtSummary_", sim_id),
        sub_ext  = make_path(experiment_dir, "Sub_ExtSummaries",  "ExtSummary_", sim_id)
    )
    # Topology and extinction summaries
    topology_df  <- build_topology(params$M)
    summary_full <- sim_all_ext(params)
    summary_sub  <- sim_sub_ext(params, ext_threshold)
    # Save
    save_feather(output_subset,  paths$out)
    save_feather(params$M,       paths$A)
    save_feather(topology_df,    paths$topology)
    save_feather(summary_full,   paths$full_ext)
    save_feather(summary_sub,    paths$sub_ext)
    cat(">> Simulation", sim_id, "completed.\n")
    list(
        id            = sim_id,
        na_ct         = any(is.na(output)),
        tts_output    = find_ts(output),
        ext_threshold = ext_threshold,
        tts_ext       = max(summary_full$time_stability)
    )
}

# Test line:
# wrapper(index=1, df_params, ext_threshold= 1e-06, experiment_dir)
#----------------------------------------------
# Section: Parallelize-code and get the summary of simulations
tictoc::tic("Section 2: Run simulations and extinctions using the parallel package")

library(parallel)
ncore = max(1, detectCores() - 1, na.rm = TRUE)
cat('>> The number of cores to use are: ', ncore, '\n')
ext_threshold = 1e-06
cat('>> The extinction threshold is: ', ext_threshold, '\n')

results_summary <- mclapply(
  #seq_len(10),
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
arrow::write_feather(x = as.data.frame(result_df), sink = info_path)
tictoc::toc() 

#--------------------------------------------
# Section: Write TXT file
elapsed = Sys.time() - start
summary_text <- paste0(
  "Extinction impact analysis summary\n",
  "=========================\n",
  "Method used: Boost keystone column by K \n",
  "Time for solver: 1000 \n",
  "Number of cores: ", ncore, "\n",
  "Extinction threshold: ", ext_threshold, "\n",
  "Communities simulated: ", nrow(df_params), "\n",
  'Elapsed time: ', round(elapsed,5), ' minutes \n'
)

writeLines(summary_text, file.path(experiment_dir, "simulation_notes.txt"))