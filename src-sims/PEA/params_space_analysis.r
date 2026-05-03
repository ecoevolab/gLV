#  Parameter Exploratory analysis
#
# Description: 
# This script is for recreating the parameter space we can simulate.

#------------------------------------------
# Generate directories
#------------------------------------------
tictoc::tic("Section 0: Total running time")

tictoc::tic("Section 1: Time for Parameter Generation")
#' Indicate directories paths
result_dir <- "/mnt/data/sur/users/mrivera/Data/PEA"    # directory to save      
if (!dir.exists(result_dir)) {
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE) # Create directory if it doesn't exist
}        
params_path <- file.path(result_dir, "simulation_params.tsv")     # Parameters-TSV

cat(">> The experiment path is:", result_dir,"\n", sep="")

#------------------------------------------
# Section: Generate data with different interactions
#-------------------------------------------
library(dplyr)
generate_params <- function (){
  p_neg = seq(0, 1, by = 0.1)       # negative-interactions
  p_noint = seq(0, .9, by = 0.1)     # null-interactions                    
  n_species = rep(c(10,30,50,100), each = 10) # number of species               
  dt <- data.table::CJ(n_species, p_neg, p_noint ) 
  dt[, p_pos := 1 - p_neg]

  # Add Columns with data table operator `:=`
  n_total <- nrow(dt)
  all_seeds <- sample.int(3e6L, 3L *n_total, replace = FALSE)
  
  dt[, `:=`(
    id = ids::random_id(n = n_total, bytes = 3),
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
cat(">> The number of simulations are:", nrow(params_df),"\n", sep=" ")
tictoc::toc() # For section 1

#------------------------------------------
# Section: Declare functions
#-------------------------------------------
# Source functions to:
codes = list.files("/mnt/data/sur/users/mrivera/gLV/src-sims/FUN", full.names=TRUE)

invisible(lapply(codes, function(file) {
  cat(">> Sourcing function:", file, "\n")
  capture.output(source(file))
  # No need for return() here
}))

#------------------------------------------
# Section: Run simulations
#-------------------------------------------
library(arrow)
library(data.table)
library(parallel)

tictoc::tic("Section 2: Run simulations and extinctions using the parallel package")
wrapper <- function(index, df_params) {
  #-----------------------------
  # Section: Generate parameters and run simulation
  row = df_params[index, ] 
  sim_id = row$id
  # Generate parameters
  params <- gen_training_params(row)            # Generate-parameters
  output <- solve_gLV(times = 1000, params)     # Run-simulation
  na_count <- any(is.na(output))                     # simulation-NAs
  #-----------------------------
  cat(">> Simulation ", sim_id, " completed.\n")  
  return(list(id = sim_id, na_count = na_count))
}

# Assign cores
ncore = max(1, detectCores() - 1, na.rm = TRUE)
cat('>> The number of cores to use are: ', ncore, '\n')

# Parallelize
results_summary <- mclapply(
  seq_len(nrow(params_df)),
  # seq_len(200),
  wrapper,
  df_params       = params_df,
  mc.cores        = ncore
)

# Generate information file
result_df <- data.table::rbindlist(results_summary, use.names = TRUE)
info_path <- file.path(result_dir, "params_space_smry.feather")  # Information-TSV
arrow::write_feather(x = as.data.frame(result_df), sink = info_path)
tictoc::toc()