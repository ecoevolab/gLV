# 11-March-2026
#
# Description: This script is to redo extinctions in case of erasing by error.
# 


# Load dataset
experiment_dir = '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone'
parmeters_path = paste0(experiment_dir, '/simulation-params.tsv')
df_params = data.table::fread(parmeters_path)

# Directories 
output_dir =  paste0(experiment_dir, '/RawOutputs')
extinctions_dir = paste0(experiment_dir, '/ExtSummaries') 

#-----------------------------
# Section: Load functions
codes = list.files("/mnt/data/sur/users/mrivera/gLV/src-sims/FUN", full.names=TRUE)

lapply(codes, function(file){
  cat(">> Sourcing function: ", file, "\n")
  capture.output(source(file))
  return()
})

# Wrap code 
wrapper <- function(index, path_core = NULL) {
    row = df_params[index, ]  
    # Section: Generate parameters and run simulation
    id = row$id
    path_out = paste0(output_dir, '/RawOutput_', id, '.feather')
    out = arrow::read_feather(path_out, col_select = 21)[[1]]
    # Section: Generate extinctions
    params = gen_Kboost_params(row)
    params$x0 = out
    summary_exts = sim_all_ext(params, path_core = path_core)
    preds_path <- file.path(extinctions_dir, paste0("ExtSummary_", id, ".feather"))
    arrow::write_feather(x = summary_exts, sink = preds_path)
    # cat(paste0('>> Extinctions generated for: ', id, '\n'))
}

library(parallel)
results <- mclapply(
    seq_len(nrow(df_params)),           # iterate over row indices
    wrapper,
    path_core = NULL,       # passed as additional argument
    mc.cores = detectCores() - 1
)