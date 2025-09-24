# Code for generating Interactions matrices

# Paths for parameters and Interactions directory
par_df = data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Exp06-D29-Apr.tsv", sep="\t")
inter_dir = "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Exp06-D29-Apr/Interacts"

source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/src/generate-params.R") # regenerate parameters




# Wrapper for one row
process_row <- function(i) {
  params <- regenerate(i)          # Generate parameters
  M <- params$M                    # Interaction matrix
  path <- file.path(inter_dir, paste0("Intrs_", params$id, ".feather"))
  arrow::write_feather(as.data.frame(M), path)
  return(NULL)
}

# Apply the function in parallel over rows
n_cores <- parallel::detectCores() - 1
parallel::mclapply(1:nrow(par_df), function(i) process_row(par_df[i, ]), mc.cores = n_cores)


tictoc::tic("Starting interactions generation...")
# Step 1: Run simulations in parallel
results <- parallel::mclapply(1:nrow(par_df), function(i) {
  params <- regenerate(par_df[i, ])
  list(id = params$id, M = params$M)
}, mc.cores = n_cores)

# Step 2: Save in sequence
for (res in results) {
  path <- file.path(inter_dir, paste0("Intrs_", res$id, ".feather"))
  cat(path, "\n")
  arrow::write_feather(as.data.frame(res$M), path)
}
tictoc::toc()