# 11-March-2026
#
# Description: This script is to recalculate relative abundance.
# 

# Declare directories
outputs_dir = '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone/RawOutputs'

# Generate directories
base_dir <- '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone'

dirs <- c(
  filtered_ext  = file.path(base_dir,'Filter_ExtSummaries'),
  filtered_topo = file.path(base_dir,'Filter_Topologies'),
  filtered_A    = file.path(base_dir,'Filter_Interactions')
)

lapply(dirs, dir.create)

#-----------------------------------------------------
# Section: Declare function for node statistics
# Declare function for topology
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

library(tidyr)
library(dplyr)
new_wrapper = function(index) {
  #-------------------
  # Section: Read output and relative abundance
  row = parameters_tsv[index,]
  id = row$id
  path_out = paste0(outputs_dir, '/RawOutput_',id,'.feather')
  out = arrow::read_feather(path_out, col_select=21)[[1]]
  relative = out/sum(out)
  to_filter <- which(relative > 1e-06)
  #-------------------
  # Section: Generate parameters
  params = gen_Kboost_params(row)
  filter_params = list(x0 = out[to_filter], # for extinctions
    M = params$M[to_filter,to_filter],
    mu = params$mu[to_filter], 
    id = id, n = length(to_filter)
  )
  #----------------------
  # Section: Generate node statistics
  topology_df = build_topology(filter_params$M)
  path_topo = paste0(dirs[["filtered_topo"]], '/topology_',id,'.feather')
  arrow::write_feather(x = as.data.frame(topology_df), sink = path_topo)
  cat(paste0('>> Node statistics filtered for: ', id, '\n'))
  #---------------------
  # Section: Generate extinctions
  extinctions_result = sim_all_ext(filter_params)
  extinctions_path = paste0(dirs[["filtered_ext"]], '/ExtSummary_',id,'.feather')
  arrow::write_feather(x = extinctions_result, sink = extinctions_path) 
  cat(paste0('>> Extinctions filtered for: ', id, '\n'))
  #---------------------
  # Section: Save interactions
  A_path = paste0(dirs[["filtered_A"]], '/A_',id,'.feather')
  arrow::write_feather(x = as.data.frame(filter_params$M), sink = A_path) 
  cat(paste0('>> Interactions filtered for: ', id, '\n'))
}

#-------------------------------------------
codes = list.files("/mnt/data/sur/users/mrivera/gLV/src-sims/FUN", full.names=TRUE)

lapply(codes, function(file){
  cat(">> Sourcing function: ", file, "\n")
  capture.output(source(file))
  return()
})

# Load parameters file
parameters_path = '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone/simulation-params.tsv'
parameters_tsv = data.table::fread(parameters_path)

# Testing line
# row = parameters_tsv[1]
# test = new_wrapper(row)

library(parallel)
results <- mclapply(seq_len(nrow(parameters_tsv)), new_wrapper, mc.cores = detectCores() - 1)