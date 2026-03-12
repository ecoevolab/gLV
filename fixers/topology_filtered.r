# 11-March-2026
#
# Description: This script is to generate new topology filtering extinct nodes.
# 

# Declare directories
A_dir = '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone/Interactions'
outs_dir = '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone/RawOutputs'
new_dir = '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone/Filter_topologies'
dir.create(new_dir)

# Simulation ids
path = '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone/simulation-params.tsv'
ids_vector = data.table::fread(path, select = 'id')[[1]]

# Declare function for topology
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
  # Convert NaN to 0
  top_df[is.nan(top_df)] <- 0
  top_df = round(top_df,3)
  return(top_df)
}

# Function to generate new topology
new_wrapper = function(id){
    #----------------------
    # Load output and relative abundance
    path_out = paste0(outs_dir, '/RawOutput_',id,'.feather')
    out = arrow::read_feather(path_out, col_select=21)[[1]]
    relative = out/sum(out)
    #----------------------
    # Section: Filter nodes
    to_filter <- which(relative > 1e-06)
    # Load interactions
    path_A = paste0(A_dir, '/A_',id,'.feather')
    A = arrow::read_feather(path_A)
    # Filter extinct nodes
    A_subset = A[to_filter,to_filter]
    topology_df = build_topology(as.matrix(A_subset))
    #---------------------
    # Section: Save new topologies
    path_topo = paste0(new_dir, '/topology_',id,'.feather')
    arrow::write_feather(x = as.data.frame(topology_df), sink = path_topo)
    cat(paste0('>> New topology generated for: ', id, '\n'))
}
library(parallel)
results <- mclapply(ids_vector, new_wrapper, mc.cores = detectCores() - 1)

new_wrapper(id)
id = ids_vector[1]