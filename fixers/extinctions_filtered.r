# 11-March-2026
#
# Description: This script is to recalculate relative abundance.
# 

# Declare directories
ext_dir = '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone/ExtSummaries'
outs_dir = '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone/RawOutputs'
new_dir = '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone/Filter_ExtSummaries'
dir.create(new_dir)

# Simulation ids
path = '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone/simulation-params.tsv'
ids_vector = data.table::fread(path, select = 'id')[[1]]

library(tidyr)
library(dplyr)
new_wrapper = function(id) {
    #-------------------
    # Section: Load extinctions
    ext_path = paste0(ext_dir, '/ExtSummary_', id, '.feather')
    table = arrow::read_feather(ext_path)
    #-------------------
    # Section: Read output and relative abundance
    path_out = paste0(outs_dir, '/RawOutput_',id,'.feather')
    out = arrow::read_feather(path_out, col_select=21)[[1]]
    relative = out/sum(out)
    #----------------------
    # Section: Filter nodes
    to_filter <- which(relative > 1e-06)
    table_subset = table[to_filter, ] %>% mutate(rel_pop_initial = relative[to_filter]) 
    # proportion of extinctions cannot be changed as number of extinctions is overall...
    new_path = paste0(new_dir, '/ExtSummary_', id, '.feather')
    arrow::write_feather(x = table_subset, sink = new_path)
    cat(paste0('>> Extinctions generated for: ', id, '\n'))
}

library(parallel)
results <- mclapply(ids_vector, new_wrapper, mc.cores = detectCores() - 1)