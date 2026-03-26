
#' Description:
#' The script is to filter the communities where at least 5 species survived.
#' Within those filtered communities we will filter the extinctions of only the survival nodes.
#' 

# Load data
training_dir = '/mnt/data/sur/users/mrivera/Training-Data'
experiments_list = list.dirs(training_dir, recursive = FALSE)

# Function to filter one simulation extinctions 
filter_wrapper = function(index, ids_list, output_dir, extinctions_dir, filtered_dir){
    #---------------------
    # Load output
    id = ids_list[index]
    output_path = file.path(output_dir, paste0('RawOutput_', id, '.feather'))
    # Last time saved is 951 which is column 21
    out = arrow::read_feather(output_path, col_select = 21)
    out_relative = out/sum(out)
    to_filter = out_relative > 1e-06
    n_filtered = sum(to_filter)
    if (!(n_filtered)) {
        cat('>> Skipping filtering of simulation', id, 'no species survived\n')
        return(list(id = id, survivors = 0))
    }
    #---------------------
    # Section: Filter extinctions
    ext_df = arrow::read_feather(file.path(extinctions_dir, paste0('ExtSummary_', id, '.feather')))
    filtered_df = ext_df[to_filter,]
    # Save filtered extinctions
    filter_path = file.path(filtered_dir, paste0('Filtered_ExtSummary_', id, '.feather'))
    arrow::write_feather(x = as.data.frame(filtered_df), sink = filter_path)
    return(list(id = id, survivors = sum(to_filter)))
}

library(parallel)
for (i in 1:length(experiments_list)) {
    #---------------------
    # Declare experiment name
    experiment = experiments_list[i]
    cat('>> Starting filtering of simulation ', experiment, ' \n')
    # Directory to save filtered extinctions
    filtered_exts_dir = file.path(experiment, 'Filtered_ExtSummaries')
    dir.create(filtered_exts_dir, showWarnings = FALSE)
    # Load parameters
    ids_list = data.table::fread(file.path(experiment, 'simulation-params.tsv'))[['id']]
    # id = ids_list[1] # Testing line
    # Output paths
    out_dir = file.path(experiment, 'RawOutputs')
    ext_dir = file.path(experiment, 'ExtSummaries')
    #---------------------
    results = mclapply(
        X        = seq_along(ids_list),
        FUN      = filter_wrapper,
        ids_list = ids_list,
        output_dir      = out_dir,
        extinctions_dir = ext_dir,
        filtered_dir    = filtered_exts_dir,
        mc.cores = parallel::detectCores() - 1  # leave 1 core free
    )
    # Convert to a data frame
    results_df <- data.table::rbindlist(results, fill = TRUE)
    save_path = file.path(experiment, 'FilterSummary.feather')
    arrow::write_feather(x= as.data.frame(results_df), sink = save_path)
    cat('>> Finished filtering of simulation ', experiment, ' \n')
}


data1 = arrow::read_feather('/mnt/data/sur/users/mrivera/Training-Data/Batch_260225_gibab/FilterSummary.feather')
head(data1)
data2 = arrow::read_feather('/mnt/data/sur/users/mrivera/Training-Data/Batch_260227_lupit/FilterSummary.feather')
head(data2)