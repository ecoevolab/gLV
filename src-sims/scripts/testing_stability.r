



# Loading data
exps = list.dirs('/mnt/data/sur/users/mrivera/Data/null_data', recursive = FALSE, full.names = TRUE)
lapply(exps, function(x) {
    dir.exists(paste0(x, '/RawOutputs/'))
})


# Generate stability analysis 
ts_analysis = function(output, threshold = 0.05, window = 10) {
    stable <- sapply(seq_len(n_cols - 1), function(i) {
        diff <- output[, i] - output[, i + 1]
        perc_change <- diff / output[, i]
        all(abs(perc_change) < threshold)
    })
    c = lapply(seq_len(length(stable)- window + 1), function(i) {
        flag = all(stable[i:(i+window-1)])
        return(flag)
    })
    if (any(c)) {
        print('>> Stability reached...')
        return(c)
    } else {
        print('>> Stability not reached...')
        return(c)
    }
}

# Plot results
 # Load files
files = list.files(paste0(exps[1], '/RawOutputs/'), full.names = TRUE)
output = as.matrix(arrow::read_feather(files[1]))
    