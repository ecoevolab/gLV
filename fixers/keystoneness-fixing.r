# 2-Mrach-2026
# 
# This code is for fixing the keyetoneness of some simulations were the index was wrong.

# For training data
train_dirs =  list.dirs('/mnt/data/sur/users/mrivera/Training-Data', recursive = FALSE)

process_file <- function(f) {
  table <- arrow::read_feather(f) %>%
    mutate(
      keystoneness     = rel_pop_initial * dissimilarity_bc,
      prop_extinctions = n_extinctions / (n() - 1)
    )
  arrow::write_feather(x = table, sink = f)
  if (!anyNA(table)) {
    cat('>> file', basename(f), 'completed...\n')
  } else {
    cat('>> file', basename(f), 'has NAs!\n')
  }
  invisible(NULL)
}

all_files <- unlist(lapply(train_dirs, function(b) {
  list.files(paste0(b, '/ExtSummaries'), full.names = TRUE)
}))


library(parallel)
mclapply(all_files, process_file, mc.cores = parallel::detectCores() - 1)