
# This code is for checking if Outputs without a perturbation
# are < 1e-06. I will count the number of species with this value per simulation.

dir = "/mnt/data/sur/users/mrivera/Experiments/c748247a-8dc2/raw-ODEs"
files = list.files(dir, full.names = TRUE)

# Testing:
file = files[1]
df =  data.frame(id = character(0), value = numeric(0), less_0 = numeric(0) )
process_file <- function(file) {
    x <- arrow::read_feather(file, col_select = 1000)
    y <- sum(x < 1e-06)
    less0 <- sum(x < 0)
    data.frame(
        id = gsub("O_(.*)\\.feather", "\\1", basename(file)), 
        value = y, 
        less_0 = less0
    )
}

# Parallelize (uses all available cores by default)
library(parallel)
results <- parallel::mclapply(files, process_file, mc.cores = parallel::detectCores())
df <- do.call(rbind, results)