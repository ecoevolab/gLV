# This code is stored on the cluster

# Set the path to your directory
functions_path <- "/mnt/atgc-d3/sur/users/mrivera/Functions"
# functions_path <- "/home/rivera/Cluster/Functions/"

# List all R files in the directory and source them
r_files <- list.files(functions_path, pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_files, source))

# Generate seeds
wd <- "/mnt/atgc-d3/sur/users/mrivera/Simulations"
# wd <-  "/home/rivera/Cluster/Simulations"
forge_seeds(10000, min = 1 , max = 20000, wd )
forge_seeds(10000, min = 20001 , max = 40000, wd )

