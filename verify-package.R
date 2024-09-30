devtools::document("~/Documents/LAB_ECO/gLV/glvpckg/R/")
devtools::install("~/Documents/LAB_ECO/gLV/glvpckg/R/")   # Reinstalls the package
devtools::load_all("~/Documents/LAB_ECO/gLV/glvpckg/R/")
devtools::check("~/Documents/LAB_ECO/gLV/glvpckg/R/")

# rm -rf /home/rivera/R/x86_64-pc-linux-gnu-library/4.4/00LOCK-glvsimulator
help(package="glvsimulator")


# Generate parameters for simulation
library(glvsimulator)

#-------------------------------------------------------------#
## Set the path to your directory
directory_path <- "~/Documents/LAB_ECO/gLV/glvpckg/R/"
r_files <- list.files(directory_path, pattern = "\\.R$", full.names = TRUE) # List all R files in the directory
lapply(r_files, source) # Source each file

# Generate seeds
# forge_seeds(n = 200, min = 2, max = 1000, wd = "~/Documents/LAB_ECO/")

wd = "~/Documents/LAB_ECO/Simulations"
# forge_seeds(n = 200, min = 2, max = 1000, wd)
seeds_path <- file.path(wd, "Seeds.tsv" )
params <- init_data(N_species = 2, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)

# Run simulation
times <- 100  # Define the number of generations
output <- run_simulation(N_species = 2, params = params, times = times)

# Generate unique ID
uniqueID <- forge_id(wd)

#----------------------- Savers----------------------------#

# Save output
output_saver(output, uniqueID, wd)

# Save parameters by seeds
params_seed_saver(N_species = 2,  C0 = 0.45, CN = 0.2, Diag_val = -0.5, params, uniqueID, wd)

# Save parameters by line
params_line_saver(params, uniqueID, wd)

#---------- Calculate all Steady States--------------------#
tolerance <- 0.005
SS_find_and_save_all(uniqueID, output, tolerance, wd)

# Apply Steady States Methods
result1 <- SS_roll_window_all(output, tolerance)
result2 <- SS_diff_means_all(output, tolerance)

#---------- Calculate individual Steady States-------------#
tolerance <- 0.05
individual_SS_find_and_save(uniqueID, output, tolerance, wd)


# Apply Steady States Methods
result1 <- all_rolling_var_SS(output, tolerance)
result2 <- all_prop_SS(output, tolerance)


#----------------------Forge tolerance---------------------#

wd = "~/Documents/LAB_ECO/Simulations"
# forge_seeds(n = 200, min = 2, max = 1000, wd)
seeds_path <- file.path(wd, "Seeds.tsv" )
params <- init_data(N_species = 2, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)

# Run simulation
times <- 100  # Define the number of generations
output <- run_simulation(N_species = 2, params = params, times = times)

# Generate unique ID
uniqueID <- forge_id(wd)


tolerance <- 0.0001
individual = TRUE

# Read the RDS file
wd <- "~/Documents/LAB_ECO/testing"

table <- load_individual_ss_rds(wd, uniqueID)

# Plot the output
out_path <- paste(wd, "/Outputs/O_", uniqueID , ".tsv", sep = "") # output path
output <- as.matrix( data.table::fread(out_path, sep = "\t") )
output_visualizer(output)
