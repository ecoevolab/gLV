devtools::document("~/Documents/LAB_ECO/gLV/glvpckg/R/")
devtools::load_all("~/Documents/LAB_ECO/gLV/glvpckg/R/")
devtools::check("~/Documents/LAB_ECO/gLV/glvpckg/R/")

# Generate parameters for simulation
library(glvsimulator)

wd = "~/Documents/LAB_ECO/"
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
tolerance <- 0.05
SS_find_and_save_all(uniqueID, output, tolerance, wd)

# Apply Steady States Methods
result1 <- SS_roll_window_all(output, tolerance)
result2 <- SS_diff_means_all(output, tolerance)

#---------- Calculate individual Steady States-------------#
tolerance <- 0.05
individual_SS_find_and_save(uniqueID, output, tolerance, wd)


# Apply Steady States Methods
result2 <- individual_diff_SS(uniqueID, output, tolerance, wd)
result1 <- individual_rolling_window_SS(uniqueID, output, tolerance, wd)


#----------------------Testing------------------------------------
# Set the path to your directory
directory_path <- "~/Documents/LAB_ECO/gLV/glvpckg/R/"

# List all R files in the directory
r_files <- list.files(directory_path, pattern = "\\.R$", full.names = TRUE)

# Source each file
lapply(r_files, source)

# Example usage:
wd <- "~/Documents/LAB_ECO/testing"
seeds_path <- file.path("~/Documents/LAB_ECO", "Seeds.tsv")

paramSettings <- list(N_species = 2,
                      C0 = 0.45,
                      CN = 0.2,
                      Diag_val = -0.5, # Parameters for data generation
                      tolerance = 0.05, # Tolerance used for Steady State search
                      times = 100) # Number of generations
individual = TRUE

# Read the RDS file
uniqueID <- forge_tolerance(paramSettings, seeds_path, wd, individual)
table <- load_individual_ss_rds(wd, uniqueID)

# Plot the output
out_path <- paste(wd, "/Outputs/O_", uniqueID , ".tsv", sep = "") # output path
output <- as.matrix( data.table::fread(out_path, sep = "\t") )
output_visualizer(output)
