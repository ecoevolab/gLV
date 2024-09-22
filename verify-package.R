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
individual_diff_SS(uniqueID, output, tolerance, wd)
