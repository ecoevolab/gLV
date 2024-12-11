
## Set the path to your directory
directory_path <- "/mnt/atgc-d3/sur/users/mrivera/Functions"
r_files <- list.files(directory_path, pattern = "\\.R$", full.names = TRUE) # List all R files in the directory
lapply(r_files, source) # Source each file

# Generate seeds
# forge_seeds(n = 200, min = 2, max = 1000, wd = "~/Documents/LAB_ECO/")

#-----------------------------------------Usage example---------------------------------------------
# wd <- "/mnt/atgc-d3/sur/users/mrivera/Simulations"
wd <- "/mnt/atgc-d3/sur/users/mrivera/NR_sim"
seeds_path <- "/mnt/atgc-d3/sur/users/mrivera/combined_seeds.tsv"

forge_directories(wd)

for (s in 1:1000) {
  
  #-----------------------Generate parameters-----------------#
  C0 <- round( x = runif(1, min = 0, max = 1), digits = 3)
  CN <- round( x = runif(1, min = 0, max = 1), digits = 3)
  N_species <- round(x = runif(n = 1,min = 5, max = 100) )
  
  # Params for simulations
  params <-   forge_data(N_species, C0, CN, Diag_val = -0.5, seeds_path, wd)

  # Generate unique ID
  uniqueID <- forge_id(wd)
  
  #----------------------- Save parameters -----------------#
  # Save parameters by seeds
  params_seed_saver(N_species,  C0, CN, Diag_val = -0.5, params, uniqueID, wd)
  
  #-----------------------Run simulations-----------------#
  output <- run_simulation(N_species, params = params, times = 800)
  output_saver(output, uniqueID, wd) # Save output
}


