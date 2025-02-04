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

#-----------------------------------------Usage example---------------------------------------------
wd = "~/Documents/LAB_ECO/Simulations"

#-----------------------Generate parameters-----------------#
# forge_seeds(n = 200, min = 2, max = 1000, wd)
C0 <- round( x = runif(1, min = 0, max = 1), digits = 3)
CN <- round( x = runif(1, min = 0, max = 1), digits = 3)

#------------------------------Parameters------------------------------------
#
        #-----------------------Generate parameters-----------------#
        # seeds_path <- "/mnt/atgc-d3/sur/users/mrivera/Simulations/combined_seeds.tsv"
        seeds_path <- file.path(wd, "Seeds.tsv" )
        params <- forge_data(N_species = 20, C0 , CN, Diag_val = -0.5, seeds_path, wd)
        
        # Generate unique ID
        uniqueID <- forge_id(wd)
        
        #----------------------- Save parameters -----------------#
        # Save parameters by seeds
        params_seed_saver(N_species = 20,  C0, CN, Diag_val = -0.5, params, uniqueID, wd)
        
        # Save parameters by line
        params_line_saver(params, uniqueID, wd)
        
        #------------------------Read parameters -----------------#
        # Read parameters by seeds
        params <- params_seed_reader(uniqueID, wd)
        
        # Read parameters by lines
        params <- params_line_reader(uniqueID, wd)
        
#-----------------------Run simulations-----------------#
# Run simulation
output <- run_simulation(N_species = 20, params = params, times = 200)
output_saver(output, uniqueID, wd) # Save output


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
result1 <- individual_prop_SS(uniqueID, output, tolerance = 0.05, wd)
result2 <- individual_raw_diff_SS(uniqueID, output, tolerance = 0.05, wd)

#----------------------Forge tolerance---------------------#



#----------------------- Visualizers------------------------#
# wd = "~/Documents/LAB_ECO/Simulations"

wd <-  "/home/rivera/Cluster/Simulations"
uniqueID <- "fe9919"
uniqueID = "7e52e2"
out_path <- paste(wd, "/Outputs/O_", uniqueID , ".tsv", sep = "") # output path
output <- as.matrix( data.table::fread(out_path, sep = "\t") )

visualize_ld(wd)
visualize_output(output)

#------------------------Extinction-------------------------#
uniqueID = "02c234"
uniqueID = "128932"
wd = "~/Documents/LAB_ECO/Simulations"
out_path <- paste(wd, "/Outputs/O_", uniqueID , ".tsv", sep = "") # output path
output <- as.matrix( data.table::fread(out_path, sep = "\t") )
tmp = run_extinction(nspecs_extinct = 2, output, tolerance = 0.0005, uniqueID, wd)
output_visualizer(tmp)




seeds_path <- file.path("~/Documents/LAB_ECO", "Seeds.tsv")
params_tmp <- init_data(N_species = 2, seeds_path = seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)

