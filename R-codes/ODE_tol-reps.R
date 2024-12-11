

#--------------------------------------Load functions-----------------------------------#
# functions_path <- "/home/rivera/Cluster/Functions"
functions_path <- "/mnt/atgc-d3/sur/users/mrivera/Functions"
  
r_files <- list.files(functions_path, pattern = "\\.R$", full.names = TRUE) # List all R files in the directory
# Source each file without printing output
invisible( lapply(r_files, function(file) {
  capture.output(source(file)) })
)

#-----------------------------------Get valid simulations------------------------------#
# Set working directory
# wd <- "/home/rivera/Cluster/testing"
wd <- "/mnt/atgc-d3/sur/users/mrivera/testing"

library(data.table)

# Filter simulations
diff_path <- file.path(wd, "Differences", "means_ld.tsv" )
table <- fread(diff_path, sep = "\t", header = TRUE, colClasses = c(ID = "character"))

# Use `complete.cases()` to filter rows directly and get valid IDs
controls <- table[complete.cases(table), ]

# Get the simulations we will repeat
raw_reps_ids <- setdiff(table$ID, controls$ID)
control_ids <- controls$ID

# Get IDs of the simulations that have been already done 
# rds_path <- file.path(wd, "Rp_tol", "ids-completed.rds")
# loaded_vector <- readRDS("my_vector.rds")

#-----------------------------Extract simulations to repeat----------------------------#
# RDS file path
rds_path <- file.path(wd, "test02", "completed_ids.rds")

# Load existing data if the file exists, otherwise use an empty vector
data <- if (file.exists(rds_path)) readRDS(rds_path) else c()

# Remove IDs that have already been completed
reps_ids <- setdiff(raw_reps_ids, data)

#------------------------------------Regenerate simulations----------------------------#

requireNamespace("deSolve")

ode_function <- function (times, params, atol, rtol) {
  
  # Define the equation
  glv_model <- function(t, x, params) {
    r <- params$Growths         # Growth rate vector
    A <- params$Interactions          # Interaction matrix
    
    # Compute dx/dt for each species
    dx <- x * (r + A %*% x)
    list(dx)
  }
  
  time_seq <- seq(0, times, by = 1)  # Define the time sequence
  
  # Get solution
  results <- deSolve::ode(y = params$Population, times = time_seq, func = glv_model, parms = params, method = "ode45",
                          rtol = rtol, 
                          atol = atol)
  
  # Remove the column of times and get the transversal, where rows are species and column generations
  ode_df <- as.matrix(t(results))[-1, -ncol(results)]
  
  # Get the proportion of each specie by column value/columnsum
  for (c in 1:ncol(ode_df)) {
    csum <- sum(ode_df[,c])
    ode_df[ , c] =  ode_df[ , c] / csum
    
  }
  
  return(ode_df)
}


#-------------------------------------------Simulations to repeat-----------------------------#

rtol_values <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1) # Relative tolerance
atol_values <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1) # Absolute tolerance

reps_fun <- function(s, wd, worker_path, rtol_values, atol_values, params_seed_reader) {
  
  # Read parameters
  params <- params_seed_reader(uniqueID = s, wd)
  
  # Loop in reverse order
  for (r in rev(rtol_values)) {
     for (a in rev(atol_values)) {
       
      # Run simulation
      new_out <- ode_function(times = 700, params, a, r)
      
      # Define paths
      tmp_path <- file.path(worker_path, paste0("tol_r", format(r, scientific = TRUE, digits = 3) ), 
                            paste0("tol_a", format(a, scientific = TRUE, digits = 3) ) )
      # Create the main directory if it doesn't exist
      if (!dir.exists(tmp_path)) {
        dir.create(tmp_path, recursive = TRUE)
      }
      
      # Save simulation
      save_path <- file.path(tmp_path, paste0("id_", s, ".tsv")  )
      utils::write.table(new_out, file = save_path, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
  }
}

#-------------------------------------------Separate IDs by core------------------------------#

library(parallel)

num_cores <- detectCores() - 1  # Use one less than the total number of cores
cat("The number of cores that will be used are: ", num_cores, "\n")

# Function to split IDs into chunks
split_ids <- function(ids, n_chunks) {
  split(ids, cut(seq_along(ids), breaks = n_chunks, labels = FALSE))
}

# For testing purposes we will sample 1000 from the pool:
# reps_ids <- sample(x = reps_ids, size = 5000)
cat("The number of IDs that will be repeated are: ", length(reps_ids), "\n")

# Split the IDs into 10 chunks
id_chunks <- split_ids(reps_ids, num_cores)
# total_ids <- sum(sapply(id_chunks, length))

#-----------------------------------------Parallelizing-------------------------------------------#

# Path to the main directory
main_dir <- file.path(wd, "test02")

# Create the main directory if it doesn't exist
if (!dir.exists(main_dir)) {
  dir.create(main_dir, recursive = TRUE)
}

# Ensure each worker has its own directory
worker_dirs <- file.path(main_dir, paste0("worker_", seq_len(num_cores)))
lapply(worker_dirs, dir.create, showWarnings = FALSE)

# Function applied
completed_ids <- mclapply(1:num_cores, function(core_id) {
  
  # Ensure the worker directories are accessible
  ids_for_core <- id_chunks[[core_id]]  # Get IDs assigned to this core
 
  # Apply the function to each ID
  lapply(ids_for_core, function(s) {
    
    # Call reps_fun for each unique ID
    reps_fun(
      s = s,
      wd = wd,
      worker_path = worker_dirs[[core_id]], # Ensure worker_dirs is properly indexed
      rtol_values = rtol_values,           # Ensure rtol_values is passed correctly
      atol_values = atol_values,           # Ensure atol_values is passed correctly
      params_seed_reader = params_seed_reader
    )
    
    # Return the ID after processing
    return(s)
  })
  
  # Return all IDs processed by this core
  #return(processed_ids)
}, mc.cores = num_cores)


#--------------------------------------Update RDS file for control-------------------------#

# RDS file path
rds_path <- file.path(main_dir, "completed_ids.rds")

# Unlist the completed results
x <- unlist(completed_ids)

# Read existing data if the RDS file exists, otherwise initialize an empty vector
data <- if (file.exists(rds_path)) { readRDS(rds_path) } else { c() }

# Combine old and new data
save_data <- c(data, x)

# Save the updated data back to the RDS file
saveRDS(save_data, rds_path)





