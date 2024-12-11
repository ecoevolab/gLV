

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
reps_ids <- setdiff(table$ID, controls$ID)
control_ids <- controls$ID

#------------------------------------Regenerate simulations----------------------------#
ode_function <- function (times, params, atol, rtol) {
  
  requireNamespace("deSolve")
  
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

rtol_value <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1) # Relative tolerance
atol_value <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1) # Absolute tolerance

for (s in reps_ids) {
  
  # Read parameters
  params <- params_seed_reader(uniqueID = s, wd)
  
  # Loop in reverse order
  for (r in rev(rtol_value)) {
    for (a in rev(atol_value)) {
      # Run simulation
      new_out <- ode_function(times = 700, params, a, r)
      
      # Define paths
      tmp_path <- file.path(wd, "Tolerances", paste0("rt_", format(r, scientific = TRUE, digits = 3)))
      tmp2_path <- file.path(tmp_path, paste0("at_", format(a, scientific = TRUE, digits = 3)), "Cases")
      
      # Create directory if it does not exist
      if (!dir.exists(tmp2_path)) {
        dir.create(tmp2_path, recursive = TRUE)
      }
      
      # Save simulation
      save_path <- file.path(tmp2_path, paste0("R_", s, ".tsv"))
      utils::write.table(new_out, file = save_path, sep = "\t", row.names = FALSE, col.names = TRUE)
      
      cat("Case ", s, " tested with atol of", a, " and rtol of", r , "\n")
    }
  }
}

#-------------------------------------------Controls----------------------------------------#


sam_controls <- sample(control_ids, length(control_ids)/4)

for (s in sam_controls) {
  
  # Read parameters
  params <- params_seed_reader(uniqueID = s, wd)
  
  # Loop in reverse order
  for (r in rev(rtol_value)) {
    for (a in rev(atol_value)) {
      # Run simulation
      new_out <- ode_function(times = 700, params, a, r)
      
      # Define paths
      tmp_path <- file.path(wd, "Tolerances", paste0("rt_", format(r, scientific = TRUE, digits = 3)))
      tmp2_path <- file.path(tmp_path, paste0("at_", format(a, scientific = TRUE, digits = 3)), "Controls")
      
      # Create directory if it does not exist
      if (!dir.exists(tmp2_path)) {
        dir.create(tmp2_path, recursive = TRUE)
      }
      
      # Save simulation
      save_path <- file.path(tmp2_path, paste0("R_", s, ".tsv"))
      utils::write.table(new_out, file = save_path, sep = "\t", row.names = FALSE, col.names = TRUE)
      
      cat("Control ", s, " tested with atol of", a, " and rtol of", r , "\n")
    }
  }
}


