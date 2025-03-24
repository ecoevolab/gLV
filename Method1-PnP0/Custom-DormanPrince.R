


# Load libraries
library(dplyr) 
library(tidyr)
require(tidyverse, lib.loc = "/mnt/atgc-d3/sur/modules/pkgs/tidyverse_mrc")
library(purrr)
library(ggplot2)
library(svglite, lib.loc = "/mnt/atgc-d3/sur/modules/pkgs/tidyverse_mrc")
library(plotly)

#======= Prepare data ======

# Define working directory
outs_dir <-  "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Exp03-D25M02/Simulate_ODE/Unified"

# Load parameters data
params_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Exp03-D25M02.tsv")
data.table::setDT(params_table)
head(params_table)

# Load Na counting table
na_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Exp03-D25M02/Nas-counting.tsv")
head(na_table)

# Filtering step
na_ids <- na_table$id[na_table$Total.NAs != 0]
df_table <-  params_table[id %in% na_ids]
head(df_table)

# Source gLV parameter generation function
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/Method1-PnP0/Forge-gLV-Parameters.R")
print(regenerate)

#======= Define gLV model for Dormand Prince ======
# Define the equation
glv_model <- function(x0, params) {
  r <- params$mu         # Growth rate vector
  A <- params$M          # Interaction matrix
  
  # Compute dx/dt for each species
  dx <- x0 * (r + A %*% x0)
  return(as.vector(dx))
}


# ==== ODE solver ====
ode.simulate <- function(times, params) {
  
  # Define the equation
  glv_model <- function(t, x0, params) {
    r <- params$mu         # Growth rate vector
    A <- params$M          # Interaction matrix
    
    # Compute dx/dt for each species
    dx <- x0 * (r + A %*% x0)
    list(dx)
  }
  
  time_seq <- seq(1, times, by = 1)  # Define the time sequence
  
  # Get solution
  results <- tryCatch(
    R.utils::withTimeout(deSolve::ode(y = params$x0, times = time_seq, func = glv_model, 
                                      parms = params,
                                      method = "ode45",
                                      rtol = 1e-06, 
                                      atol = 1e-06),
                         timeout = 600), 
    error = function(e) {
      message(">> Simulation failed... skipping")
      return(NULL) # Return NULL instead of an NA matrix
    })
  
  
  # Remove the first column (`time`) and transpose it so columns represent generations
  # Check for valid output and return transposed results
  if (!is.null(results) && ncol(results) > 1) {
    return(t(results[, -1]))
  } else {
    return(matrix(NA, nrow = nrow(params$M), ncol = times)) # Return properly shaped NA matrix
  }
}

#======= Identify failed generation ======
search_failed <- function(outs_dir, params){
  
  # Create path to directory
  na_output <- data.table::fread(file.path(outs_dir, paste0("O_", params$id, ".tsv")) )
  
  # Find at which generation NAs appear
  na_cols <- which(colSums(is.na(na_output)) > 0)
  return(na_cols[1])
}

#======= Define Dorman-Prince method ======

DP_RK <- function(glv_model, params) {
  
  iter <- 1 # Iteration counter
  epsilon <- 1e-06 # Specify error thershold
  h0 <- 0.5 # First step size
  
  # Steps vector
  # ci <- c(0, 1/5, 3/10, 4/5, 8/9, 1, 1)
  
  # Coefficients for stages
  values_aij <- c(0, 0, 0, 0, 0, 0, # K1
                  1/5, 0, 0, 0, 0, 0, # K2
                  3/40, 9/40, 0, 0, 0, 0, # K3
                  44/45, -56/15, 32/9, 0, 0, 0, # K4
                  19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, # K5
                  9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0, # K6
                  35/384, 0, 500/1113, 125/192, -2187/6784, 11/84) # K7
  
  aij <- matrix(values_aij, nrow = 7, ncol=6, byrow = TRUE)
  
  k_steps <- function(params) {
    # Calculate each step 
    x0 <- params$x0
    k1 = h0 * glv_model(x0, params)
    
    x0 <- params$x0 + (aij[2,1]*k1) 
    k2 = h0 * glv_model(x0, params)
    
    x0 <- params$x0 + (aij[3,1]*k1) + (aij[3,2]*k2)
    k3 = h0 * glv_model(x0, params)
     
    x0 <- params$x0 + (aij[4,1]*k1) + (aij[4,2]*k2) + (aij[4,3]*k3) 
    k4 = h0 * glv_model(x0, params)
     
    x0 <- params$x0 + (aij[5,1]*k1) + (aij[5,2]*k2) + (aij[5,3]*k3) + (aij[5,4]*k4)
    k5 = h0 * glv_model(x0, params)
     
    x0 <- params$x0 + (aij[6,1]*k1) + (aij[6,2]*k2) + (aij[6,3]*k3) + (aij[6,4]*k4)  + (aij[6,5]*k5)
    k6 = h0 * glv_model(x0, params)
    
    x0 <- params$x0 + (aij[7,1]*k1) + (aij[7,2]*k2) + (aij[7,3]*k3) + (aij[7,4]*k4)  
     + (aij[7,5]*k5) + (aij[7,6]*k6)
    k7 = h0 * glv_model(x0, params)
    
    # Generate list of steps
    ks <- list(k1, k2, k3, k4, k5, k6, k7)
    return(ks)
  }
  
  # Calculate steps
  ks <- k_steps(params)
  
  # Coefficients for RK4(b4)
  b4 <- c(35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0)
  weighted_ks <- mapply(function(k, b) b * k, k = ks, b= b4, SIMPLIFY = FALSE)
  y_new4 <- params$x0 + Reduce("+", weighted_ks) # Sum the weighted k vectors to get the new x
  
  # Coefficients for RK5 (b5)
  b5 <- c(5179/57600, 0, 7571/16695, 393/640, - 92097/339200, 187/2100, 1/40)
  weighted_ks <- mapply(function(k, b) b * k, k = ks, b= b5, SIMPLIFY = FALSE)
  y_new5 <- params$x0 + Reduce("+", weighted_ks) # Sum the weighted k vectors to get the new x
  
  # Update step size if failed
  # Order used for error estimation (4) + 1
  error <- abs(y_new5 - y_new4)
  
  while (any(error > epsilon)) {
    h0 <- 0.9 * h0 * (epsilon / error)^(1/5)  # Update step size
    ks <- k_steps(params)  # Calculate steps (only once)
    
    # Coefficients for RK4(b4)
    weighted_ks <- mapply(function(k, b) b * k, k = ks, b= b4, SIMPLIFY = FALSE) # Sum the weighted k vectors for RK4
    y_new4 <- params$x0 + Reduce("+", weighted_ks) # RK4 solution
    
    # Coefficients for RK5 (b5)
    weighted_ks <- mapply(function(k, b) b * k, k = ks, b= b5, SIMPLIFY = FALSE) # Sum the weighted k vectors for RK5
    y_new5 <- params$x0 + Reduce("+", weighted_ks) # RK5 solution
    
    # Update error estimation
    error <- abs(y_new5 - y_new4)
    
    iter = iter + 1 # Update iteration counter
    
  }
  
  message(paste("Simulation", params$id, "met tolerance at iteration", iter, "\n"))
  
  return(y_new5)
  
}

#======= Example step ======

# Select one column
index <- df_table[1,]

# Regenerate parameters
params <- regenerate(index)
params

# Simulate with ode solver
output1 <- ode.simulate(times = 100, params)
output2 <- map_dfc(1:100, ~ {
  params$x0 <- DP_RK(glv_model, params)  # Update x0
  as.data.frame(params$x0)  # Ensure returning a data frame
})


# Find column at wchich NA values appear
output_readed <- data.table::fread(file.path(outs_dir, paste0("O_", params$id, ".tsv")) )
na_col <- search_failed(outs_dir, params)
tmp <- as.numeric(na_col - 1)
old_x0 <- output[,1] 
new_x0 <- output[,..tmp] 
cat("The old x0 value is ", unlist(output_readed[,1]), " \nthe new value is: ", unlist(new_x0))

# Set population parameters one time before Na values appear
params$x0 <- unname(unlist(new_x0))

# Apply Dormand Prince function 100 times and bind the results as columns
results <- map_dfc(1:100, ~ {
  params$x0 <- DP_RK(glv_model, params)  # Update x0
  as.data.frame(params$x0)  # Ensure returning a data frame
})

results <- as.data.frame(results) %>%
  rename_with(~ as.character(1:ncol(results))) %>% # Assign time points 
  mutate(Species = factor(1:n())) %>%  # Assign species numbers
  pivot_longer(cols = -Species, names_to = "Time", values_to = "Frequency") %>%
  mutate(Time = as.numeric(Time))  # Convert time indices to numeric

# Plot
p1 <- ggplot(results, aes(x = Time, y = Frequency, color = Species, group = Species)) +
  geom_line() +
  geom_point() +
  labs(x = "Time", y = "Frequency", title = "Species Frequencies Over Time") +
  theme_minimal()

interactive_plot1 <- plot_ly(results, x = ~Time, y = ~Frequency, color = ~Species, 
                             type = "scatter", mode = "lines+markers") %>%
  layout(title = "Species Frequencies Over Time",
         xaxis = list(title = "Time"),
         yaxis = list(title = "Frequency"),
         hovermode = "closest")

# Save plots
htmlwidgets::saveWidget(interactive_plot1, "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/Exp03-D25M02-DPRK.html")
ggsave("/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/Exp03-D25M02-DPRK.svg", width = 6, height = 4)

#======= Work in progress ======

# Regenerate parameters and apply Dorman-Prince method
tmp <- df_table[1,]
apply(df_table, 1, function(x){
  
  # Regenerate parameters
  params <- regenerate(x)
  cat("Parameters generated for simulation: ", params$id, "\n")
  
  # Apply Dorman Prince method
})