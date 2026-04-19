#  Parameter Exploratory analysis for tolerances
#
# Description: 
# This script is for testing the effect of different tolerances on the solver.


#------------------------------------------
# Generate directories
#------------------------------------------
result_dir <- "/mnt/data/sur/users/mrivera/Data/PEA/tol_analysis"    # directory to save      
if (!dir.exists(result_dir)) {
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE) # Create directory if it doesn't exist
}        
params_path <- file.path(result_dir, "tolerances_params.tsv")     # Parameters-TSV

cat(">> The experiment path is:", result_dir,"\n", sep="")

#------------------------------------------
# Section: Generate data with different interactions
#-------------------------------------------
library(dplyr)
generate_params <- function (){
  p_noint = seq(0, .9, by = 0.1)    # null-interactions                    
  n_species = rep(30, each = 5) # number of species               
  dt <- data.table::CJ(n_species, p_noint ) 

  # Add Columns with data table operator `:=`
  n_total <- nrow(dt)
  all_seeds <- sample.int(3e6L, 3L *n_total, replace = FALSE)
  
  dt[, `:=`(
    p_neg = 1,
    id = ids::random_id(n = n_total, bytes = 3),
    x0_seed = all_seeds[1:n_total],
    mu_seed = all_seeds[(n_total+1):(n_total*2)],
    A_seed = all_seeds[(2 * n_total + 1):(3 * n_total)]
  )]
  
  return(dt)
}

params_df <- generate_params()
# Verify if ids are unique and in case they are, save the parameters.
while (nrow(params_df) != length(unique(params_df$id))) {
    params_df <- generate_params() # Repeat function
}

data.table::fwrite(x = params_df, file = params_path, sep = "\t", quote = FALSE, row.names = FALSE) # Save parameters
message("\nParameteres generated and saved at path:\n", params_path, "\n")
cat(">> The number of simulations are:", nrow(params_df),"\n", sep=" ")


#------------------------------------------
# Section: Declare functions
#-------------------------------------------
# Source functions to:
codes = list.files("/mnt/data/sur/users/mrivera/gLV/src-sims/FUN", full.names=TRUE)

invisible(lapply(codes, function(file) {
  cat(">> Sourcing function:", file, "\n")
  capture.output(source(file))
  # No need for return() here
}))


solve_gLV <- function(times, params, atol, rtol) {
  
  # Define the gLV model
  glv_model <- function(t, x0, params) {
    x0[x0 < 1e-6] <- 0 # Ignore the effect of species with population below a threshold
    
    # dx/dt = X *(r + A * X) + D
    dx <- (x0 * (params$mu + params$M %*% x0)) + 1e-6
    list(dx)
  }
  
  time_seq <- seq_len(times)  # Times to simulate
  
  # Try solving the system with a timeout
  results <- tryCatch(
    R.utils::withTimeout(
      deSolve::ode(
        y = params$x0,
        times = time_seq,
        func = glv_model,
        parms = params,
        method = "ode45",
        rtol = rtol,
        atol = atol
      ),
      timeout = 600
    ),
    error = function(e) {
    message(">> Error details: ", e$message)
    message(">> Error class: ", class(e))
    print(str(params))  # Check params structure
    return(NULL)
  }
  )
  
  # Process results if valid
  if (!is.null(results) && ncol(results) > 1) {
    tmp <- results[, -1] |>  # remove time column
      t() |>                 # transpose
      as.data.frame() |>     # convert to df
      # Populations that went extinct (no effect)
      dplyr::mutate(across(everything(), ~ replace(., . < 1e-8, 0))) 
    
    return(tmp)
  }
  
  # If results not valid: return NA matrix 
  matrix(NA, nrow = nrow(params$M), ncol = times)
}


ref_solver <- function(times, params) {
  
  # Define the gLV model
  glv_model <- function(t, x0, params) {
    x0[x0 < 1e-6] <- 0 # Ignore the effect of species with population below a threshold
    
    # dx/dt = X *(r + A * X) + D
    dx <- (x0 * (params$mu + params$M %*% x0)) + 1e-6
    list(dx)
  }
  
  time_seq <- seq_len(times)  # Times to simulate
  
  # Try solving the system with a timeout
  results <- tryCatch(
    R.utils::withTimeout(
      deSolve::ode(
        y = params$x0,
        times = time_seq,
        func = glv_model,
        parms = params,
        method = "ode45"
      ),
      timeout = 600
    ),
    error = function(e) {
    message(">> Error details: ", e$message)
    message(">> Error class: ", class(e))
    print(str(params))  # Check params structure
    return(NULL)
  }
  )
  
  # Process results if valid
  if (!is.null(results) && ncol(results) > 1) {
    tmp <- results[, -1] |>  # remove time column
      t() |>                 # transpose
      as.data.frame() |>     # convert to df
      # Populations that went extinct (no effect)
      dplyr::mutate(across(everything(), ~ replace(., . < 1e-8, 0))) 
    
    return(tmp)
  }
  
  # If results not valid: return NA matrix 
  matrix(NA, nrow = nrow(params$M), ncol = times)
}


#------------------------------------------
# Section: Run simulations
#-------------------------------------------
library(arrow)
library(data.table)
library(parallel)

tols_list <- 10^seq(-1, -6, by = -1) # Tolerances to test
tol_df <- expand.grid(atol = tols_list, rtol = tols_list)  # 36 rows
tol_df$name <- paste0("C", 1:nrow(tol_df)) # Add name column

# Generate combination directory
data.table::fwrite(x = tol_df, file = file.path(result_dir, "combinations_tols.tsv"), sep = "\t", quote = FALSE, row.names = FALSE) 

wrapper <- function(index, params_df) {
    #-----------------------------
    # Section: Generate parameters and run simulation
    row = params_df[index, ] 
    sim_id = row$id
    params <- gen_training_params(row)            # Generate-parameters
    # Generate reference simulation with default tolerances
    ref_out <- ref_solver(times = 1000, params = params)     
    ref_vec <- as.vector(ref_out[,ncol(ref_out)]) # Reference population
    # Test different tolerances
    final_list <- lapply(1:nrow(tol_df), function(i) {
        tol_row <- tol_df[i, ]
        output <- solve_gLV(times = 1000, params = params, atol = tol_row$atol, rtol = tol_row$rtol)     # Run-simulation
        last <- output[, ncol(output)]     # Get last time point
        na_count <- any(is.na(output))     # simulation-NAs
        return(list(id = sim_id, combination = tol_row$name, final = last, na_count = na_count))
    })
    #-----------------------------
    # Perform RMSE
    rmse_df <- do.call(rbind, lapply(final_list, function(x) { 
        data.frame(
            combination = x$combination,
            rmse        = sqrt(mean((x$final - ref_vec)^2)),
            na_count    = x$na_count
        )
    }))
    cat(">> Simulation ", sim_id, " completed.\n")  
    return(rmse_df)
}

# Assign cores
ncore = max(1, detectCores() - 1, na.rm = TRUE)
cat('>> The number of cores to use are: ', ncore, '\n')

# Parallelize
results_summary <- mclapply(
  seq_len(nrow(params_df)),
  # seq_len(200),
  wrapper,
  params_df       = params_df,
  mc.cores        = ncore
)

# Join all results
all_rmse <- bind_rows(results_summary)
#------------------------------------------
# Section: Summary examples
#-------------------------------------------
library(dplyr)

summary_rmse <- all_rmse %>%
    group_by(combination) %>%
    summarise(
        mean_rmse = mean(rmse),
        sd_rmse   = sd(rmse),
        max_rmse  = max(rmse),
        mean_na  = mean(na_count)
    ) %>%
    left_join(tol_df, by = c("combination" = "name")) %>%
    arrange(atol, rtol)
summary_rmse$atol <- factor(summary_rmse$atol)
summary_rmse$rtol <- factor(summary_rmse$rtol)
print(summary_rmse, n = 40)
arrow::write_feather(x = summary_rmse, sink = file.path(result_dir, "tolerance_summary.feather"))

#------------------------------------------
# Section: Generate heatmap of the parameter space
#-------------------------------------------
library(ggplot2)

custom_theme <- function(){
  theme(
    plot.title      = element_text(hjust = 0.5, size = 13, face = "bold", margin = margin(b = 10)),
    plot.subtitle   = element_text(hjust = 0.5, size = 9, color = "grey40", margin = margin(b = 10)),
    strip.text      = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", color = NA),
    axis.title      = element_text(size = 10),
    axis.text       = element_text(size = 8, color = "grey30"),
    legend.title    = element_text(size = 9, face = "bold"),
    legend.text     = element_text(size = 8),
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 1.5),
    plot.margin     = margin(10, 15, 10, 10)
  )
}

p <- ggplot(summary_rmse, aes(x = atol, y = rtol, fill = mean_na)) +
    # Add tiles with heatmap colors
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradientn(colours = c("#fee8c8", "#e34a33"), values = c(0,1), limits  = c(0, 1), labels  = scales::percent) +
    # Custom theme adjustments
    custom_theme() +
    # Add labels
    labs(
        title    = "Proportion of failed simulations by solver tolerances",
        subtitle = "Solver used: ode45",
        x        = "Absolute tolerance (atol)",
        y        = "Relative tolerance (rtol)",
        fill     = "Failed\nsimulations"
    ) 

# Save plot in best quality
ggsave(filename = file.path(result_dir, "tolerance_heatmap.png"), plot = p, width = 10, height = 8, dpi = 300)