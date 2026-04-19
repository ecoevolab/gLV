#  Parameter Exploratory analysis
#
# Description: 
# This script is for recreating the parameter space we can simulate.

#------------------------------------------
# Generate directories
#------------------------------------------
tictoc::tic("Section 0: Total running time")

tictoc::tic("Section 1: Time for Parameter Generation")
#' Indicate directories paths
result_dir <- "/mnt/data/sur/users/mrivera/Data/PEA"    # directory to save      
if (!dir.exists(result_dir)) {
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE) # Create directory if it doesn't exist
}        
params_path <- file.path(result_dir, "simulation_params.tsv")     # Parameters-TSV

cat(">> The experiment path is:", result_dir,"\n", sep="")

#------------------------------------------
# Section: Generate data with different interactions
#-------------------------------------------
library(dplyr)
generate_params <- function (){
  p_neg = seq(0, 1, by = 0.1)       # negative-interactions
  p_noint = seq(0, .9, by = 0.1)     # null-interactions                    
  n_species = rep(c(10,30,50,100), each = 10) # number of species               
  dt <- data.table::CJ(n_species, p_neg, p_noint ) 
  dt[, p_pos := 1 - p_neg]

  # Add Columns with data table operator `:=`
  n_total <- nrow(dt)
  all_seeds <- sample.int(3e6L, 3L *n_total, replace = FALSE)
  
  dt[, `:=`(
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
tictoc::toc() # For section 1

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

#------------------------------------------
# Section: Run simulations
#-------------------------------------------
library(arrow)
library(data.table)
library(parallel)

tictoc::tic("Section 2: Run simulations and extinctions using the parallel package")
wrapper <- function(index, df_params) {
  #-----------------------------
  # Section: Generate parameters and run simulation
  row = df_params[index, ] 
  sim_id = row$id
  # Generate parameters
  params <- gen_training_params(row)            # Generate-parameters
  output <- solve_gLV(times = 1000, params)     # Run-simulation
  na_count <- any(is.na(output))                     # simulation-NAs
  #-----------------------------
  cat(">> Simulation ", sim_id, " completed.\n")  
  return(list(id = sim_id, na_count = na_count))
}

# Assign cores
ncore = max(1, detectCores() - 1, na.rm = TRUE)
cat('>> The number of cores to use are: ', ncore, '\n')

# Parallelize
results_summary <- mclapply(
  seq_len(nrow(params_df)),
  # seq_len(200),
  wrapper,
  df_params       = params_df,
  mc.cores        = ncore
)

# Generate information file
result_df <- data.table::rbindlist(results_summary, use.names = TRUE)
info_path <- file.path(result_dir, "params_space_smry.feather")  # Information-TSV
arrow::write_feather(x = as.data.frame(result_df), sink = info_path)

#------------------------------------------
# Section: Generate summary per parameter combination
#-------------------------------------------
library(dplyr)

summary_df <- params_df %>% 
select( n_species, p_neg, p_noint, p_pos, id) %>%   # Select columns
filter(id %in% result_df$id) %>%                    # Filter by ids in results
left_join(result_df, by = "id") %>%                 # Joint na counts
group_by(n_species, p_neg, p_noint) %>%             # group by parameters
summarise(na_count_mean= mean(na_count == TRUE, na.rm = TRUE), .groups = "drop") 

tictoc::toc() 

#------------------------------------------
# Section: Generate heatmap of the parameter space
#-------------------------------------------
library(ggplot2)
library(plotly)
library(tidyr)

# Custom theme
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

# Create plot
p <- ggplot(summary_df, aes(x = p_neg, y = p_noint, fill = na_count_mean)) +
    facet_wrap(~ n_species, labeller = label_both) +
    # Add tiles with heatmap colors
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradientn(colours = c("#fee8c8", "#e34a33"), labels = scales::percent, values = c(0,1)) +
    # Scale axes 
    scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
    # theme_minimal(base_size = 11, base_family = "serif") +
    # Custom theme adjustments
    custom_theme() +
    # Add labels
    labs(
        title    = "Proportion of failed simulations by parameter combination",
        subtitle = "Arrows indicate the increase in null(gray) and negative (red) interactions",
        x        = "Prop. negative interactions (off-diagonal and not null)",
        y        = "Prop. null interactions (off-diagonal)",
        fill     = "Failed\nsimulations"
    ) +
    # Highlighted row
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.95, ymax = 1.05, fill = alpha("#65e7eb77", 0.2), color = "white", linewidth = 0.6) +
    # Negative interaction arrow
    annotate("segment", x = 0, xend = 0.9, y = 0, yend = 0, colour = "#c0392b", linewidth = 0.9, arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + 
    # Null interaction arrow
    annotate("segment", x = 0, xend = 0, y = 0, yend = 0.9, colour = "#808080", linewidth = 0.9, arrow = arrow(length = unit(0.2, "cm"), type = "closed"))
    

# Save plot in best quality
ggsave(filename = file.path(result_dir, "params_space_heatmap.png"), plot = p, width = 10, height = 8, dpi = 300)
