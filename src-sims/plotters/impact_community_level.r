"""
Code for disproportion analysis.

We want to compare how the removal survival nodes impact the community when the community is composed of all nodes or only survival nodes.
"""


#-----------------------------------------
# Section: Generate data
#-----------------------------------------
# Load which communities have extinctions.
dir = '/mnt/data/sur/users/mrivera/clean_controls/91074c4e25b4'
ids_extinct = arrow::read_feather(paste0(dir,'/simulation_summary.feather')) %>% filter(ext_performed) %>% pull(id)
parameters = data.table::fread(file.path(dir, 'simulation_params.tsv')) %>% filter(id %in% ids_extinct)

# Load data
full_impact_path = file.path(dir, 'Full_ExtSummaries')
sub_impact_path = file.path(dir, 'Sub_ExtSummaries')

library(dplyr)

generate_df <- function(i) {
    #-----------------------
    # Simulation 
    row <- parameters[i, ]
    id <- row$id
    key <- as.numeric(row$key)
    p_noint <- as.numeric(row$p_noint)
    #-----------------------
    # Full-community impact
    cols_to_drop = c('specie','n_extinctions', 'time_stability', 'pop_initial')
    full_path <- paste0(full_impact_path, paste0("/ExtSummary_", id, ".feather"))
    full_df <- arrow::read_feather(full_path) %>%  mutate(keystone = if_else(specie == key, TRUE, FALSE)) %>% select(-all_of(cols_to_drop))
    # Sub-community impact
    sub_path <- paste0(sub_impact_path, paste0("/ExtSummary_", id, ".feather"))
    sub_df <- arrow::read_feather(sub_path) %>%  mutate(keystone = if_else(specie == key, TRUE, FALSE)) %>% select(-all_of(cols_to_drop))
    #------------------------
    list(full = full_df, sub = sub_df)
}

generate_df(1)

# Parallelize
library(parallel)
results <- mclapply(1:length(ids_extinct), generate_df, mc.cores = detectCores() - 1)

# Generate full tables
library(purrr)
library(tidyr)
full_master <- map_dfr(results, "full")
sub_master  <- map_dfr(results, "sub")

#-----------------------------------------
# Section: Plot impact
#-----------------------------------------
plots_dir = '/mnt/data/sur/users/mrivera/Plots'
experiment_plots = file.path(plots_dir, basename(dir))
if(!(dir.exists(experiment_plots))) {dir.create(experiment_plots)}

# Generate plot
library(ggplot2)
library(ggridges)

# Full impact plot
p1 <- full_master %>% select(-rel_pop_initial) %>%
    pivot_longer(-keystone, names_to = "variable", values_to = "value") %>%
    ggplot(aes(x = keystone, y = value, fill = keystone)) +
    geom_boxplot(outliers = FALSE) +
    facet_wrap(~variable, scales = "free_y") +
    labs(title = "Full community") +
    theme_bw()
# Save plot
full_plot_path = file.path(experiment_plots, "full_impact_distr.png")
ggsave(full_plot_path, plot = p1, width = 12, height = 5)
# Sub-community plot
p2 <- sub_master %>% select(-rel_pop_initial) %>%
    pivot_longer(-keystone, names_to = "variable", values_to = "value") %>%
    ggplot(aes(x = keystone, y = value, fill = keystone)) +
    geom_boxplot(outliers = FALSE) +
    facet_wrap(~variable, scales = "free_y") +
    labs(title = "Sub community") +
    theme_bw()
# Save plot
sub_plot_path = file.path(experiment_plots, "sub_impact_distr.png")
ggsave(sub_plot_path, plot = p2, width = 12, height = 5)
