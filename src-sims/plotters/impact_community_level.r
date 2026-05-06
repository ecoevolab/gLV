"""
Code for disproportion analysis.

We want to compare how the removal survival nodes impact the community when the community is composed of all nodes or only survival nodes.
"""

# Import libraries
library(dplyr)
library(parallel)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggridges)

#-----------------------------------------
# Section: Generate data
#-----------------------------------------
# Load which communities have extinctions.
dir = '/mnt/data/sur/users/mrivera/Data/clean_controls/Kboost_proof_0c7b0f'
ids_extinct = arrow::read_feather(paste0(dir,'/simulation_summary.feather')) %>% pull(id)
sim_params = data.table::fread(file.path(dir, 'simulation_params.tsv')) %>% filter(id %in% ids_extinct)

# Load data
full_impact_path = file.path(dir, 'Full_ExtSummaries')
sub_impact_path = file.path(dir, 'Sub_ExtSummaries')



generate_df <- function(i) {
    #-----------------------
    # Simulation 
    row <- sim_params[i, ]
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
results <- mclapply(1:length(ids_extinct), generate_df, mc.cores = detectCores() - 1)
# Generate full tables
full_master <- map_dfr(results, "full")
sub_master  <- map_dfr(results, "sub")

#-----------------------------------------
# Section: Plot impact
#-----------------------------------------
save_dir = '/mnt/data/sur/users/mrivera/thesis_plots'

plot_impact <- function(data, title, filename, width = 12, height = 5) {
  p <- data %>%
    select(-rel_pop_initial) %>%
    pivot_longer(-keystone, names_to = "variable", values_to = "value") %>%
    ggplot(aes(x = keystone, y = value, fill = keystone)) +
    geom_boxplot(outliers = FALSE) +
    facet_wrap(~variable, scales = "free_y") +
    labs(title = title) +
    theme_bw()
  ggsave(filename, plot = p, width = width, height = height)
  invisible(p)
}

plot_impact(data = full_master, title = "Impacto en la comunidad completa", filename = file.path(save_dir, "ctrl_full_impact.png") )
plot_impact(data = sub_master,  title = "Impacto en la comunidad de sobrevivientes",  filename = file.path(save_dir, "ctrl_sub_impact.png") )
