


#-------------------------------------
# Section: Load data
#-------------------------------------
dir = '/mnt/data/sur/users/mrivera/Data/Controls/NegCtrl-V1'
ext_dir = paste0(dir, "/ExtSummaries")
out_dir = paste0(dir, "/RawOutputs")
sim_params = data.table::fread(paste0(dir, "/simulation-params.tsv"), select = c("id", "key", "p_noint"))

library(dplyr)
library(parallel)

process_row <- function(i) {
  row <- sim_params[i, ]
  id <- row$id
  key <- as.numeric(row$key)
  p_noint <- as.numeric(row$p_noint)
  
  # Generate path
  ext_file_path <- file.path(ext_dir, paste0("ExtSummary_", id, ".feather"))
  # Generate summary
  arrow::read_feather(ext_file_path) %>%
    rename(dissimilarity = dissimilarity_bc) %>%
    mutate(
      is_keystone = row_number() == key,
      p_noint = p_noint,
      specie = row_number()
    ) %>%
    select(-pop_initial,-p_noint, -time_stability)
}

t = paste0('ExtSummary_', sim_params$id, '.feather')
t = t[t !%in% list.files(ext_dir)]
x = process_row(1)
colnames(x)
results <- mclapply(1:nrow(sim_params), process_row, mc.cores = detectCores() - 1)
full_summary <- do.call(rbind, results)

#-------------------------------------
# Section: Metric distributions
#-------------------------------------
library(ggplot2)

# Define common theme and colors
theme_custom <- theme_minimal(base_size = 12) +
    theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        axis.title = element_text(size = 10),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
        strip.text = element_text(face = "bold", size = 10) # facet title style
    )

# Only plot these variables
variables_to_plot <- c("keystoneness", "dissimilarity", "prop_extinctions")

# Generate plot
plt <- full_summary %>%
    select(is_keystone, all_of(variables_to_plot)) %>%
    pivot_longer(cols = all_of(variables_to_plot), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, levels = variables_to_plot)) %>%
    ggplot(aes(x = is_keystone, y = value, fill = is_keystone)) +
    geom_boxplot(
        alpha = 0.6, outlier.size = 0.5, outlier.alpha = 0.2
    ) +
    scale_y_continuous(trans = "log1p") +
    scale_x_discrete(labels = c("Non-keystone", "Keystone")) +
    scale_fill_manual(values = c("FALSE" = "#4575b4", "TRUE" = "#d73027"), labels = c("Non-keystone", "Keystone")) +
    facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = variable_labels)) +
    labs(
        title = "Distribution of metrics by keystone status",
        x = "Species type",
        y = "Value (log1p scale)",
        fill = "Species type"
    ) +
    theme_custom

ggsave(file.path(dir, 'all_metrics_distribution.png'), plot = plt, width = 12, height = 5)

#-------------------------------------
# Section: Metrics vs Relative Abundance for keystone species
#-------------------------------------
library(hexbin)
colors = scale_color_manual(values = c("Non-keystone" = "#4575b4", "Keystone" = "#d73027"))

plt <- full_summary %>%
    filter(is_keystone) %>%
    select(rel_pop_initial, keystoneness, dissimilarity, prop_extinctions) %>%
    pivot_longer(
        cols = c(keystoneness, dissimilarity, prop_extinctions),
        names_to = "variable",
        values_to = "value"
    ) %>%
    mutate(variable = factor(variable, levels = variables_to_plot)) %>%
    ggplot(aes(x = rel_pop_initial, y = value)) +
    geom_point(alpha = 0.4, size = 0.8, color = "#4575b4") +
    scale_y_continuous(trans = "log1p") +
    scale_x_continuous(trans = "log1p") +
    facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = variable_labels)) +
    labs(
        title = "Keystone species: relative abundance vs. community metrics",
        x = "Relative abundance (log1p scale)",
        y = "Metric value (log1p scale)"
    ) +
    theme_custom +
    theme(legend.position = "none")

ggsave(file.path(dir, 'keystone_relative_metrics.png'), plot = plt, width = 12, height = 5)

#-------------------------------------
# Section: Metrics distribution removing extinct species
#-------------------------------------

# Generate plot
plt <- full_summary %>%
    select(is_keystone, rel_pop_initial, keystoneness, dissimilarity, prop_extinctions) %>%
    filter(rel_pop_initial > 1e-06) %>%  # Remove extinct species
    pivot_longer(cols = all_of(variables_to_plot), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, levels = variables_to_plot)) %>%
    ggplot(aes(x = is_keystone, y = value, fill = is_keystone)) +
    geom_boxplot(
        alpha = 0.6, outlier.size = 0.5, outlier.alpha = 0.2
    ) +
    scale_y_continuous(trans = "log1p") +
    scale_x_discrete(labels = c("Non-keystone", "Keystone")) +
    scale_fill_manual(values = c("FALSE" = "#4575b4", "TRUE" = "#d73027"), labels = c("Non-keystone", "Keystone")) +
    facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = variable_labels)) +
    labs(
        title = "Distribution of metrics by keystone status for surviving species",
        x = "Species type",
        y = "Value (log1p scale)",
        fill = "Species type"
    ) +
    theme_custom

ggsave(file.path(dir, 'alive_metrics_distribution.png'), plot = plt, width = 12, height = 5)