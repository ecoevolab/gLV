# =============================================================================
# Extinctions impact analysis: full community vs sub-community
# =============================================================================
#
# Interaction matrix method: Column boosting with fixed diagonal and column K
#
# =============================================================================

#-------------------------------------
# Section: Load data
#-------------------------------------
out_path = '/mnt/data/sur/users/mrivera/Data/clean_controls/ImpactAnalysis_9ee93e'
full_dir  <- paste0(out_path, "/Full_ExtSummaries")
sub_dir <- paste0(out_path, "/Sub_ExtSummaries")
out_dir  <- paste0(out_path, "/RawOutputs")
sim_params = data.table::fread(paste0(out_path, "/simulation_params.tsv"), select = c("id", "key", "p_noint"))

library(dplyr)
library(parallel)

process_row <- function(i, ext_dir) {
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

results <- mclapply(1:nrow(sim_params), ext_dir = full_dir, process_row, mc.cores = detectCores() - 1)
full_summary <- do.call(rbind, results) 
results <- mclapply(1:nrow(sim_params), ext_dir = sub_dir, process_row, mc.cores = detectCores() - 1)
sub_summary <- do.call(rbind, results) 


#-------------------------------------
# Section: Metric distributions
#-------------------------------------
library(ggplot2)
library(tidyr)
# Define common theme and colors
theme_custom <- theme_minimal(base_size = 12) +
    theme(
        legend.position = "none",                                                    
        plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        axis.title = element_text(size = 10),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
        plot.subtitle = element_text(face = "italic", size = 9, hjust = 0.5),
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "grey92", color = NA)             
    )

# Only plot these variables
variables_to_plot <- c("keystoneness", "dissimilarity", "prop_extinctions")
variable_labels <- c(keystoneness    = "Keystoneness", dissimilarity   = "Dissimilarity", prop_extinctions = "Prop. Extinctions")

make_boxplot <- function(data, title, subtitle) {
    data %>%
        pivot_longer(cols = all_of(variables_to_plot), names_to = "variable", values_to = "value") %>%
        mutate(variable = factor(variable, levels = variables_to_plot)) %>%
        ggplot(aes(x = is_keystone, y = value, fill = is_keystone)) +
        geom_jitter(aes(color = is_keystone), width = 0.15, size = 0.3, alpha = 0.15, show.legend = FALSE) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA) +
        scale_y_continuous(trans = "log1p") +
        scale_x_discrete(labels = c("Non-keystone", "Keystone")) +
        scale_fill_manual(values = c("FALSE" = "#4575b4", "TRUE" = "#d73027")) +
        scale_color_manual(values = c("FALSE" = "#4575b4", "TRUE" = "#d73027")) +
        facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = variable_labels)) +
        labs(title = title, subtitle = subtitle, x = "Species type", y = "Value (log1p scale)") +
        theme_custom
}

# Plot 1: Metric distributions
method <- "Column boosting method with sub-community extinction impact"
make_boxplot(
    data = sub_summary %>% select(is_keystone, all_of(variables_to_plot)),
    title = "Distribution of metrics by keystone status",
    subtitle = method
) %>% ggsave(file.path(out_path, 'sub_metrics_distribution.png'), plot = ., width = 12, height = 5)

# Plot 3: Metrics distribution removing extinct species
method <- "Column boosting method with full-community extinction impact"
make_boxplot(
    data = full_summary %>% select(is_keystone, rel_pop_initial, all_of(variables_to_plot)) %>% filter(rel_pop_initial > 1e-06),
    title = "Distribution of metrics by keystone status for surviving species",
    subtitle = method
) %>% ggsave(file.path(out_path, 'full_metrics_distribution.png'), plot = ., width = 12, height = 5)