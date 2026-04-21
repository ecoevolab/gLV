# =============================================================================
# Control simulations — cascade keystoneness 
# =============================================================================
#
# Interaction matrix setup:
#   - Diagonal:          fixed at -0.5
#   - Cascade interactions cell: -U(0,1) * 10
#   - Other off-diagonal and not cascade specific interactions: sparse mix of zeros and negatives from U(0,1)
#
# =============================================================================
#-------------------------------------
# Section: Load data
#-------------------------------------
out_path <- '/mnt/data/sur/users/mrivera/Data/Controls/Boosted_keystone'
ext_dir  <- paste0(out_path, "/ExtSummaries")
out_dir  <- paste0(out_path, "/RawOutputs")
sim_params = data.table::fread(paste0(out_path, "/simulation-params.tsv"), select = c("id", "key", "p_noint"))

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

results <- mclapply(1:nrow(sim_params), process_row, mc.cores = detectCores() - 1)
full_summary <- do.call(rbind, results)

#-------------------------------------
# Section: Metric distributions
#-------------------------------------
library(ggplot2)

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
make_boxplot(
    data = full_summary %>% select(is_keystone, all_of(variables_to_plot)),
    title = "Distribution of metrics by keystone status",
    subtitle = "Column boosting method"
) %>% ggsave(file.path(out_path, 'all_metrics_distribution.png'), plot = ., width = 12, height = 5)
print(paste0('plot saved at: ', file.path(out_path, 'all_metrics_distribution.png')))

# Plot 3: Metrics distribution removing extinct species
make_boxplot(
    data = full_summary %>% select(is_keystone, rel_pop_initial, all_of(variables_to_plot)) %>% filter(rel_pop_initial > 1e-06),
    title = "Distribution of metrics by keystone status for surviving species",
    subtitle = "Column boosting method"
) %>% ggsave(file.path(out_path, 'alive_metrics_distribution.png'), plot = ., width = 12, height = 5)
print(paste0('plot saved at: ', file.path(out_path, 'alive_metrics_distribution.png')))
#-------------------------------------
# Section: Metrics vs Relative Abundance for keystone species
#-------------------------------------
library(ggtext)  

threshold <- 1e-06

# Generate plot data
counts <- plot_data %>%
    group_by(variable, low_abundance) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(color = ifelse(low_abundance, "#d73027", "#4575b4")) %>%
    group_by(variable) %>%
    summarise(
        label = paste0(
            "<span style='color:", color, "'>&#9679;</span> ", n,   # colored circle with count
            collapse = "<br>"
        ),
        .groups = "drop"
    ) %>%
    mutate(x_pos = Inf, y_pos = Inf)

# Plotting
plt <- plot_data %>%
    ggplot(aes(x = rel_pop_initial, y = value, color = low_abundance)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_richtext(                                                                  
        data = counts,
        aes(x = x_pos, y = y_pos, label = label),
        hjust = 1.1, vjust = 1.5, size = 3,
        fill = "white", color = "grey30",
        label.r = unit(0.25, "lines"),                                              # rounded corners
        label.padding = unit(c(0.3, 0.4, 0.3, 0.4), "lines"),
        inherit.aes = FALSE
    ) +
    scale_color_manual(
        values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
        labels = c("TRUE" = paste0("< ", threshold), "FALSE" = paste0("≥ ", threshold)),
        name = "Relative abundance"
    ) +
    scale_y_continuous(trans = "log1p") +
    scale_x_continuous(trans = "log1p") +
    facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = variable_labels)) +
    labs(
        title = "Keystone species: relative abundance vs. community metrics",
        subtitle = "Column boosting method",
        x = "Relative abundance (log1p scale)",
        y = "Metric value (log1p scale)"
    ) +
    theme_custom +
    theme(legend.position = "right")

ggsave(file.path(out_path, 'keystone_relative_metrics.png'), plot = plt, width = 12, height = 5)
print(paste0('plot saved at: ', file.path(out_path, 'keystone_relative_metrics.png')))