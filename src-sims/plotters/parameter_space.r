
#------------------------------------------
# Section: Generate summary per parameter combination
#-------------------------------------------
library(dplyr)

save_dir <- "/mnt/data/sur/users/mrivera/thesis_plots"        # directory to save   
params_df <- data.table::fread('/mnt/data/sur/users/mrivera/Data/PEA/simulation_params.tsv')
result_df <- arrow::read_feather('/mnt/data/sur/users/mrivera/Data/PEA/params_space_smry.feather')

# Load data
summary_df <- params_df %>% 
select( n_species, p_neg, p_noint, p_pos, id) %>%   # Select columns
filter(id %in% result_df$id) %>%                    # Filter by ids in results
left_join(result_df, by = "id") %>%                 # Joint na counts
group_by(n_species, p_neg, p_noint) %>%             # group by parameters
summarise(na_count_mean= mean(na_count == TRUE, na.rm = TRUE), .groups = "drop") 

# Add missing combinations with NA values for the heatmap
missing_rows <- expand.grid(
    p_noint      = 1,
    p_neg        = unique(summary_df$p_neg),
    n_species    = unique(summary_df$n_species)
)
missing_rows$na_count_mean <- NA  # will get the special color

# Add to summary_df
summary_df <- rbind(summary_df, missing_rows)

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
    # panel.border    = element_rect(color = "black", fill = NA, linewidth = 1.5),
    plot.margin     = margin(10, 15, 10, 10)
  )
}

p <- ggplot(summary_df, aes(x = p_neg, y = p_noint, fill = na_count_mean)) +
    facet_wrap(~ n_species, labeller = as_labeller(c("10" = "n = 10", "30" = "n = 30", "50" = "n = 50", "100" = "n = 100"))) +
    # Add tiles with heatmap colors
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradientn(colours = c("#fee8c8", "#e34a33"), labels = scales::percent, values = c(0,1), na.value = "#2c3e50") +
    # Scale axes 
    scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
    # theme_minimal(base_size = 11, base_family = "serif") +
    # Custom theme adjustments
    custom_theme() +
    # Add labels
    labs(
        title    = "Proporción de simulaciones fallidas por combinación de parámetros",
        subtitle = "Las flechas indican el aumento en interacciones nulas (gris) y negativas (rojo)",
        x        = "Prop. interacciones negativas (fuera de la diagonal y no nulas)",
        y        = "Prop. interacciones nulas (fuera de la diagonal)",
        fill     = "Simulaciones\nfallidas"
    ) +
    # Not simulated row
    annotate("text", x = 0.5, y = 1.0, label = "Excluido", color = "white", size = 3.5, fontface = "italic") +
    # Selected column of parameters
    geom_rect(xmin = 0.95, xmax = 1.05, ymin = -0.05, ymax = 0.95, fill = NA, color = "#2ecc71", linewidth = 0.4) +
    annotate("text", x = 1, y = 0.5, label = "Elegido", color = "#2ecc71", size = 3.5, fontface = "italic", angle = 90) +
    # geom_rect(xmin = 0.95, xmax = 1.05, ymin = -Inf, ymax = Inf,  color = "#a8d5b5", linewidth = 0.6) +
    # Negative interaction arrow
    annotate("segment", x = 0, xend = 1, y = 0, yend = 0, colour = "#c0392b", linewidth = 0.9, arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + 
    # Null interaction arrow
    annotate("segment", x = 0, xend = 0, y = 0, yend = 0.9, colour = "#808080", linewidth = 0.9, arrow = arrow(length = unit(0.2, "cm"), type = "closed")) # +
  # Selected examples
    # geom_point(aes(x = 0, y = 0), shape = 21, size = 5, fill = NA, color = "#6c5ce7", stroke = 0.6) +
    # geom_point(aes(x = 1, y = 0.9), shape = 21, size = 5, fill = NA, color = "#6c5ce7", stroke = 0.6)



# Save plot in best quality
ggsave(filename = file.path(save_dir, "params_space_heatmap.png"), plot = p, width = 10, height = 8, dpi = 300)
