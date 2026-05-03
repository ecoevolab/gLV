# Plotter for solver tolerances analysis
#
# Description:
# Communities were simulated at different solver tolerances to test wether solver parameter election would affect simulation outcome.
#
# May 02, 26


# Load data
dir = '/mnt/data/sur/users/mrivera/Data/PEA/tol_analysis'
data = arrow::read_feather(file.path(dir, 'tolerance_summary.feather'))

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

p <- ggplot(data, aes(x = atol, y = rtol, fill = mean_na)) +
    # Add tiles with heatmap colors
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradientn(colours = c("#fee8c8", "#e34a33"), values = c(0,1), limits  = c(0, 1), labels  = scales::percent) +
    # Custom theme adjustments
    custom_theme() +
    # Add labels
    labs(
    title    = "Proporción de simulaciones fallidas por tolerancias del solver",
    subtitle = "Solver utilizado: ode45",
    x        = "Tolerancia absoluta (atol)",
    y        = "Tolerancia relativa (rtol)",
    fill     = "Simulaciones\nfallidas"
    )

# Save plot in best quality
save_dir = '/mnt/data/sur/users/mrivera/thesis_plots'
ggsave(filename = file.path(save_dir, "pea_solver_tols.png"), plot = p, width = 10, height = 8, dpi = 300)