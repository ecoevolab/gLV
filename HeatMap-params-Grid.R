
#------------Load master table#------------#
data <- read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/D13M02Y25/NAcount-D13M02Y25.tsv", sep = "\t", header = TRUE)

params_table <- read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/D13M02Y25.tsv", sep = "\t", header = TRUE)

# Cargar paquetes
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)

# Remover columna de IDs y modificar valores > 0
df_mutate <- data %>%
  mutate(across(-1, ~ ifelse(. > 0, 1, .))) %>% # Do 0 or 1 if simulations had NA values
  # Join parameters by id
  left_join(params_table %>% select(id, p_neg, p_noint, n_species), by = c("TSV_ID" = "id")) %>% 
  mutate(across(c(p_neg, p_noint), ~ round(.x, 2))) # Make double decimals the parameters

# Sumar y contar valores por grupo
df_summed <- df_mutate %>%
  group_by(p_neg, p_noint) %>% # Group by parameters
  summarise(
    NA_Simulations = sum(NA_Counts, na.rm = TRUE),  # Count the sum of NA simulations
    Total_Sims = n(),  # Count total simulations
    .groups = "drop"
  ) %>%
  mutate(NA_SimsProps = NA_Simulations / Total_Sims) # Compute proportions

# Create ggplot object with interactive text
p <- ggplot(df_summed, aes(x = p_neg, y = p_noint, fill = NA_SimsProps, 
                           text = paste("Total Sims:", Total_Sims,
                                        "<br>p_neg:", p_neg,
                                        "<br>p_noint:", p_noint,
                                        "<br>NA Sims:", NA_Simulations,
                                        "<br>Proportions:", NA_SimsProps
                                        ))) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Proportions HeatMap of NA Counts by Parameters ODE method", 
       x = "Negative Probability", 
       y = "No Interaction Probability", 
       fill = "NAs") +
  theme(plot.title = element_text(hjust = 0.5))  # Center title

# Convert to interactive plotly with hover text
p_interactive <- ggplotly(p, tooltip = "text") %>%
  layout(hoverlabel = list(bgcolor = "white"))

# Save as interactive HTML if needed
htmlwidgets::saveWidget(p_interactive, "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/Params/D13M02Y25-C0CN.html", selfcontained = FALSE)
