
#------------Load master table#------------#
data <- read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Unified-exp01/NA-RawUnified-D10M02Y24.tsv", sep = "\t", header = TRUE)

params_table <- read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Params-exp01.tsv", sep = "\t", header = TRUE)

# Cargar paquetes
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)

# Remover columna de IDs y modificar valores > 0
df_mutate <- data %>%
  mutate(across(-1, ~ ifelse(. > 0, 1, .)))

# Seleccionar última columna y renombrarla
df_selected <- df_mutate %>%
  select(1, last_col()) %>%
  rename("NA-counts" = last_col())

#--------------------------
# Unir con params_table y redondear columnas
df_joined <- df_selected %>%
  left_join(params_table %>% select(ID_simulation, Prob_neg, Prob_0, N_specs), 
            by = c("TSV_ID" = "ID_simulation")) %>%
  mutate(across(c(Prob_neg, Prob_0), round, 1))

# Code for grid method
df_joined <- df_selected %>%
  left_join(params_table %>% select(id, p_neg, p_noint, n_species), 
            by = c("TSV_ID" = "id")) %>%
  mutate(across(c(p_neg, p_noint), ~ round(.x, 1)))

#-------------------------------
# Sumar y contar valores por grupo
df_summed <- df_joined %>%
  group_by(Prob_neg, Prob_0) %>%
  summarise(NA_Simulations = sum(`NA-counts`, na.rm = TRUE), .groups = "drop")

counts_df <- df_joined %>%
  count(Prob_neg, Prob_0, name = "Total_Sims", sort = TRUE)

# Calcular `adjusted_NA` y unir conteo
df_result <- df_summed %>%
  left_join(counts_df, by = c("Prob_neg", "Prob_0")) %>%
  mutate(NA_SimsProps = NA_Simulations / Total_Sims)

p <- ggplot(df_result, aes(x = Prob_neg, y = Prob_0, fill = NA_SimsProps)) +
  geom_tile() +
  geom_text(aes(label = Total_Sims), color = "black", size = 6) +  # Add num_counts as text
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Proportions HeatMap of NA Counts by parameters", x = "Negative Probability", y = "No Interaction Probability", fill = "NAs") +
  theme(plot.title = element_text(hjust = 0.5))  # Center title
 
# Convert ggplot to an interactive plot
interactive_p <- ggplotly(p)

# Save as an interactive HTML file
library(htmlwidgets)
saveWidget(interactive_p, "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/Params/Exp01-C0CN", selfcontained = FALSE)

