
#------------Load master table#------------#
data_miasim <- read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/D13M02Y25-miaSim/NAcount-D13M02Y25-miaSim.tsv", sep = "\t", header = TRUE)

data_ode <-  read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/D13M02Y25/NAcount-D13M02Y25.tsv", sep = "\t", header = TRUE)

params_table <- read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/D13M02Y25.tsv", sep = "\t", header = TRUE)

# Cargar paquetes
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)

#--------------------------------------HeatMap miaSim method-----------------------#

# Remover columna de IDs y modificar valores > 0
df_miasim <- data_miasim %>%
  mutate(across(-1, ~ ifelse(. > 0, 1, .))) %>% # Do 0 or 1 if simulations had NA values
  # Join parameters by id
  left_join(params_table %>% select(id, p_neg, p_noint, n_species), by = c("TSV_ID" = "id")) %>% 
  mutate(across(c(p_neg, p_noint), ~ round(.x, 2))) # Make double decimals the parameters

# Sumar y contar valores por grupo
summed_miasim <- df_miasim %>%
  group_by(p_neg, p_noint) %>% # Group by parameters
  summarise(
    NA_Simulations = sum(NA_Counts, na.rm = TRUE),  # Count the sum of NA simulations
    Total_Sims = n(),  # Count total simulations
    .groups = "drop"
  ) %>%
  mutate(NA_SimsProps = NA_Simulations / Total_Sims) # Compute proportions

# Create ggplot object with interactive text
p1 <- ggplot(summed_miasim, aes(x = p_neg, y = p_noint, fill = NA_SimsProps, 
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
  theme(plot.title = element_text(hjust = 0.5))  # Center title

# Convert to interactive plotly with hover text
p1_interactive <- ggplotly(p1, tooltip = "text") %>%
  layout(hoverlabel = list(bgcolor = "white"))

#--------------------------------------HeatMap ODE method-----------------------#

# Remover columna de IDs y modificar valores > 0
df_ode <- data_ode %>%
  mutate(across(-1, ~ ifelse(. > 0, 1, .))) %>% # Do 0 or 1 if simulations had NA values
  # Join parameters by id
  left_join(params_table %>% select(id, p_neg, p_noint, n_species), by = c("TSV_ID" = "id")) %>% 
  mutate(across(c(p_neg, p_noint), ~ round(.x, 2))) # Make double decimals the parameters

# Sumar y contar valores por grupo
summed_ode <- df_ode %>%
  group_by(p_neg, p_noint) %>% # Group by parameters
  summarise(
    NA_Simulations = sum(NA_Counts, na.rm = TRUE),  # Count the sum of NA simulations
    Total_Sims = n(),  # Count total simulations
    .groups = "drop"
  ) %>%
  mutate(NA_SimsProps = NA_Simulations / Total_Sims) # Compute proportions

# Create ggplot object with interactive text
p2 <- ggplot(summed_ode, aes(x = p_neg, y = p_noint, fill = NA_SimsProps, 
                                text = paste("Total Sims:", Total_Sims,
                                             "<br>p_neg:", p_neg,
                                             "<br>p_noint:", p_noint,
                                             "<br>NA Sims:", NA_Simulations,
                                             "<br>Proportions:", NA_SimsProps
                                ))) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  theme(plot.title = element_text(hjust = 0.5))  # Center title

# Convert to interactive plotly with hover text
p2_interactive <- ggplotly(p2, tooltip = "text") %>%
  layout(hoverlabel = list(bgcolor = "white"))

#-----------------------------------#
# Combine both interactive plots in a single figure
combined_plot <- subplot(p1_interactive, p2_interactive, nrows = 1, shareY = TRUE, titleX = TRUE) %>%
  layout(
    annotations = list(
      list(x = 0.2, y = 1, text = "miaSim Method", showarrow = FALSE, xref='paper', yref='paper', font=list(size=16)),
      list(x = 0.8, y = 1, text = "ODE Method", showarrow = FALSE, xref='paper', yref='paper', font=list(size=16))
    ),
    xaxis = list(title = "Negative Probability"),
    xaxis2 = list(title = "Negative Probability"),
    yaxis = list(title = "No Interaction Probability"),
    yaxis2 = list(title = "Negative Probability"),
    title = list(title = "Some general title")
  )

# Save as interactive HTML if needed
htmlwidgets::saveWidget(combined_plot, "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/Params/combined testing.html", selfcontained = FALSE)










































library(ggplot2)
library(patchwork)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(plotly)

# Crear datos de ejemplo (matrices de 10x10)
set.seed(42)
heat1 <- matrix(runif(100), nrow=10)
heat2 <- matrix(runif(100), nrow=10)

# Convertir matrices a data frames en formato largo
df1 <- melt(heat1)
df2 <- melt(heat2)

# Crear heatmaps con ggplot2
p1 <- ggplot(df1, aes(Var1, Var2, fill=value)) + 
  geom_tile(alpha=0.5) + 
  scale_fill_viridis_c() 

p2 <- ggplot(df2, aes(Var1, Var2, fill=value)) + 
  geom_tile(alpha=0.5) + 
  scale_fill_viridis_c() 

p1_interactive <- ggplotly(p1) # %>% layout(title = "Heatmap 1")
p2_interactive <- ggplotly(p2) # %>% layout(title = "some title")

# Display both heatmaps side by side with titles
fig <- subplot(p1_interactive, p2_interactive, nrows = 1, titleY = TRUE, titleX = TRUE, margin = 0.1 )
fig <- fig %>% layout(annotations = list(
  list(x = 0.2, y = 1, text = "ODE method", showarrow = FALSE, xref='paper', yref='paper', xanchor = "center",  
       yanchor = "bottom", font=list(size=16)),
  list(x = 0.8, y = 1, text = "miaSim method", showarrow = FALSE, xref='paper', yref='paper', xanchor = "center",  
       yanchor = "bottom", font=list(size=16))
))
fig
