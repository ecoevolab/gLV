
#------------Load master table#------------#
data <- read.delim("/mnt/atgc-d3/sur/users/mrivera/testing/Sims_Exp02/NA-counts-exp02.tsv", sep = "\t", header = TRUE)

params_table <- read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Params-exp01.tsv", sep = "\t", header = TRUE)

# Remover columna de IDs
ids <- data[,1]


# Modify all columns where values are greater than 0
library(dplyr)
df_mutate <- data %>%
  mutate(across(-1, ~ ifelse(. > 0, 1, .)))


# Subset data
df_selected <- df_mutate %>%
  select(1, ncol(df_mutate))
df_selected <- df_selected |> rename("NA-counts" = `X1e.01.1e.01`)

# Join the two data frames by the 'id' column
df_joined <- df_selected %>%
  left_join(params_table %>% select(ID_simulation, Prob_neg, Prob_0, N_specs), 
            by = c("TSV_ID" = "ID_simulation"))

df_joined <- df_joined |> 
  mutate(across(c(Prob_neg, Prob_0), ~ round(.x, 1)))

df_summed <- df_joined |> 
  group_by(Prob_neg, Prob_0) |> 
  summarise(Total_NA_counts = sum(`NA-counts`, na.rm = TRUE))

# Create HeatMap
library(tidyr)
library(ggplot2)
library(plotly)

p <- ggplot(df_summed, aes(x = Prob_neg, y = Prob_0, fill = `Total_NA_counts`)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Heatmap of NA Counts", x = "Negative Probability", y = "No Interaction Probability", fill = "NAs") +
  annotation_custom(segmentsGrob(c(0.3, -0.1), c(-0.085, 0.28), 
                                 c(1, -0.1), c(-0.085, 1), gp = gpar(lwd = 2),
                                 arrow = arrow(length = unit(2.5, 'mm'))))

# Convert ggplot to an interactive plot
interactive_p <- ggplotly(p)

# Save as an interactive HTML file
library(htmlwidgets)
saveWidget(interactive_p, "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/HeatMap_exp02-old.html", selfcontained = FALSE)

