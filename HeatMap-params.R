
#------------Load master table#------------#
data <- read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Unified-exp01/NA-PropsUnified-D10M02Y24.tsv", sep = "\t", header = TRUE)

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

# Count the ocurrances per group
counts_df <-  df_joined %>% count(Prob_neg, Prob_0, sort = TRUE,  name = "num_counts")

# Join the count dataframe back to df_joined
df_result <- df_summed %>% left_join(counts_df, by = c("Prob_neg", "Prob_0")) %>%
  mutate(adjusted_NA = Total_NA_counts / num_counts)

# Create HeatMap
library(tidyr)
library(ggplot2)
library(plotly)

p <- ggplot(df_result, aes(x = Prob_neg, y = Prob_0, fill = `adjusted_NA`)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Heatmap of NA Counts", x = "Negative Probability", y = "No Interaction Probability", fill = "NAs")

# Convert ggplot to an interactive plot
interactive_p <- ggplotly(p)

# Save as an interactive HTML file
library(htmlwidgets)
saveWidget(interactive_p, "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/PropsHeatMap-Exp01-Params.html", selfcontained = FALSE)

