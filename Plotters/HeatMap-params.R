
#------------Load master table#------------#
counts <- read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Unified-exp01/NA-RawUnified-D10M02Y24.tsv", sep = "\t", header = TRUE)

params_table <- read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Params-exp01.tsv", sep = "\t", header = TRUE)

# Cargar paquetes
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)

# Remove IDs column, modify values greater than 0, and rename the last column
df_result <- counts  %>%
  mutate(across(-1, ~ ifelse(. > 0, 1, 0))) %>%   # Make values 0's or 1's depending in the abscence or presence of NAs.
  select(1, 2) %>% # Select last column
  rename("NA-counts" = last_col()) %>% # Rename last column
  # Join with params_table
  left_join(params_table %>% select(ID_simulation, Prob_neg, Prob_0, N_specs), by = c("TSV_ID" = "ID_simulation")) %>%
  mutate(across(c(Prob_neg, Prob_0), ~ round(.x, 1))) %>%  # Make double decimals the parameters
  group_by(Prob_neg, Prob_0) %>%  # Group by parameters
  summarise(
    Total_Sims = n(),  # Count the number of simulations
    NA_Simulations = sum(`NA-counts`, na.rm = TRUE),  # Sum the "NA-counts"
    .groups = "drop"  # Drop the grouping after summarising
  ) %>% 
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
saveWidget(interactive_p, "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/Params/Exp01-C0CN.html", selfcontained = FALSE)

