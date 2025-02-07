
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

# Join the two data frames by the 'id' column
df_joined <- df_selected %>%
  left_join(params_table %>% select(ID_simulation, Prob_neg, Prob_0, N_specs), 
            by = c("TSV_ID" = "ID_simulation"))

# Assuming df_joined is your data frame
library(tidyr)

df_long <- df_joined %>%
  select(TSV_ID, Prob_neg, Prob_0, N_specs) %>%
  pivot_longer(cols = c("Prob_neg", "Prob_0"),
               names_to = "Probability_Type",
               values_to = "Value")


# Create HeatMap
library(ggplot2)
p <- ggplot(df_long, aes(x = Prob_neg, y = Prob_0, fill = Value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(
    plot.title = element_text(hjust = 0.5)  # This centers the title
  ) +
  labs(title = "NA by parameters", x = "Negative probability", y = "No interaction prob")

# Convert ggplot to an interactive plot
interactive_p <- ggplotly(p)

