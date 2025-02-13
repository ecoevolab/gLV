# Code for separating the NA counts table


#--------------------------Obtener los conteos------------------#

counts_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Unified/NA-RawUnified-D10M02Y24.tsv"
  
counts_table <- data.table::fread(counts_path)

# Remover columna de IDs
counts_table <- counts_table[,-1]

# Modify all columns where values are greater than 0
library(dplyr)
df_mutate <- counts_table %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, .)))

sumatory <- colSums(df_mutate)
df <- data.frame(matrix(sumatory, ncol = 6, byrow = TRUE ))  # Assuming 12 rows (you can change if needed)

vector <- names(sumatory)
tols <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1) 

# Rows correspond to absolute tolerance
rownames(df) <-   format(tols, scientific = TRUE, digits = 3)

# Columns correspond to relative tolerance
colnames(df) <-  format(tols, scientific = TRUE, digits = 3)

# Optional: Guardar la tabla
# counts_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/RawHeatMap-NA.tsv"
# data.table::fwrite(df, file = counts_path, sep = "\t" )

#--------------------------Convertirlo a Heat Map------------------#

# Reshape the data frame into long format
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(plotly)

#  df <-data.table::fread("/home/rivera/Cluster/testing/Sims_Exp02/HeatMap-counts.tsv")

df_long <- df %>%
  tibble::rownames_to_column(var = "Tolerance_R") %>%
  pivot_longer(cols = -Tolerance_R, names_to = "Tolerance_A", values_to = "Value")

# Create the heatmap with ggplot2
p <- ggplot(df_long, aes(x = Tolerance_A, y = Tolerance_R, fill = Value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(
    plot.title = element_text(hjust = 0.5)  # This centers the title
  ) +
  labs(title = "Raw ODE data", x = "Absolute tolerance", y = "Relative tolerance")

# Convert ggplot to an interactive plot
interactive_p <- ggplotly(p)

# Save as an interactive HTML file
saveWidget(interactive_p, "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/RawHeatMap-NA.html", selfcontained = FALSE)