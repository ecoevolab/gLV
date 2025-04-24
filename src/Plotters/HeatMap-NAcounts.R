# Code for separating the NA counts table


#--------------------------Obtener los conteos------------------#

dirpath <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/Unified-exp01"

props_counts <- read.delim( file.path(dirpath, "NA-PropsUnified-D10M02Y24.tsv"), sep = "\t", header = TRUE)
raw_counts <- read.delim( file.path(dirpath, "NA-RawUnified-D10M02Y24.tsv"), sep = "\t", header = TRUE)
  
library(dplyr)

# Remove IDs column and convert to 0's or 1's (for both data frames)
df1_props <- props_counts[,-1] %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))

df2_raw <- raw_counts[,-1] %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))

# Compute column sums
df1_sum <- colSums(df1_props)
df2_sum <- colSums(df2_raw)

# Convert to data frame
df1_sum <- data.frame(matrix(df1_sum, ncol = 6, byrow = TRUE )) 
df2_sum <- data.frame(matrix(df2_sum, ncol = 6, byrow = TRUE )) 

# Format row and column names for both data frames
tols <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
formatted_tols <- format(tols, scientific = TRUE, digits = 3)

# Format row and column names for both data frames
formatted_tols <- format(tols, scientific = TRUE, digits = 3)
rownames(df1_sum) <- formatted_tols
rownames(df2_sum) <- formatted_tols
colnames(df1_sum) <- formatted_tols
colnames(df2_sum) <- formatted_tols

#--------------------------Convertirlo a Heat Map------------------#

# Reshape the data frame into long format
library(ggplot2)
library(reshape2)
library(tidyr)
library(plotly)

# Function to create long-format data and plot the heatmap
create_heatmap <- function(df_sum) {
  # Convert to long data format
  df_long <- df_sum %>%
    tibble::rownames_to_column(var = "Tolerance_R") %>%
    pivot_longer(cols = -Tolerance_R, names_to = "Tolerance_A", values_to = "Value")
  
  # Create heatmap plot
  p <- ggplot(df_long, aes(x = Tolerance_A, y = Tolerance_R, fill = Value,
                           text = paste("NA_Count:", Value,
                                        "<br>Abs_tol:", Tolerance_A,
                                        "<br>Rel_tol:", Tolerance_R))) +
    geom_tile() +
    scale_fill_distiller(palette = "RdBu", direction = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  
  # Convert to interactive plot
  p_interactive <- ggplotly(p, tooltip = "text") %>%
    layout(hoverlabel = list(bgcolor = "white"))
  
  return(p_interactive)
}

# Create interactive heatmaps for both dataframes
p1_interactive <- create_heatmap(df1_sum)
p2_interactive <- create_heatmap(df2_sum)

# Combine both interactive plots in a single figure
combined_plot <- subplot(p1_interactive, p2_interactive, nrows = 1, shareY = TRUE, titleX = TRUE) |>
  add_annotations(
    x = c(.25, .75),
    y = 1,
    xref = "paper",
    yref = "paper",
    text = c("Proportions method, D10M02Y24 data", "Raw ODE method, D10M02Y24 data"),
    showarrow = F,
    xanchor = "center",
    yanchor = "bottom",
    font = list(size = 12)
  )

# Save as interactive HTML if needed
htmlwidgets::saveWidget(combined_plot, "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/Outputs/Raw-vs-Props-D10M02Y24.html", selfcontained = FALSE)
