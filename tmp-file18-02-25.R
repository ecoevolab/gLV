
#------------Load master table#------------#
data_miasim <- read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/D13M02Y25-miaSim/NAcount-D13M02Y25-miaSim.tsv", sep = "\t", header = TRUE)

data_ode <-  read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/D13M02Y25/NAcount-D13M02Y25.tsv", sep = "\t", header = TRUE)

params_table <- read.delim("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/D13M02Y25.tsv", sep = "\t", header = TRUE)

#--------------------------------------HeatMap miaSim method-----------------------#
library(dplyr)
library(ggplot2)
library(plotly)

# Function to process and summarize data
process_data <- function(data, params_table) {
  data %>%
    mutate(across(-1, ~ if_else(. > 0, 1, .))) %>% # Convert to 0's or 1's
    # Join parameters
    left_join(params_table %>% select(id, p_neg, p_noint), by = c("TSV_ID" = "id")) %>%
    mutate(across(c(p_neg, p_noint), ~ round(.x, 2))) %>%
    group_by(p_neg, p_noint) %>%
    summarise(
      NA_Simulations = sum(NA_Counts, na.rm = FALSE),
      Total_Sims = n(),
      NA_SimsProps = NA_Simulations / Total_Sims,
      .groups = "drop"
    )
}

# Apply function to both datasets
summed_miasim <- process_data(data_miasim, params_table)
summed_ode <- process_data(data_ode, params_table)

# Function to create a heatmap
create_heatmap <- function(data) {
  ggplot(data, aes(x = p_neg, y = p_noint, fill = NA_SimsProps,
                   text = paste("Total Sims:", Total_Sims,
                                "<br>p_neg:", p_neg,
                                "<br>p_noint:", p_noint,
                                "<br>NA Sims:", NA_Simulations,
                                "<br>Proportions:", NA_SimsProps))) +
    geom_tile() +
    scale_fill_distiller(palette = "RdBu", direction = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5))
}

# Generate heatmaps
p1 <- create_heatmap(summed_miasim)
p2 <- create_heatmap(summed_ode)

# Convert to interactive plots
p1_interactive <- ggplotly(p1, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))
p2_interactive <- ggplotly(p2, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))

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

