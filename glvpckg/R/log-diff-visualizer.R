

visualizer_logdiff <- function(wd) {
  
  # Log differences path
  mdiff_path <- file.path(wd, "Differences", paste0("means_ld", ".tsv"))
  data <- read.table(mdiff_path, sep = "\t", header = FALSE)
  
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("plotly package is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required but not installed.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("tidyr package is required but not installed.")
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("reshape2 package is required but not installed.")
  }
  
  # Calculate species number and generations
  gens <- ncol(data)
  colnames(data) <- c("ID", seq_len(gens - 1) )
  
  #-------------------- Interactive Plot using Plotly -------------#
  melted_df <- reshape2::melt(data)
  colnames(melted_df) <- c("ID", "Generation", "Value")
  
  df_long <- reshape2::melt(data, id.vars = "ID", variable.name = "Generation", value.name = "Value")
  
  # Create the plot
  alternate_colors <- c("5fe0a0" = "blue", "44d2c3" = "red")
  # Create the plot
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = as.numeric(Generation), y = Value, group = ID)) +
    ggplot2::geom_line(color = "#666666") +  # Set all lines to grey
    # ggplot2::scale_color_brewer(palette = "Set1") +  # Use a Brewer palette for distinct colors
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Difference means", x = "Generations", y = "Population") +
    ggplot2::geom_hline(yintercept = log(1.001), color = "red", linetype = "solid") + 
    ggplot2::geom_hline(yintercept = -log(1.001), color = "red", linetype = "solid") +
    ggplot2::coord_cartesian(xlim = c(50, max(as.numeric(df_long$Generation))), ylim = c(-log(1.01), log(1.01))) +  # Set x and y limits
    ggplot2::scale_x_continuous(breaks = seq(from = 50, to = max(as.numeric(df_long$Generation)), by = 100))  # Set x-axis ticks every 100
  
  # Make the ggplot interactive using ggplotly
  interactive_plot <- plotly::ggplotly(p)
  
  return(interactive_plot)
}
  



























