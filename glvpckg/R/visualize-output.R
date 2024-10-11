#' Plot the Output of All Species
#'
#' This function plots the change in population over time (generations).
#'
#' @param output A data frame or matrix containing the simulation output to be visualized. The columns represent generations, and the rows represent species.
#' 
#' @return A plot where the x-axis represents generations and the y-axis represents population. Each line corresponds to a species within the simulation.
#' 
#' @import ggplot2
#' @import plotly
#' @import tidyr
#' @import utils
#'
#' @examples
#' # Example usage
#' wd <- "~/Documents/LAB_ECO/testing"
#' uniqueID <- "bcfg45"
#' 
#' out_path <- paste(wd, "/Outputs/O_", uniqueID, ".tsv", sep = "") # Output path
#' output <- as.matrix(data.table::fread(out_path, sep = "\t"))
#' visualize_output(output)
#' @export


visualize_output <- function(output) {
  
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
  requireNamespace("utils", quietly = TRUE)
  
  # Declare global variables
  utils::globalVariables(c("Generation", "Value", "Species"))
  
  # Calculate species number and generations
  specs <- nrow(output)
  gens <- ncol(output)
  
  # Assign row and column names
  rownames(output) <- paste("Species", 1:specs, sep = "")
  colnames(output) <- seq_len(gens)
  
  #-------------------- Interactive Plot using Plotly -------------#
  melted_df <- reshape2::melt(output)
  colnames(melted_df) <- c("Species", "Generation", "Value")
  
  # Create the plot
  p <- ggplot2::ggplot(data = melted_df, ggplot2::aes(x = Generation, y = Value, color = Species, group = Species)) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Line Plot of Comparisons vs. Values by Species", x = "Generations", y = "Population") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # Make the ggplot interactive using ggplotly
  interactive_plot <- plotly::ggplotly(p)
  
  return(interactive_plot)
}


