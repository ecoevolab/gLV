
#---------------------------------------Graficador de todas las especies-----------------------------------------

All_visualizer <- function(ID, wd) {
  
  #library(dplyr)
  #library(purrr)
  library(plotly)
  library(reshape2)
  library(ggplot2)
  library(data.table)
  
  out_path <- paste(wd, "/Outputs/O_", ID , ".tsv", sep = "") # output path
  out <- as.matrix( fread(out_path, sep = "\t") )
  
  # Calculate species number and generations
  specs <- ncol(out)
  gens <- nrow(out)
  
  # Assign row and column names
  rownames(out) <- paste("specie", 1:specs, sep = "")
  colnames(out) <- seq(1, gens)
  
  
  #-------------------- Interactive Plot using Plotly -------------#
  melted_df <- melt(out)
  colnames(melted_df) <- c("Species", "Generation", "Value")
  
  # Create the plot
  p <- ggplot(melted_df, aes(x = Generation, y = Value, color = Species, group = Species)) +
    geom_line() +
    theme_minimal() +
    labs(title = "Line Plot of Comparisons vs. Values by Species", x = "Generations", y = "Population") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Make the ggplot interactive using ggplotly
  interactive_plot <- ggplotly(p)
  
  return(interactive_plot)
}

#--------------------------------How to run it-----------------------#
ID <- "e4cffa"
wd <-"/home/rivera/Cluster/"
All_visualizer (ID, wd)


#--------------------------------------------Graficador de los promedios-----------------------------------------------
Mean_visualizer <- function(ID, wd) {
  
    #library(dplyr)
    #library(purrr)
    library(plotly)
    library(reshape2)
    library(ggplot2)
    library(data.table)
    
  #-------------------------------------Read table and format it-------------------------#
    out_path <- paste(wd, "/Outputs/O_", ID , ".tsv", sep = "") # output path
    out <- as.matrix( fread(out_path, sep = "\t") )
    col_means <- colMeans(out, na.rm = TRUE)  # Compute column means
   
    # Create a data frame with 'Generation' and 'Value'
    col_means_df <- data.frame(
      Generation = seq(1, length(col_means)),  # Assign generation numbers (1, 2, 3, ...)
      Value = col_means  # Assign column means as the 'Value'
    )
    
    #-------------------- Interactive Plot using Plotly -------------#
    # Create the plot
    p <- ggplot(col_means_df, aes(x = Generation, y = Value)) +
      geom_line(color = "red")
      theme_minimal() +
      labs(title = "Line Plot of Comparisons vs. Values by Species", x = "Generations", y = "Population") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Make the ggplot interactive using ggplotly
    interactive_plot <- ggplotly(p)
    
    return(interactive_plot)
}

#--------------------------------How to run it-----------------------#
ID <- "e4cffa"
wd <-"/home/rivera/Cluster/"
Mean_visualizer (ID, wd)

