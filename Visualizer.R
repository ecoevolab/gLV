
#---------------------------------------Graficador de todas las especies-----------------------------------------
#library(dplyr)
#library(purrr)
library(plotly)
library(reshape2)
library(ggplot2)

ID <- "51ee5f"

setwd("~/Documents/LAB_ECO") # Set Working Directory
out_path <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output path
out <- as.matrix( data.table::fread(out_path, sep = "\t") )

specs <- ncol(out)
gens <- nrow(out)

# Assign row and column names
rownames(out) <- paste("specie", 1:specs, sep = "")
colnames(out) <- seq(1, gens)


#-------------------- Graficar-------------#
melted_df <- reshape2::melt(out)
colnames(melted_df) <- c("Species", "Generation", "Value")

# Create the plot
ggplot(melted_df, aes(x = Generation, y = Value, color = Species, group = Species)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Line Plot of Comparisons vs. Values by Species", x = "Comparisons", y = "Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#---------------------------------------Graficador de todas las especies-----------------------------------------