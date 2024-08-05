
setwd("~/Documents/LAB_ECO") # Set Working Directory

# Carpeta que contiene los Outputs
ID <- "2c9a4c"
out <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output


# Leer tabla
library(readr)
test <- read.csv(out, sep = "\t")

# Buscar estabilidad
gens <- round(ncol(test)) # Obtener el numero de comparaciones
specs <- round(nrow(test)) # Obtener el numero de comparaciones
Stb_vec <- numeric() # Vector numerico vacio
Stb_mat <- rbind()
  
for (s in 1:specs) {
  for (g in 1:(gens-1)) {
      x <- test[s,g+1] - test[s,g]
      Stb_vec <- c(Stb_vec, x)
  }
  Stb_mat <- rbind(Stb_mat, Stb_vec) # Crear matriz de estabilidad
}

# Elevar al cuadrado la matriz
Stb_mat <- Stb_mat^2

# Dar formato
labels_r <- paste("specie", 1:specs, sep = "")
rownames(Stb_mat) <- labels_r

Com_lab <- paste("Cr_", 2:(gens), "-", 1:(gens-1), sep = "")
colnames(Stb_mat) <- Com_lab


#---------------------------Individual plot-----------------------------------------
library(plotly)
library(reshape2)
library(ggplot2)

melted_df <- melt(Stb_mat)
colnames(melted_df) <- c("Species", "Comparison", "Value")

# Create the plot
ggplot(melted_df, aes(x = Comparison, y = Value, color = Species, group = Species)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Line Plot of Comparisons vs. Values by Species", x = "Comparisons", y = "Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############,,
ln_mat <- log10(Stb_mat)

melted_df <- melt(ln_mat)
colnames(melted_df) <- c("Species", "Comparison", "Value")

# Create the plot
ggplot(melted_df, aes(x = Comparison, y = Value, color = Species, group = Species)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Line Plot of Comparisons vs. Values by Species", x = "Comparisons", y = "Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#---------------------------Sum plot-----------------------------------------

# Sumar valores de todas las columnas
ln_mat <- log(Stb_mat)
Stb_sum <- colSums(ln_mat)
df <- data.frame(
  Comparison = names(Stb_sum),  # Column names are used as x-axis labels
  Value = Stb_sum
)

# Create the plot
ggplot(df, aes(x = Comparison, y = Value)) +
  geom_line() +
  geom_point() +  # Add points to the plot
  theme_minimal() +
  labs(title = "Sum of Values by Comparison", x = "Comparison", y = "Sum of Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


