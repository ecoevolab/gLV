

st_search = function(ID, out, tol, individual) {
  
  # Leer tabla
  library(readr)
  test <- read.csv(out, sep = "\t")
  
  # Buscar estabilidad
  gens <- ncol(test) # Obtener el numero de comparaciones
  specs <- nrow(test) # Obtener el numero de comparaciones
  Stb_vec <- numeric() # Vector numerico vacio
  Stb_mat <- rbind()
    
  for (s in 1:specs) {
    for (g in 1:(gens-1)) {
        x <- test[s,g+1] - test[s,g]
        Stb_vec <- c(Stb_vec, x)
    }
    Stb_mat <- rbind(Stb_mat, Stb_vec) # Crear matriz de estabilidad
  }
  
  Stb_mat <- Stb_mat^2 # Elevar al cuadrado la matriz
  ln_mat <- ifelse(Stb_mat == 0, NA, log(Stb_mat)) # Remover ceros del df
  
  # Dar formato
  labels_r <- paste("specie", 1:specs, sep = "")
  rownames(Stb_mat) <- labels_r
  #Com_lab <- paste("Cr_", 2:(gens), "-", 1:(gens-1), sep = "")
  colnames(Stb_mat) <- seq(1, gens-1)
  
    if (individual) {
      
      # Inicializar variables
      steady  <- FALSE
      row_index <- 1
      sss_vector <- c() #Search for Steady State
      
      while (!steady) {
        
        # Define the condition 
        condition <- ln_mat[row_index, ] < log(tol)
        
        # Find the first column where the condition is met
        first_col <- which(condition)[1]
        
        # Add first column to a vector
        sss_vector <- c(sss_vector, first_col)
        
        # Next row
        row_index <- row_index + 1
        
        # Exit condition
        if (row_index > nrow(ln_mat)) {
          steady <- TRUE
        }
      }
      
      return(list(Stb_mat = Stb_mat, sss_vector = sss_vector))
      
    } else {
      Stb_mean <- colMeans(ln_mat) # Obtener promedios 
      condition <- Stb_mean < log(tol)    # Define the condition 
      first_col <- which(condition)[1] # Find the first column where the condition is met
      
      return(list(Stb_mat = Stb_mat, Stb_mean = Stb_mean, Fc = first_col ))
    }
  
    
}


# Revisar si el df contiene valores menores a 0
any_non_positive <- any(Stb_mat < 0) 
print(any_non_positive)

#----------------------------------------------Search for individual steady state--------------------------------------------------
library(dplyr)
library(purrr)
library(plotly)
library(reshape2)
library(ggplot2)

setwd("~/Documents/LAB_ECO") # Set Working Directory
#ID <- "2c9a4c"
ID <- "e03e07" #ID output
out <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output


individual  <- TRUE
tol <- 5
result <- st_search(ID,out,tol,individual)
  
Stb_mat <- result$Stb_mat
sss_vector <- result$sss_vector

  #-------------------- Graficar-------------#
  melted_df <- melt(Stb_mat)
  colnames(melted_df) <- c("Species", "Comparison", "Value")
  
  # Create the plot
  ggplot(melted_df, aes(x = Comparison, y = Value, color = Species, group = Species)) +
    geom_line() +
    theme_minimal() +
    labs(title = "Line Plot of Comparisons vs. Values by Species", x = "Comparisons", y = "Values") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
#----------------------------------------------Search for ALL steady state--------------------------------------------------

individual  <- FALSE
tol <- 5
result <- st_search(ID,out,tol,individual)
  
Stb_mat <- result$Stb_mat
Stb_mean <- result$Stb_mean
St_time <- result$Fc  

  #-------------------- Graficar-------------#
  # Create a data frame from the vector
  df <- data.frame(Comparison = 1:length(Stb_mean), Value = Stb_mean)

  # Create the plot
  ggplot(df, aes(x = Comparison, y = Value)) +
    geom_line() +
    geom_point() +  # Add points to the plot
    theme_minimal() +
    labs(title = "Sum of Values by Comparison", x = "Comparison", y = "Sum of Values") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlim(0,500)
