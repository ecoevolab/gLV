
#----------------------------------------------Search for individual steady state--------------------------------------------------

SS_individual = function(ID, out, tol, individual) {
  
  # Leer tabla
  library(readr)
  test <- read.csv(out, sep = "\t")
  
  # Buscar estabilidad
  gens <- ncol(test) # Times
  specs <- nrow(test) # Species number
  Stb_vec <- numeric() # Empty numerical vector
  Stb_mat <- rbind() # Empty matrix
  
  for (s in 1:specs) {
    for (g in 1:(gens-1)) {
      x <- test[s,g+1] - test[s,g] # Difference t+1-t
      Stb_vec <- c(Stb_vec, x) # Save difference
    }
    Stb_mat <- rbind(Stb_mat, Stb_vec) # Create steady state matrix
  }
  
  Stb_mat <- Stb_mat^2 # Square matrix
  ln_mat <- ifelse(Stb_mat == 0, NA, log(Stb_mat)) # Remove 0's from df
  
  # Dar formato
  rownames(Stb_mat) <- paste("specie", 1:specs, sep = "") # Assign rownames
  #Com_lab <- paste("Cr_", 2:(gens), "-", 1:(gens-1), sep = "")
  colnames(Stb_mat) <- seq(1, gens-1) # Assign colnames
  
  steady  <- FALSE # Flag
  row_index <- 1 # Row to start
  sss_vector <- c() # Empty vector
  
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
  
}

library(dplyr)
library(purrr)
library(plotly)
library(reshape2)
library(ggplot2)

setwd("~/Documents/LAB_ECO") # Set Working Directory
#ID <- "2c9a4c"
ID <- "e03e07" #ID output
out_path <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output


individual  <- TRUE
tol <- 5
result <- st_search(ID,out_path,tol,individual)
  
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

  
SS_all <- function(ID, out_path, tol) {
    
    # Load required package
    library(data.table)
    
    # Read table
    test <- as.matrix( fread(out_path, sep = "\t") )
    
    #------------------------Get differences------------------------#
    # Extract the number of columns and rows
    gens <- ncol(test) # Times
    specs <- nrow(test) # Species number
    
    # Compute differences and square them
    Stb_mat <- rbind() # Empty matrix
    for (r in 1:specs) {
      v <- as.vector(test[r,])
      Stb_mat <- rbind(Stb_mat, diff(v)^2) 
    }
    
    # Apply log transformation, replace 0 with NA
    ln_mat <- log(ifelse(Stb_mat == 0, NA, Stb_mat))
    
    # Assign row and column names
    rownames(ln_mat) <- paste("specie", 1:specs, sep = "")
    colnames(ln_mat) <- seq(1, gens-1)
    
    #---------------------------Create data frame------------------------#
    
    Stb_mean <- colMeans(Stb_mat, na.rm = TRUE) #Column means
    first_col <- which(Stb_mat < log(tol), arr.ind = TRUE)[1, 2] # Generation where the mean<tolerance
    
    tmp_df <- data.frame(
      ID = ID,
      Steady_generation = first_col,
    )
    
    #---------------------------Save generation-------------------------#
    SS_path <- paste("./Scan/SSS_all", ".tsv", sep = "") #Parameters
    exist <- file.exists(SS_path) # Bandera
    
    if (!exist) { # File doesnt exist
      
      file.create(SS_path) # Make file
      write.table(tmp_df, file = SS_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
    } else {
      
      SS_table <- read.delim(SS_path, sep = "\t", header = TRUE) # Read table
      Join_CPr <- rbind(SS_table, tmp_df) # Join tables
      write.table(Join_CPr, file = SS_path, sep = "\t", row.names = FALSE, col.names = TRUE) # Save
    }
      
    return(list(Stb_mat = Stb_mat, Stb_mean = Stb_mean, Fc = first_col))
}
  
tol <- 0.5
out_path <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output
result <- SS_all(ID,out_path,tol)
  
#--------------------Results------------------------------#
Stb_mat <- result$Stb_mat
Stb_mean <- data.frame(result$Stb_mean)
St_time <- result$Fc  

###########################################################################################################
  #--------------Save Deltas----------------#
  unlink("~/Documents/LAB_ECO/Scan/*", recursive = TRUE, force = TRUE) # Remove all data from Scan 
    
  
  st_all_save = function(ID, tol, Stb_mean) {
    
      library(readr)
      
      stb_path <- paste("./Scan/Stb_ALL", ".tsv", sep = "") # output
      exist <- file.exists(stb_path) # Bandera
      
      if (!exist) { # File doesnt exist
        
        system(paste("touch", stb_path)) # Make file
        
        tmp <- data.frame(id = ID, tolerance = tol) 
        tmp_df <- cbind(tmp, time = t(Stb_mean)) # Add mean info
        
        
        write.table(tmp_df, file = stb_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
      } else {
        
        Rd_table <- read.delim(stb_path, sep = "\t", header = TRUE) # Read table
        
        tmp <- data.frame(id = ID, tolerance = tol)
        tmp_df <- cbind(tmp, time = t(Stb_mean)) # Add mean info
        rownames(tmp_df) <- NULL # Eliminate row names
        Join_tmp <- rbind(Rd_table, tmp_df) # Join tables
        
        write.table(Join_tmp, file = stb_path, sep = "\t", row.names = FALSE, col.names = TRUE) # Save
      }
  }

  st_all_save(ID, tol, Stb_mean)

#-------------------- Graficar--------------------------------------------------------------
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
    
    
    
#---------------------------------EXTRA---------------------------------
    
    # Revisar si el df contiene valores menores a 0
    any_non_positive <- any(Stb_mat < 0) 
    print(any_non_positive)
    
    
    
    
    
    
    
    
    
    
    
    #--------------------- Save Deltas-----------------------------#
    library(readr)
    stb_path <- paste("./Scan/Stb_ALL", ".tsv", sep = "") # output
    exist <- file.exists(stb_path) # Bandera
    
    if (!exist) { # File doesnt exist
      
      system(paste("touch", stb_path)) # Make file
      
      tmp <- data.frame(id = ID, tolerance = tol) 
      tmp_df <- cbind(tmp, time = t(Stb_mean)) # Add mean info
      
      write.table(tmp_df, file = stb_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
      
    } else {
      
      Rd_table <- read.delim(stb_path, sep = "\t", header = TRUE) # Read table
      
      tmp <- data.frame(id = ID, tolerance = tol)
      tmp_df <- cbind(tmp, time = t(Stb_mean)) # Add mean info
      rownames(tmp_df) <- NULL # Eliminate row names
      
      Join_tmp <- rbind(Rd_table, tmp_df) # Join tables
      write.table(Join_tmp, file = stb_path, sep = "\t", row.names = TRUE, col.names = TRUE) # Save
    }