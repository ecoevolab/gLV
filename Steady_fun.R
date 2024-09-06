
#----------------------------------------------Search for individual steady state--------------------------------------------------

SS_individual = function(ID, out, tol, individual) {
  
  #-----------------------Read table--------------------#
  library(readr)
  tmp <- read.csv(out, sep = "\t")
  specs <- nrow(tmp) # Species number
  Stb_vec <- numeric() # Empty numerical vector
  tol <- log(tol^2)
  
  #----------------------------Search Stability-----------------------#
  for (s in 1:specs) {
    V_spec <- as.vector(tmp[s,]) # Species row
    V_spec <- diff(v)^2 # Get differences and square them
    ln_Vspec <- ifelse(V_spec == 0, NA, log(V_spec) )
                       
    Stab_col <- which( ln_Vspec < log(tol), arr.ind = TRUE)[1] # Generation where the mean<tolerance
    Stb_vec <- c(Stb_vec, Stab_col) # Add vector
  }
  
  return(Stb_vec)
  
}


    
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
      
    return(list(Stb_mat = Stb_mat, 
                Stb_mean = Stb_mean, 
                Fc = first_col)
           )
}

#------------------------------------------Moving average INDIVIDUAL----------------------------------------------------

Rwindow_indiv <- function(ID, tol) {
  # Load necessary packages
  library(data.table)
  library(zoo)
  
  #-----------------------------Read table----------------------------#
  out_path <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output
  tmp <- fread(out_path, sep = "\t")  # Read TSV file
  tmp <- as.data.table(tmp)  # Ensure tmp is a data.table
  specs <- nrow(tmp)  # Number of species
  times <- ncol(tmp)  # Number of generations
  Stb_vec <- numeric(specs)  # Pre-allocate vector
  
  #----------------------------Search Stability-----------------------#
  window_size <- round(times * 0.05)  # Calculate window size
  
  # Function to calculate stability for one row
  find_stability <- function(row) {
    Row_tmp <- rollmean(row, window_size, fill = NA, align = "right")  # Calculate moving average
    first_stable <- which(Row_tmp < tol)[1]  # Find first stable point
    return(first_stable)
  }
  
  # Apply the stability function to each species row
  Stb_vec <- apply(tmp, 1, find_stability)
  names(Stb_vec) <- paste0("Specie", seq_len(specs))
  return(Stb_vec)
}

#--------------------------------#
# Read the data
out_path <- paste("./Outputs/O_", ID, ".tsv", sep = "")  # Output path
tmp <- fread(out_path, sep = "\t")  # Read TSV file

times <- ncol(tmp)  # Number of generations
window_size <- round(times * 0.05)  # Calculate window size

V_spec <- as.vector(unlist(tmp[2, ])) # Convert to vector

result1 <- Rwindow_indiv(ID,tol)
x <- result1[2]
y <- x + window_size

if (mean(V_spec[x:y]) <= tol) {
  print("Test succeded")
}

#---------------------------------------------------Moving average ALL--------------------------------------------------
Rwindow_ALL <- function(ID, tol) {
  # Load necessary packages
  library(data.table)
  library(zoo)
  
  #-----------------------------Read table----------------------------#
  out_path <- paste("./Outputs/O_", ID , ".tsv", sep = "")  # Output path
  tmp <- fread(out_path, sep = "\t")  # Read TSV file
  tmp <- as.data.table(tmp)  # Ensure tmp is a data.table
  times <- ncol(tmp)  # Number of generations
  
  #----------------------------Calculate Column Means-----------------#
  col_means <- colMeans(tmp, na.rm = TRUE)  # Compute column means
  
  #----------------------------Apply Moving Average---------------------#
  window_size <- round(times * 0.05)  # Calculate window size
  
  # Function to calculate moving average
  moving_avg <- rollmean(col_means, window_size, fill = NA, align = "right")
  
  # Find stability based on moving average
  Stb_vec <- which(moving_avg < tol)[1]  # Find the first generation where moving average < tolerance
  
  # Assign name to the stability vector
  names(Stb_vec) <- "Stability_Point"
  
  return(list(All_vec = Stb_vec,
              Means = col_means
              )
  )
}


result <- Rwindow_ALL(ID, tol)
x <- result$All_vec

if (result$Means[x] <= tol) {
  print("Test succeded")
}
