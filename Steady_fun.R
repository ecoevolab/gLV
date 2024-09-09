
#----------------------------------------------Search for individual steady state--------------------------------------------------

SS_individual = function(ID, tol, wd) {
  
  #-----------------------Read table--------------------#
  library(data.table)
  out_path <- paste(wd, "./Outputs/O_", ID , ".tsv", sep = "") # output
  out <- fread(out_path, sep = "\t")  # Read TSV file

  specs <- nrow(out) # Species number
  Stb_vec <- numeric() # Empty numerical vector
  tol <- log(tol^2) # Apply transformation to tolerance
  
  #----------------------------Search Stability-----------------------#
  for (s in 1:specs) {
    V_spec <- as.numeric(out[s,]) # Species row
    V_spec <- diff(V_spec)^2 # Get differences and square them
    V_spec <- ifelse(V_spec == 0, NA, log(V_spec) )
                       
    Stab_col <- which( V_spec < tol, arr.ind = TRUE)[1] # Generation where the mean<tolerance
    Stb_vec <- c(Stb_vec, Stab_col) # Add vector
  }
  
  return(Stb_vec)
  
}

                    #-----------------------Moving average INDIVIDUAL---------------------------#

Rwindow_individual <- function(ID, tol, wd) {
  
  # Load necessary packages
  library(data.table)
  library(zoo)
  
  #-----------------------------Read table----------------------------#
  out_path <- paste(wd, "./Outputs/O_", ID , ".tsv", sep = "") # output
  out_table <- as.data.table ( fread(out_path, sep = "\t") ) # Read TSV file
  
  specs <- nrow(out_table)  # Species
  times <- ncol(out_table)  # Generations
  Stb_vec <- numeric(specs)  # Pre-allocate vector
  
  #----------------------------Search Stability-----------------------#
  window_size <- round(times * 0.1)  # Calculate window size
  
  # Function to calculate stability for one row
  find_stability <- function(row, window_size, tol) {
    moving_avg <- rollmean(row, window_size, fill = NA, align = "left") # Calculate moving average
    first_stable <- which(moving_avg < tol)[1]  # Find the first generation where moving average < tolerance
    return(first_stable)
  }
  
  # Apply the stability function to each species row
  for (s in 1:specs) {
    row <- as.numeric( out_table[s,] )
    Stb_vec[s] <- find_stability(row, window_size, tol)  # Assign directly to pre-allocated vector
  }
  
  names(Stb_vec) <- paste0("Specie", seq_len(specs))
  
  return(Stb_vec)
}

    
#----------------------------------------------Search for ALL steady state--------------------------------------------------
  
diff_SS <- function (out, tol) {
  
  #------------------------Get differences------------------------#
  gens <- ncol(out) # Times
  specs <- nrow(out) # Species number
  tol <- log(tol^2) # Apply transformation to tolerance
  
  # Compute differences and square them
  Stb_mat <- rbind() # Empty matrix
  for (r in 1:specs) {
    v <- as.vector(out[r,])
    Stb_mat <- rbind(Stb_mat, diff(v)^2) 
  }
  
  # Apply log transformation, replace 0 with NA
  ln_mat <- log(ifelse(Stb_mat == 0, NA, Stb_mat))
  
  # Assign row and column names
  rownames(ln_mat) <- paste("specie", 1:specs, sep = "")
  colnames(ln_mat) <- seq(1, gens-1)
  
  #---------------------------Create data frame------------------------#
  Stb_mean <- colMeans(Stb_mat, na.rm = TRUE) #Column means
  first_col <- which(ln_mat < tol, arr.ind = TRUE)[1] # Generation where the mean<tolerance
  
  return(list(Means = Stb_mean,
              Stable = first_col)
  )
}

#---------------------------------Moving average ALL-------------------------#

Rwindow_ALL <- function(out, tol) {
  
  library(zoo)
  times <- ncol(out)  # Number of generations
  
  #----------------------------Calculate Column Means-----------------#
  col_means <- colMeans(out, na.rm = TRUE)  # Compute column means
  
  #----------------------------Apply Moving Average---------------------#
  window_size <- round(times * 0.1)  # Calculate window size
  
  # Function to calculate moving average
  # "left" covers following rows 
  moving_avg <- rollmean(col_means, window_size, fill = NA, align = "left")
  
  # Find stability based on moving average
  Stb_vec <- which(moving_avg < tol)[1]  # Find the first generation where moving average < tolerance
  
  return(list(Stable = Stb_vec,
              Means = col_means
  )
  )
}

#-------------------------------------------Saver-------------------------------------------------------------

All_SS_save <- function (ID, tol, params, wd) {
  
  #-------------------------Read table----------------------------#
  # Load required package
  library(data.table)
  
  # Read table
  out_path <- paste(wd, "/Outputs/O_", ID , ".tsv", sep = "") # output
  out <- as.matrix( fread(out_path, sep = "\t") )
  
  gens <- ncol(out) # Times
  specs <- nrow(out) # Species number
  
  #--------------------------Test functions---------------------#
  result1 <- diff_SS(out, tol) # Differences^2 method
  result2 <- Rwindow_ALL(out, tol) # Rolling window method
  
  SS_df <- data.frame(
    'ID' = ID,
    'Population_seed' = params$Semilla[1],
    'Interaction_seed' = params$Semilla[2],
    'Growth_seed' = params$Semilla[3],
    'Generations' = gens,
    'Diff_ALL' = ifelse(is.na(as.numeric(result1$Stable)), "Not found", as.numeric(result1$Stable)),
    'Rwindow_ALL' = ifelse(is.na(as.numeric(result2$Stable)), "Not found", as.numeric(result2$Stable)) , 
    'Tolerance' = tol,
    'Individual' = FALSE
  )
  
  #---------------------------Save generation-------------------------#
  SS_path <- paste(wd, "/Scan/SS_all_new", ".tsv", sep = "") #Parameters
  exist <- file.exists(SS_path) # Bandera
  
  if (!exist) { # File doesnt exist
    
    file.create(SS_path) # Make file
    write.table(SS_df, file = SS_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
  } else {
    
    SS_table <- read.delim(SS_path, sep = "\t", header = TRUE) # Read table
    Join_ss <- rbind(SS_table, SS_df) # Join tables
    write.table(Join_ss, file = SS_path, sep = "\t", row.names = FALSE, col.names = TRUE) # Save
  }
  
  return(list(Window_method = result2$Stable,
              Diff_method = result1$Stable)
  )

}


#--------------------------Testing--------------------------
# ID <- "e4cffa"
# tol <- 0.05
# wd <-"/home/rivera/Cluster/"
# 
# res1 <- SS_individual(ID, tol, wd)
# res2 <- Rwindow_individual(ID, tol, wd)
# cat("Roll window method found the Steady State at generation", R_ss$Roll_window)
# cat("Differences^2 method found the Steady State at generation", R_ss$Method_dif)


  





