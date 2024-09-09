


#----------------------------------------------------Function to save Output----------------------------------------------------------

save_test = function(N_specs, C0, CN, V_diag, output, params){
  
  library(ids)
  
  # Generar el ID
  ID <- ids::random_id(1, 3)
  
  # Especificar el path del output
  out_path <- paste("./test/Outputs/O_", ID , ".tsv", sep = "") # output
  
  # --------------------------------------------Save OUTPUT---------------------------------------------------------#
  
  # Revisar si un archivo con ese ID existe
  exist <- file.exists(out_path) # flag
  
  while (exist) { # TRUE-> file exist, change ID...
    
    # Generar el ID nuevo
    ID <- ids::random_id(1, 3)
    out_path <- paste("./test/Outputs/O_", ID , ".tsv", sep = "") # output
    exist <- file.exists(out_path) # No existe = FALSE
  }
  
  file.create(out_path) # Create files
  write.table(output, file = out_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  # --------------------------------------------Save Seeds---------------------------------------------------------#
  # Primera semilla = Poblaciones
  # Segunda semilla = Interactions
  # Tercera semilla = Growth rates
  
  seeds_df <- data.frame(
    ID_simulation = ID,
    N_specs = N_specs,
    Prob_0 = C0,
    Prob_neg = CN,
    Diagonal = V_diag,
    Population_seed = params$Semilla[1],
    Interacs_seed = params$Semilla[2],
    Growth_seed = params$Semilla[3]
  )
  
  # Revisar si un archivo con ese ID existe
  S_path <- paste("./test/Parameters/Seeds_save", ".tsv", sep = "") #Parameters
  exist <- file.exists(S_path) # flag
  
  if (!exist) { # File doesnt exist
    
    file.create(S_path) # Make file
    write.table(seeds_df, file = S_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
  } else {
    
    Og_table <- read.delim(S_path, sep = "\t", header = TRUE) # Read table
    Join_tables <- rbind(Og_table, seeds_df) # Join tables
    write.table(Join_tables, file = S_path, sep = "\t", row.names = FALSE, col.names = TRUE) # Save
  }
  return(ID)
}


#--------------------------------------------------------Function to generate data--------------------------------------

generate <- function(N,seeds,C0,CN, Val_diag){
  
  #------------------Populations-----------------------------#
  S_p <- sample(seeds, 1)
  Pobl <- vector("numeric", length = N)
  set.seed(S_p)
  for (i in 1:N) { # Generate initial populations
    Pobl[i] <- runif(1, min = 0.1, max = 1)
  }
  
  #--------------------Interactions-------------------------#
  
  S_i <- sample(seeds, 1) # Set seed
  set.seed(S_i)
  
  vec <- numeric(N * N)
  counter <- 0
  
  # Vectorized computation
  while (counter < N * N) {
    P_neg <- rbinom(n = 1, size = 1, prob = CN) 
    tmp <- rbinom(n = 1, size = 1, prob = 1 - C0) * ifelse(P_neg != 0, runif(1, min=0, max=1), -runif(1, min=0, max=1))
    
    # Store result directly at the appropriate index
    vec[counter + 1] <- tmp
    counter <- counter + 1
  }
  
  inter <- matrix(vec, nrow = N, ncol = N) # Hacerlo matriz
  # diag(inter) <- rnorm(N, mean = 0, sd = 1) # DIAGONAL
  diag(inter) <- Val_diag
  
  #------------------------Grow rates---------------------#
  S_g <- sample(seeds, 1)
  Grow <- vector("numeric", length = N)
  set.seed(S_g)
  for (i in 1:N){
    Grow[i] <- runif(1, min = 0.001, max = 1)
  }
  
  seed = c(S_p, S_i, S_g)
  
  return(list(Interacs = inter,
              Growth = Grow,
              Population = Pobl,
              Seeds = seed)
  )
}



#------------------------------------------------Steady state----------------------------------------------------------

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

find_SS <- function (ID, tol, params) {
  
  #------------------------Read output------------------------#
  # Load required package
  library(data.table)
  
  # Read table
  out_path <- paste("./test/Outputs/O_", ID , ".tsv", sep = "") # output
  out <- as.matrix( fread(out_path, sep = "\t") )
  times <- ncol(out)
  
  #----------------------Steady State search-----------------#
  result1 <- diff_SS(out, tol)
  result2 <- Rwindow_ALL(out,tol)
  
  
  #---------------------------Create DATA FRAME-------------#
  SS_df <- data.frame(
    'ID' = ID,
    'Population_seed' = params$Semilla[1],
    'Interaction_seed' = params$Semilla[2],
    'Growth_seed' = params$Semilla[3],
    'Generations' = times,
    'Diff_ALL' = ifelse(is.na(as.numeric(result1$Stable)), "Not found", as.numeric(result1$Stable)),
    'Rwindow_ALL' = ifelse(is.na(as.numeric(result2$Stable)), "Not found", as.numeric(result2$Stable)), 
    'Tolerance' = tol,
    'Individual' = FALSE
  )
  
  
  #-------------------------Save Data Frame---------------#
  SS_path <- paste("./test/Scan/SS_all", ".tsv", sep = "") #Parameters
  exist <- file.exists(SS_path) # Bandera
  
  if (!exist) { # File doesnt exist
    
    file.create(SS_path) # Make file
    write.table(SS_df, file = SS_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
  } else {
    
    SS_table <- read.delim(SS_path, sep = "\t", header = TRUE) # Read table
    Join_ss <- rbind(SS_table, SS_df) # Join tables
    write.table(Join_ss, file = SS_path, sep = "\t", row.names = FALSE, col.names = TRUE) # Save
  }
  
  return(list(Window_SS = result2$Stable, 
              Diff_SS = result1$Stable)
         )
}

#----------------------------------------------Simulation function-------------------------------------------------
test <- function (N, C0, CN, times, tol, individual, V_diag, n_sim) {
  
  #--------------------------Load all seeds-------------------------------#
  library(data.table)
  seeds <- as.matrix(fread("./Seeds.tsv", sep = "\t"))
  
  #-----------------------------Generate data----------------------------------#
  res <- generate(N, seeds, C0, CN, V_diag)
  
  params <- list(
    alpha = res$Interacs, # Interactions
    r = res$Growth, # Growth rates
    Pobl = res$Population, # Population
    Semilla = res$Seeds # Seeds
  )
  
  cat("Datos generados...\n")
  cat("El numero de semillas utilizadas son:", params$Semilla, "\n")
  cat("Population seed", params$Semilla[1], "\n")
  cat("Growth rates", params$r, "\n")
  cat("Interaction specie1-specie3--->", params$alpha[1, 3], "\n")
  
  #----------------------------Simulate ------------------------------------#
  glvmodel <- miaSim::simulateGLV(
    n_species = N, A = params$alpha, x0 = params$Pobl, 
    growth_rates = params$r, t_start = 0, t_store = times, 
    t_end = times, migration_p = 0, stochastic = FALSE, norm = TRUE
  )
  
  output <- glvmodel@assays@data@listData[["counts"]]
  cat("Primer renglon columnas 1:4:", output[1, 1:4], "\n")
  
  #-------------------------- Save results---------------------------------#
  ID <- save_test(N, C0, CN, V_diag, output, params)
  cat("El numero de ID es:", ID, "\n")
  
  
  #-------------------------find steady state------------------------------#
  tmp <- find_SS(ID, tol, params)
  cat("Simulation number 1 DONE \n\n")
  
  # Repeat loop for new populations
  for (sim in 2:n_sim) {
    
    S_p <- sample(seeds, 1)  # New population seed
    Pobl <- runif(N, min = 0.1, max = 1)  # Generate new populations
    
    cat("New population seed", S_p, "\n")
    cat("Growth rates", params$r, "\n")
    cat("Interaction specie1-specie3--->", params$alpha[1, 3], "\n")
    
    # Update only the population seed and populations
    params$Semilla[1] <- S_p
    params$Pobl <- Pobl
    
    # Simulate again with updated populations
    glvmodel <- miaSim::simulateGLV(
      n_species = N, A = params$alpha, x0 = params$Pobl, 
      growth_rates = params$r, t_start = 0, t_store = times, 
      t_end = times, migration_p = 0, stochastic = FALSE, norm = TRUE
    )
    
    output <- glvmodel@assays@data@listData[["counts"]]
    cat("Primer renglon columnas 1:4:", output[1, 1:4], "\n")
    
    # Save results and find steady state
    ID <- save_test(N, C0, CN, V_diag, output, params)
    cat("El numero de ID es:", ID, "\n")
    tmp <- find_SS(ID, tol, params)
    
    # Print simulation number
    cat("Simulation number", sim, "DONE\n\n")
  }
}


#--------------------------------------------------------Testing-----------------------------------------------------
N <- 2 # Number of species
C0 <- 0.45 # Prob. interaction =0
CN <- 0.2 # Prob. interaction <0
times <- 50 # Generations
tol <- 0.05 # Tolerance
V_diag <- -0.5
n_sim = 30 # Number of Simulations

# setwd("~/Documents/LAB_ECO")
wd <- "~/Cluster/"
setwd(wd)
test (N, C0, CN, times, tol, individual, V_diag, n_sim)


# Primera semilla = Poblaciones
# Segunda semilla = Interactions
# Tercera semilla = Growth rates



