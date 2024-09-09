
#-----------------------------------------------Seed saver function-----------------------------------------------------

Ms_save = function(N_specs, C0, CN, V_diag, output, params){
  
  library(ids)
  
  # Generar el ID
  ID <- ids::random_id(1, 3)
  
  # Especificar el path del output
  out_path <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output
  
  # --------------------------------------------Save OUTPUT---------------------------------------------------------#
  
  # Revisar si un archivo con ese ID existe
  exist <- file.exists(out_path) # flag
  
  while (exist) { # TRUE-> file exist, change ID...
    
    # Generar el ID nuevo
    ID <- ids::random_id(1, 3)
    out_path <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output
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
  S_path <- paste("./Parameters/Seeds_save", ".tsv", sep = "") #Parameters
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


#-----------------------------------------------Data generation function------------------------------------------------

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


#-----------------------------------------------Steady State function---------------------------------------------------

SS_all <- function(ID, tol, wd) {
  
  #-------------------------Read table----------------------------#
  # Load required package
  library(data.table)
  
  # Read table
  out_path <- paste(wd, "Outputs/O_", ID , ".tsv", sep = "") # output
  out <- as.matrix( fread(out_path, sep = "\t") )
  
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
  
  return(list(Method_dif = first_col,
              Dif_means = Stb_mean)
  )
}


#-------------------------------------Code profiling function-----------------------------------------------------------
CPr_sim <- function(N, C0, CN, times, tol, individual, V_diag, wd) { 
  
  #--------------------------Load all seeds-------------------------------#
  library(data.table)
  seeds_path <- paste(wd, "./Seeds.tsv", sep="")
  seeds <- as.matrix( fread(seeds_path, sep = "\t") ) # Read table
  
  #-----------------------------Generate data----------------------------------#
  CPr_generate <- system.time({
    res <- generate(N, seeds, C0, CN, V_diag) # Results
    
    params <- list(
      alpha = res$Interacs, # Interactions
      r = res$Growth, # Grow rates
      #alpha = matrix(unlist(res[[1]]) , nrow = N, ncol = N) 
      Pobl = res$Population, # Population
      Semilla = res$Seeds #Semillas
    )
  })
  cat("Datos generados...", "\n")
  
  # Primera semilla = Poblaciones
  # Segunda semilla = Interactions
  # Tercera semilla = Growth rates
  cat("El numero de semillas utilizadas son:", params$Semilla, "\n")
  
  #----------------------------------------Simulate----------------------------#
  CPr_glv <- system.time({
    glvmodel <- miaSim::simulateGLV(n_species = N, 
                                    A = params$alpha, # interaction matrix
                                    x0 = params$Pobl, # Initial abundances
                                    growth_rates = params$r, # Growth rates
                                    t_start = 0, 
                                    t_store = times, 
                                    t_end= times, 
                                    migration_p = 0,
                                    stochastic = FALSE, # Ignorar ruido
                                    norm = TRUE) # FALSE=conteo, TRUE=proporciones
    
    output <- glvmodel@assays@data@listData[["counts"]]
  })
  cat("Primer renglon colmnas 1:4:", output[1,1:4], "\n")
  
  #------------------Save Simulation results-----------------------------------#
  CPr_Ms_save <- system.time({
    ID <- Ms_save(N, C0, CN, V_diag, output, params)
  })
  cat("El numero de ID es:", ID, "\n")
  
  #---------------------Search for steady state-----------------------------------#
  CPr_SS <- system.time({
    result <- SS_all(ID, tol, wd)
  })
  cat("Steady state search DONE",  "\n")
  
  #--------------------Save Code profiling times---------------------------------#
  CPr_df <- data.frame(
    ID = ID ,
    Specs = N , 
    Generations = times,
    'Generation_time-s' = round(CPr_generate[3], 6) ,
    'Simulation_time-s' = round(CPr_glv[3], 6) ,
    'Ms_save_time-s' = round(CPr_Ms_save[3], 6),
    'SS_time-s' = round(CPr_SS[3], 6),
    Tolerance = tol,
    Individual = individual
  )
  
  library(readr)
  CPr_path <- paste(wd, "./Scan/CPr_time", ".tsv", sep = "") # output
  exist <- file.exists(CPr_path) # Flag
  
  if (!exist) { # File doesnt exist
    
    file.create(CPr_path) # Make file
    write.table(CPr_df, file = CPr_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
    cat("Code profiling times SAVED",  "\n") 
  } else {
    
    CPr_table <- read.delim(CPr_path, sep = "\t", header = TRUE) # Read table
    Join_CPr <- rbind(CPr_table, CPr_df) # Join tables
    write.table(Join_CPr, file = CPr_path, sep = "\t", row.names = FALSE, col.names = TRUE) # Save
    cat("Code profiling times SAVED",  "\n") 
  }
}
