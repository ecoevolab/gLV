#----------------------------------------------------Functions----------------------------------------------------------

# New function to save results
save = function(output,params, Pobl, Semilla){
  
  library(ids)
  
  # Generar el ID
  ID <- ids::random_id(1, 3)
  
  # Especificar el path del output
  out_path <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output
  pms_path <- paste("./Parameters/P_", ID , ".tsv", sep = "") #Parameters
  
  # Revisar si un archivo con ese ID existe
  exist <- file.exists(out_path) # Bandera
  while (exist) { # TRUE-> EXISTE
    
    # Generar el ID nuevo
    ID <- ids::random_id(1, 3)
    pms_path <- paste("./Parameters/P_", ID , ".tsv", sep = "") #Parameters
    exist <- file.exists(pms_path) # Revise si existe el archivo
    
  }
  
  # Create files
  file.create(out_path)
  file.create(pms_path)
  
  # Guardar la tabla OUTPUT
  write.table(output, file = out_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  # Separate parameters
  alp <- params$alpha
  gr <- params$r
  
  # Assign names
  name <- paste("S_", 1:nrow(alp))
  rownames(alp) <- name
  colnames(alp) <- name
  names(gr) <- name
  names(Pobl) <- name
  names(Semilla) <- c("Seed_pop", "Seed_inter", "Seed_grow")
  
  #--------------------------------------Save interactions--------------------------#
  cat("Interactions", file = pms_path)
  cat("\n", file = pms_path, append = TRUE)
  write.table(alp, file = pms_path, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms_path, append = TRUE)
  
  #--------------------------------------Save growth rates--------------------------#
  cat("Grow rates", file = pms_path, append = TRUE)
  cat("\n", file = pms_path, append = TRUE)
  write.table(gr, file = pms_path, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms_path, append = TRUE)
  
  #--------------------------------------Save Initial Populations-------------------#
  cat("Initial populations", file = pms_path, append = TRUE)
  cat("\n", file = pms_path, append = TRUE)
  write.table(Pobl, file = pms_path, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms_path, append = TRUE)
  
  #--------------------------------------Save seeds---------------------------------#
  cat("Seeds", file = pms_path, append = TRUE)
  cat("\n", file = pms_path, append = TRUE)
  write.table(Semilla, file = pms_path, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms_path, append = TRUE)
  
  
  return(ID)
}

#--------------------------------------------------------Function to generate data--------------------------------------#

generate <- function(N,seeds,C0,CN){
  
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
  diag(inter) <- -0.5
  
  #------------------------Grow rates---------------------#
  S_g <- sample(seeds, 1)
  Grow <- vector("numeric", length = N)
  set.seed(S_g)
  for (i in 1:N){
    Grow[i] <- runif(1, min = 0.001, max = 1)
  }
  
  seed = c(S_p, S_i, S_g)
  result <- list(Interacs = inter,
                 Growth = Grow,
                 Population = Pobl,
                 Seeds = seed)
  return(result)
}


#------------------------------------- New function to find steady state---------------------------------------------#


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
  first_col <- which(Stb_mat < log(tol), arr.ind = TRUE)[1] # Generation where the mean<tolerance
  
  tmp_df <- data.frame(
    ID = ID,
    Steady_generation = first_col
  )
  
  #---------------------------Save generation-------------------------#
  SS_path <- paste("./Scan/SS_all", ".tsv", sep = "") #Parameters
  exist <- file.exists(SS_path) # Bandera
  
  if (!exist) { # File doesnt exist
    
    file.create(SS_path) # Make file
    write.table(tmp_df, file = SS_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
  } else {
    
    SS_table <- read.delim(SS_path, sep = "\t", header = TRUE) # Read table
    Join_ss <- rbind(SS_table, tmp_df) # Join tables
    write.table(Join_ss, file = SS_path, sep = "\t", row.names = FALSE, col.names = TRUE) # Save
  }
  
  return(list(Stb_mat = Stb_mat, Stb_mean = Stb_mean, Fc = first_col))
}


#------------------------------------- New function for code profiling---------------------------------------------#
CPr_sim <- function(N, C0, CN, times,seeds_n, tol, indivdual) { 
  
  CPr_time1 <- system.time({
    seeds <- random::randomNumbers(n=seeds_n, min=1, max=100)
  })
  
  #-----------------------------Generate data----------------------------------#
  CPr_time2 <- system.time({
    res <- generate(N,seeds,C0, CN) # Results
    
    # V_inter <- unlist(res[[1]]) 
    params <- list(
      r = unlist(res[2]), # Grow rates
      alpha = matrix(unlist(res[[1]]) , nrow = N, ncol = N) # Interaction
    )
    Pobl <- unlist(res[3])
    Semilla <- unlist(res[4])
  })
  cat("Datos generados...", "\n")
  cat("El numero de semillas es de:", Semilla, "\n")
  
  #----------------------------------------Simulate----------------------------#
  CPr_time3 <- system.time({
    #library(miaViz)
    interacs <- params$alpha
    glvmodel <- miaSim::simulateGLV(n_species = N, 
                                    A = params$alpha, # interaction matrix
                                    x0 = Pobl, # Initial abundances
                                    growth_rates = params$r, # Growth rates
                                    t_start = 0, 
                                    t_store = times, 
                                    t_end=times, 
                                    migration_p = 0,
                                    stochastic = FALSE, # Ignorar ruido
                                    norm = TRUE) # FALSE=conteo, TRUE=proporciones
    
    output <- glvmodel@assays@data@listData[["counts"]]
  })
  cat("Primer renglon colmnas 1:4:", output[1,1:4], "\n")
  
  #------------------Save Simulation results-----------------------------------#
  CPr_time4 <- system.time({
    ID <- save(output,params,Pobl,Semilla)
  })
  cat("El numero de ID es:", ID, "\n")
  
  #---------------------Search for steady state-----------------------------------#
  CPr_time5 <- system.time({
    out_path <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output
    result <- SS_all(ID,out_path,tol)
  })
  cat("Steady state search DONE",  "\n")
  
  #--------------------Save Code profiling times---------------------------------#
  CPr_df <- data.frame(
    ID = ID ,
    Specs = N , 
    Generations = times,
    Seed_num = seeds_n,
    'Seeds_time-s' = round(CPr_time1[3], 6) , 
    'GenDat_time-s' = round(CPr_time2[3], 6) ,
    'Simulation_time-s' = round(CPr_time3[3], 6) ,
    'Save_time-s' = round(CPr_time4[3], 6),
    'SSS_time-s' = round(CPr_time5[3], 6),
    Tolerance = tol,
    Individual = individual
  )
  
  library(readr)
  CPr_path <- paste("./Scan/CPr_time", ".tsv", sep = "") # output
  exist <- file.exists(CPr_path) # Flag
  
  
  cat("Saving code profiling times...",  "\n") 
  if (!exist) { # File doesnt exist
    
    file.create(CPr_path) # Make file
    write.table(CPr_df, file = CPr_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
  } else {
    
    CPr_table <- read.delim(CPr_path, sep = "\t", header = TRUE) # Read table
    Join_CPr <- rbind(CPr_table, CPr_df) # Join tables
    write.table(Join_CPr, file = CPr_path, sep = "\t", row.names = FALSE, col.names = TRUE) # Save
  }
}

#--------------------------------------------------------Testing-----------------------------------------------------#
N <- 50 # Number of species
C0 <- 0.45 # Prob. interaction =0
CN <- 0.2 # Prob. interaction <0
times <- 3000 # Generations
seeds_n <- 300 # Number of possible seeds
tol <- 0.05 # Tolerance
individual <- FALSE

counter <- 1 
STOP <- FALSE # Flag
ct <- 1 # Number of simulation
setwd("~/Documents/LAB_ECO")

cat("\n" , "\n", strrep("#", 25))
current_date <- Sys.Date()
formatted_date <- format(current_date, "%d %B %Y")
cat("Simulation started at:", formatted_date , strrep("#", 25), "\n")

while (!STOP) {
  
  
  # Run the simulation function
  suppressWarnings(CPr_sim(N, C0, CN, times, seeds_n, tol, individual))
  
  cat("Simulation number", ct ,"DONE", "\n")
  cat("-----------------------------------------", "\n")
  
  if (counter==3) {  # Check if `times` has been updated 3 times
    times <- times + 100  # Add 100 generations
    counter <- counter + 1  # Increment `times` update counter
  } else {
    N <- N + 50  # Add 50 species
    counter <- 0  # Reset `times` update counter
  }
  
  # Stop simulation when there are at least 500 species and 10,000 generations
  if (N==1000) {
    STOP <- TRUE
  }
  
  ct <- ct +1
}