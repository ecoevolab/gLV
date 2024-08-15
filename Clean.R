# This Rscript is made for the clean version of GLV

#----------------------------------------------------Functions----------------------------------------------------------

# New function to save results
save = function(output,params, Pobl, Semilla){
  
  library(ids)
  
  # Generar el ID
  ID <- ids::random_id(1, 3)
  
  # Especificar el path del output
  out_path <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output
  pms <- paste("./Parameters/P_", ID , ".tsv", sep = "") #Parameters
  
  # Revisar si un archivo con ese ID existe
  exist <- file.exists(out_path) # Bandera
  while (exist==TRUE) { # False-> NO EXISTE
    
    # Generar el ID nuevo
    ID <- ids::random_id(1, 3)
    pms_a <- paste("./Parameters/P_", ID , ".tsv", sep = "") #Parameters
    exist <- file.exists(pms) # Revise si existe el archivo
    
  }
  
  # Generar los archivos
  system(paste("touch", out_path))
  system(paste("touch", pms))
  
  # Guardar la tabla OUTPUT
  write.table(output, file = out_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  # Separar los parametros
  alp <- params$alpha
  gr <- params$r
  
  # Asignar nombres 
  name <- paste("S_", 1:nrow(alp))
  rownames(alp) <- name
  colnames(alp) <- name
  names(gr) <- name
  names(Pobl) <- name
  names(Semilla) <- c("Seed_pop", "Seed_inter", "Seed_grow")
  
  #--------------------------------------Save interactions--------------------------#
  cat("Interactions", file = pms)
  cat("\n", file = pms, append = TRUE)
  write.table(alp, file = pms, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  
  #--------------------------------------Save growth rates--------------------------#
  cat("Grow rates", file = pms, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  write.table(gr, file = pms, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  
  #--------------------------------------Save Initial Populations-------------------#
  cat("Initial populations", file = pms, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  write.table(Pobl, file = pms, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  
  #--------------------------------------Save seeds---------------------------------#
  cat("Seeds", file = pms, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  write.table(Semilla, file = pms, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  
  
  return(ID)
}

#Generar datos
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
  result <- list(inter, Grow, Pobl, seed)
  return(result)
}

#---------------------------------------------------Read parameters------------------------------------------------#
Read_params <- function(ID){
  
  library(readr)
  
  # Generate path
  path <- paste("Parameters/P_",ID,".tsv",sep = "")
  
  # Leer tsv como un vector de caracteres
  tsv_lines <- readLines(path)
  
  # Obtener lineas donde comienzan mis datos
  Interactions_line <- grep("Interactions", tsv_lines)
  Grows_line <- grep("Grow rates", tsv_lines)
  Population_line <- grep("Initial population", tsv_lines)
  Seed_line <- grep("Seeds", tsv_lines)
  last_line <- grep("Seed_grow", tsv_lines)
  
  # Obtener tabla de interacciones
  tmp <- tsv_lines[(Interactions_line + 1):Grows_line-1]
  Interacs <- read.table(text = tmp, sep = "\t", skip=1, header=TRUE, row.names = 1)
  
  # Obtener tabla de grows rates
  tmp <- tsv_lines[(Grows_line + 1):Population_line-1]
  Grows <- read.table(text = tmp, sep = "\t", skip=1, header=TRUE, row.names = 1)
  
  # Obtener tablas de poblaciones iniciales
  tmp <- tsv_lines[(Population_line + 1):Seed_line-1]
  Pobl <- read.table(text = tmp, sep = "\t", skip=1, header=TRUE, row.names = 1)
  
  # Obtener tablas de semillas
  tmp <- tsv_lines[(Seed_line):last_line]
  Seed <- read.table(text = tmp, sep = "\t", skip=1, header=TRUE, row.names = 1)
  
  result <- list(Interacs, Grows, Pobl, Seed)
  return(result)
  
}

#------------------------------------- New function to find steady state---------------------------------------------#
st_search = function(ID, out, tol, individual) {
  
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
  
  if (individual) {
    
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
    
  } else {
    
    Stb_mean <- colMeans(ln_mat) # Obtener promedios 
    condition <- Stb_mean < log(tol)    # Define the condition 
    first_col <- which(condition)[1] # Find the first column where the condition is met
    return(list(Stb_mat = Stb_mat, Stb_mean = Stb_mean, Fc = first_col ))
  }
}


#------------------------------------- New function for code profiling---------------------------------------------#
CPr_sim <- function(N, C0, CN, times,seeds_n, tol, individual) { 
  
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
  
  #------------------Save Simulation results-----------------------------------#
  CPr_time4 <- system.time({
    ID <- save(output,params,Pobl,Semilla)
  })
  
  #---------------------Search for steady state-----------------------------------#
  CPr_time5 <- system.time({
    out_path <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output
    result <- st_search(ID,out_path,tol,individual)
  })
  
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
  
  if (!exist) { # File doesnt exist
    
    system(paste("touch", CPr_path)) # Make file
    write.table(CPr_df, file = CPr_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
  } else {
    
    CPr_table <- read.delim(CPr_path, sep = "\t", header = TRUE) # Read table
    Join_CPr <- rbind(CPr_table, CPr_df) # Join tables
    write.table(Join_CPr, file = CPr_path, sep = "\t", row.names = FALSE, col.names = TRUE) # Save
  }
  
  
  
}



    #----------------------------- New function to simulate extinctions---------------------------#

  
  


#--------------------------------------------------------Testing--------------------------------------------------------

N <- 50 # Number of species
C0 <- 0.45 # Prob. interaction =0
CN <- 0.2 # Prob. interaction <0
times <- 3000 # Generations
seeds_n <- 300
tol <- 5
individual  <- TRUE
setwd("~/Documents/LAB_ECO") # Set Working Directory

CPr_sim(N,C0,CN,times,seeds_n,tol,individual)



unlink("~/Documents/LAB_ECO/Scan/CPr_time.tsv", recursive = TRUE, force = TRUE) # Remove all data from Scan 




#########################################################################################################################
# Generar las semillas posibles
setwd("~/Documents/LAB_ECO") # Set Working Directory

for (simulations in 1:15) {
  
  library (random)
  seeds <- randomNumbers(n=300, min=1, max=100)
  
  N <- 20 # Number of species
  C0 <- 0.45 # Prob. interaction =0
  CN <- 0.2 # Prob. interaction <0
  times <- 3000 # Generations
  
  # Generar datos y extraerlos
  res <- generate(N,seeds,C0, CN) # Results
  
  # V_inter <- unlist(res[[1]]) 
  params <- list(
    r = unlist(res[2]), # Grow rates
    alpha = matrix(unlist(res[[1]]) , nrow = N, ncol = N) # Interaction
  )
  Pobl <- unlist(res[3])
  Semilla <- unlist(res[4])
  
  #----------------------------------------Simulate----------------------------#
  library(miaSim)
  #library(miaViz)
  interacs <- params$alpha
  glvmodel <- simulateGLV(n_species = N, 
                          A = params$alpha, # interaction matrix
                          x0 = Pobl, # Initial abundances
                          growth_rates = params$r, # Growth rates
                          t_start = 0, 
                          t_store = times, 
                          t_end=times, 
                          migration_p = 0,
                          stochastic = FALSE, # Ignorar ruido
                          norm = TRUE) # FALSE=conteo, TRUE=proporciones
  
  out <- glvmodel@assays@data@listData[["counts"]]
  
  #------------------Save Simulation results-----------------------------------#
  ID <- save(out,params,Pobl,Semilla)
  
  #------------------Search for ALL steady state--------------------------------#
  individual  <- FALSE
  tol <- 5
  
  result <- st_search(ID,out_path,tol,individual) # Search for steady state
  Stb_mean <- data.frame(result$Stb_mean) # Mean of deltas 
  st_all_save(ID, tol, Stb_mean)
  St_time <- result$Fc # Generation where Mean_Delta < log(tolerance)
  
  if (simulations==1) {
    St_time_vec <- c()
  } else {
    St_time_vec <- c(St_time_vec,St_time)
  }
 
}


#########################################################################################################################
out <- glvmodel@assays@data@listData[["counts"]]
miaViz::plotSeries(glvmodel, "time")

####
# Compare
negative_count <- sum(interacs == 0)
count_intrf = sum(interacs > 0)
count_test = sum(interacs < 0)

  
            
