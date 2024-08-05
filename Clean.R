# This Rscript is made for the clean version of GLV

#----------------------------------------------------Functions----------------------------------------------------------

# New function to save results
save = function(output,params, Pobl, Semilla){
  
  library(ids)
  
  # Generar el ID
  ID <- ids::random_id(1, 3)
  
  # Especificar el path del output
  out <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output
  pms <- paste("./Parameters/P_", ID , ".tsv", sep = "") #Parameters
  
  # Revisar si un archivo con ese ID existe
  exist <- file.exists(out) # Bandera
  while (exist==TRUE) { # False-> NO EXISTE
    
    # Generar el ID nuevo
    ID <- ids::random_id(1, 3)
    pms_a <- paste("./Parameters/P_", ID , ".tsv", sep = "") #Parameters
    exist <- file.exists(pms) # Revise si existe el archivo
    
  }
  
  # Generar los archivos
  system(paste("touch", out))
  system(paste("touch", pms))
  
  # Guardar la tabla OUTPUT
  write.table(output, file = out, sep = "\t", row.names = FALSE, col.names = TRUE)
  
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
  
  #--------------------------------------Save interactions-------------------------------------------------------
  cat("Interactions", file = pms)
  cat("\n", file = pms, append = TRUE)
  write.table(alp, file = pms, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  
  #--------------------------------------Save growth rates-------------------------------------------------------
  cat("Grow rates", file = pms, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  write.table(gr, file = pms, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  
  #--------------------------------------Save Initial Populations-------------------------------------------------------
  cat("Initial populations", file = pms, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  write.table(Pobl, file = pms, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  
  #--------------------------------------Save seeds-------------------------------------------------------
  cat("Seeds", file = pms, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  write.table(Semilla, file = pms, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  
  
  return(ID)
}

#Generar datos
generate <- function(N,seeds,C0,CN){
  
  #----------------------------------------------Populations------------------------------------------------------------
  # Generate initial populations
  S_p <- sample(seeds, 1)
  Pobl <- vector("numeric", length = N)
  set.seed(S_p)
  for (i in 1:N) {
    Pobl[i] <- runif(1, min = 0.1, max = 1)
  }
  
  #----------------------------------------------Interactions-----------------------------------------------------------
  # Set seed
  S_i <- sample(seeds, 1)
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
  
  #----------------------------------------------Grow rates-----------------------------------------------------------
  # Generate growth rates
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

#Leer parametros
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

# New function to find steady state


# New function to simulate extinctions
extinct <- function(ID, t_ext, S_ext, t_gens) {
  
  # Alert
  if (t_ext==0){
    print(paste("The time of extinction (", t_ext, ")", "hast to be greater than 0"))
    return()
  }
  
  #----------------------------------------------Read parameters--------------------------------------------------------
  R_params <- Read_params(ID)
  r <- unlist(R_params[2]) # Grows
  alpha <- as.matrix(R_params[[1]]) # Interactions
  N <- ncol(alpha)

  # Read output
  file <- paste("./Outputs/O_",ID,".tsv", sep = "")
  output <- read.csv(file, sep = "\t")
  
  # Alert 2 
  if (t_ext>ncol(output)){
    print(paste("The time of extinction (", t_ext, ")", "is not valid it is greater than the number of times simulated(", ncol(output) ,")"))
    return()
  }
  
  #----------------------------------------------Change parameters------------------------------------------------------
  # Columns-> times 
  # Rows -> species
  if (t_ext-1>0){ 
    Pobl <- output[,t_ext-1] # The extinction time is >1
  } else {
    Pobl <- output[,t_ext] # The extinction time is 1
    }
 
  # Asign 0 to the specie in the population
  Pobl[S_ext] <- 0
  
  #-----------------------------------------------------Simulate--------------------------------------------------------
  library(miaSim)
  glvmodel <- simulateGLV(n_species = N, 
                          A = alpha, # interaction matrix
                          x0 = Pobl, # Initial abundances
                          growth_rates = r, # Growth rates
                          t_start = 0, 
                          t_store = t_gens, 
                          t_end=t_gens, 
                          migration_p = 0,
                          stochastic = FALSE, # Ignorar ruido
                          norm = TRUE) # FALSE=conteo, TRUE=proporciones
  
  out <- glvmodel@assays@data@listData[["counts"]]
  
  #---------------------------------------Save new simulation-----------------------------------------------------------
  f_out <- paste("./Extinction/E_", ID ,"_S0", S_ext, "-T0", t_ext, ".tsv", sep = "") # output
  
  # Revisar si un archivo con ese ID existe
  exist <- file.exists(f_out) # Bandera
  if (exist==TRUE) { 
    print(paste("A file with the same conditions already exist, exiting..."))
    return()
  }
  
  # Generar los archivos
  system(paste("touch", f_out))
  
  # Guardar la tabla OUTPUT
  write.table(out, file = f_out, sep = "\t", row.names = FALSE, col.names = TRUE)
}



#--------------------------------------------------------testing--------------------------------------------------------

# Generar las semillas posibles
setwd("~/Documents/LAB_ECO") # Set Working Directory
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
library(miaViz)
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
miaViz::plotSeries(glvmodel, "time")

####
# Compare
negative_count <- sum(interacs == 0)
count_intrf = sum(interacs > 0)
count_test = sum(interacs < 0)


              #----------------------------------------Save-----------------------------------#
ID <- save(out,params,Pobl,Semilla)
  
              #----------------------------------------Extinction-----------------------------#
t_ext=3
S_ext=5
t_gens=25
extinct(ID, t_ext, S_ext, t_gens)
