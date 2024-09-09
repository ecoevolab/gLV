#-------------------------------------Read Parameters Lines method---------------------------------------
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

#-------------------------------------Read Parameters Seeds method---------------------------------------
Read_seeds <- function(ID){
  
  S_path <- paste("./Parameters/Seeds_save", ".tsv", sep = "") #Parameters
  tmp_table <- read.delim(S_path, sep = "\t", header = TRUE) # Read table
  
  pattern <- paste0("\\b", ID, "\\b") # regular expressions
  Selec_row <- tmp_table[grepl(pattern, tmp_table$ID_simulation), ]
  
  # Get initial values
  N_specs <- Selec_row$N_specs
  C0_t <- Selec_row$Prob_0
  CN_t <- Selec_row$Prob_neg
  
  #------------------Populations-----------------------------#
  R_pobl <- vector("numeric", length = N_specs)
  set.seed(Selec_row$Population_seed)
  for (i in 1:N_specs) { # Generate initial populations
    R_pobl[i] <- runif(1, min = 0.1, max = 1)
  }
  
  #--------------------Interactions-------------------------#
  
  vec <- numeric(N_specs * N_specs) # Empty matrix
  counter <- 0
  
  set.seed(Selec_row$Interacs_seed)
  while (counter < N_specs * N_specs) {
    P_neg <- rbinom(n = 1, size = 1, prob = CN_t) 
    tmp <- rbinom(n = 1, size = 1, prob = 1 - C0_t) * ifelse(P_neg != 0, runif(1, min=0, max=1), -runif(1, min=0, max=1))
    
    # Store result directly at the appropriate index
    vec[counter + 1] <- tmp
    counter <- counter + 1
  }
  
  inter <- matrix(vec, nrow = N, ncol = N) # Hacerlo matriz
  # diag(inter) <- rnorm(N, mean = 0, sd = 1) # DIAGONAL
  diag(inter) <- Selec_row$Diagonal
  
  # # Comparar si son iguales
  # son_iguales <- all(params$alpha == inter)
  # print(son_iguales)
  
  #------------------------Grow rates---------------------#
  Grow <- vector("numeric", length = N_specs)
  set.seed(Selec_row$Growth_seed)
  for (i in 1:N_specs){
    Grow[i] <- runif(1, min = 0.001, max = 1)
  }

  return(list(R_interacs = inter,
              R_Growth = Grow,
              R_Population = R_pobl)
         
  )
}

#-------------------------------------Testing----------------------------------------------------------------
# Run all the previous functions
Rdr_params <- Read_seeds(ID)


# Comparar si son iguales
son_iguales <- all(params$Pobl == Rdr_params$R_Population)
print(son_iguales)

