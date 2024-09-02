
#----------------------------------------------------First save function--------------------------------------------
All_Save = function(output,params){
  
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
  Pobl <- params$Pobl
  Semilla <- params$Semilla
  
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




#----------------------------------------------------Second function save--------------------------------------------
Ms_save = function(N_specs, C0, CN, V_diag, output,params){
  
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


#----------------Testing-------------------#

setwd("~/Documents/LAB_ECO")
ID = Ms_save(N, C0, CN, V_diag, output, params)
All_Save(output,params)

# Test if saving the seeds would return the same values

test <- vector("numeric", length = N)
set.seed(params$Semilla[1])
for (i in 1:N) { # Generate initial populations
  test[i] <- runif(1, min = 0.1, max = 1)
}
# Comparar si son iguales
son_iguales <- all(params$Pobl == R_pobl)
print(son_iguales)





