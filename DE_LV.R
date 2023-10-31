# Diferential equations solvers
# Limpiar memoria
rm(list=ls()) # Clear the memory

#---------------Definir ecuacion------------------------------------------------

GLV = function(t,state,params){

  with(as.list(c(state, params)),{ # Access the names of the variables in the function Pobl
    
    # Extraer numero de especies
    n <- length(state)
    
    # Vector para guardar cambio de especies
    dN <- numeric(length = n)
    
    # Calculate the rate of change for each species
    for (i in 1:n) {
      result <- (params$r[i] * state[i]) + (state[i] * sum(params$alpha[i, ] * state))
      dN[i] <- result
    }
    return(list(dN)) # Return
  })
}

#-----------------------------Guardar resultados--------------------------------

save = function(output,params, Pobl, Semilla){
  
  library(ids)
  
  # Generar el ID
  ID <- ids::random_id(1, 3)
  
  # Set Working Directory
  setwd("~/Documents/LAB_ECO")
  
  # Especificar el path del output
  out <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output
  pms <- paste("./Parameters/A_", ID , ".tsv", sep = "") #Parameters
  
  # Revisar si un archivo con ese ID existe
  exist <- file.exists(pms) # Bandera
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
  
  # Guardar la INTERACCIONES
  cat("Interactions", file = pms)
  cat("\n", file = pms, append = TRUE)
  write.table(alp, file = pms, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  
  # Guardar los CRECIMIENTOS
  cat("Grow rates", file = pms, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  write.table(gr, file = pms, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  
  # Guardar POBLACIONES INICIALES
  cat("Initial populations", file = pms, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  write.table(Pobl, file = pms, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  
  # Guardar Semillas
  cat("Seeds", file = pms, append = TRUE)
  cat("\n", file = pms, append = TRUE)
  write.table(Semilla, file = pms, sep = "\t", col.names = NA, append = TRUE)
  cat("\n", file = pms, append = TRUE)
 
  
  return(ID)
}




#------------------------------Generar datos------------------------------------
generate <- function(N,seeds){
  
  # Sample from the vector 'a' 1 element.
  S_p <- sample(seeds, 1)
  S_i <- sample(seeds, 1)
  S_g <- sample(seeds, 1)
  
  # Vector de especies POBLACIONES INCIALES
  Pobl <- vector("numeric", length = N)
  set.seed(S_p)
  for (i in 1:N) {
    Pobl[i] <- round(abs(rnorm(1, 10, 5)))
  }
  
  # Hacer tabla de interacciones
  inter <- matrix(nrow = N, ncol = N)
  set.seed(S_i)
  for (r in 1:N) { # column
    for (c in 1:N){ # row
      
      # Diagonal
      if (r==c) { 
        #inter [r,c] <- -runif(1)
        inter [r,c] <- rnorm(1, mean = 0, sd = 1)
      } else {
        inter [r,c] <- rnorm(1, mean = 0, sd = 1)
      }
    }
  }
  
  # Vector de crecimiento de especies
  Grow <- vector("numeric", length = N)
  set.seed(S_g)
  for (i in 1:N){
    Grow[i] <- sample.int(n=100, size=1)
  }
  
  seed = c(S_p, S_i, S_g)
  
  result <- list(inter, Grow, Pobl, seed)
  return(result)
}


#----------------------Population graphs-------------------------------------------------

Pgraph <- function(output, specs){
  
  library(ggplot2)
  library(reshape2)
  
  # Graficar resultados
  data <- as.data.frame(output)
  data <- data[specs]
  data_merge <- melt(data, id.vars = "time", variable.name = "Population", value.name = "Value")
  
  # Load the reshape2 package
  library(reshape2)
  
  # Number
  ggplot(data_merge, aes(x = time , y =  Value, color = Population)) +   # Define data and aesthetics
    geom_point() +      # Add a scatter plot layer
    scale_x_continuous() +  # Explicitly specify the x-axis scale
    scale_y_continuous() +  # Explicitly specify the y-axis scale
    labs(title = "Population change", x = "Time", y = "Population")
}

#----------------------------------Graficar redes nodos-------------------------------
red <- function(params){
  
  # Graficar redes
  library(igraph)
  
  # Sample data frame with nodes (columns)
  data <- params$alpha
  
  
  # Convert the data frame to an adjacency matrix
  graph <- graph_from_adjacency_matrix(data, mode = "directed", weighted = TRUE)
  E(graph)$weight <- abs(E(graph)$weight) #make weight +
  
  #Method 2
  #plot(graph, layout = layout_with_fr(graph), edge.label = E(graph)$weight, edge.label.color = "black")
  
  # Opcion 2
  # Set the minimum and maximum edge widths
  min_width <- 1  # Minimum edge width
  max_width <- 5  # Maximum edge width
  
  # Normalize edge weights between min_width and max_width
  edge_weights <- E(graph)$weight
  normalized_weights <- (edge_weights - min(edge_weights)) / (max(edge_weights) - min(edge_weights))
  edge_widths <- min_width + (max_width - min_width) * normalized_weights
  
  
  # Plot the graph with customized edge widths
  plot(graph, edge.width = edge_widths)
  
}

#----------------------------Calcular negativos output--------------------------

scan_neg <- function(output, params, ID){
  
  # Contar AL MENOS UN NEGATIVO
  start_col <- 2  # Primera columna es la segunda, ya que primer es tiempo
  end_col <- ncol(output)  # Ultima columna
  has_negative <- apply(output[,start_col:end_col], 2, function(col) any(col < 0))
  least_1neg <- sum(has_negative) #Numero de columnas con al menos un negativo
  
  # Contar cuantas especies TERMINAN EN NEGATIVO
  last_row <- output[nrow(output), ]
  end_n <- sum(last_row < 0)
  
  # Contar numero de INTERACCIONES NEGATIVAS
  neg_int <- sum(params$alpha < 0)
  
  # Contar numero de INTERACCIONES POSITIVAS
  pos_int <- sum(params$alpha > 0)
  
  # Guardar SCANEO 
  setwd("~/Documents/LAB_ECO") # Set Working Directory
  scn1 <- "./Scan/Scan_res.tsv" #Scan path
 
  if (file.exists(scn1)!=TRUE) { # False-> NO EXISTE
     system(paste("touch", scn1))
  }
  
  # Generar data frame
  Vec_scan <- data.frame(ID, least_1neg, end_n, neg_int, pos_int)
  
  #Si el archivo esta vacio añadir nombres de columnas
  if (file.info(scn)$size == 0) {
    write.table(Vec_scan, file = scn, sep = "\t", col.names = TRUE, row.names = FALSE, append = TRUE)
    
  } else {
      write.table(Vec_scan, file = scn, sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
  }
}

#------------------------Plot interactions-negative end-------------------------

Scan_plot <- function() {
  
  setwd("~/Documents/LAB_ECO") # Set Working Directory
  scn <- "./Scan/Scan_res.tsv"
  data <- read.table(scn, header = TRUE, sep = "\t")
  
  library(ggplot2) 
  
  ggplot(data, aes(x=neg_int, y=end_n)) +
    geom_point() +                 # Add points
    labs(x = "# Negative interactions",  y = "# species end negative")   
  
}

#---------------------------Generar simulacion----------------------------------
# rm(list=setdiff(ls(), c("GLV", "save", "generate", "Pgraph", "red")))

# Generar las semillas posibles
library (random)
seeds <- randomNumbers(n=100, min=1, max=200)

# Introducir el numero de especies
N <- 3

# Generar datos y extraerlos
res <- generate(N,seeds)
V_inter <- unlist(res[1])
params <- list(
  r = unlist(res[2]), # Grow rates
  alpha = matrix(V_inter, nrow = N, ncol = N) # Interaction
)
Pobl <- unlist(res[3])
Semilla <- unlist(res[4])


# Evaluar la euacion
library(deSolve)
times <- seq(0, 20, by = 1)
output <- ode(y = Pobl, times = times, func = GLV, parms = params,  method = "lsoda", atol = 1e+1, rtol = 1e+1,
              maxsteps = 100)

# Guardar resultados
ID <- save(output, params, Pobl, Semilla)

# Contar negativos
scan_neg(output, params, ID)

# Plot interacciones negativas (x) especies que terminan en netivas (y)
Scan_plot ()   

# Grafica poblacion
esp <- c("time", "1","3")
Pgraph(output, esp)

# Graficar redes
red(params)

#-----------------------------EJEMPLO ARTICULO----------------------------------


# Load library
library(deSolve)

# Set time points for simulation
times <- seq(0, 20, by = 1)

# Set the number of species
n <- 4

# Solve the differential equations
Pobl <- c(1.2, 0.3 , 2 , 0.001)

params <- list(
  r = c(0.044, 0.216 , 0.116 , 0.2),  # Growth rates for each species
  alpha = matrix(c(-0.08, 0.02, 0.08, 0,
                   -0.04, -0.08, 0.04, 0,
                   -0.16, 0.16, -0.08, 0,
                   0, 0, 0, -0.1
                  ), nrow = n, byrow = TRUE)  # Interaction coefficients
)

output <- ode(y = Pobl, times = times, func = GLV, parms = params,  method = "lsoda", atol = 1e+1, rtol = 1e+1,
              maxsteps = 100)

Semilla <- c(12,52,3)

# Grafica poblacion
esp <- c("time", "1","3")
Pgraph(output, esp)

# Graficar redes
red(params)

# Guardar resultados ANTES QUE SCAN
save(output, params, Pobl, Semilla)

# Contar negativos
scan <- scan_neg(output, params, ID)

#------------------------------Limpiar carpetas---------------------------------

# Outputs
unlink("~/Documents/LAB_ECO/Outputs/*", recursive = TRUE, force = TRUE)

# Parameters
unlink("~/Documents/LAB_ECO/Parameters/*", recursive = TRUE, force = TRUE)

# Scan
unlink("~/Documents/LAB_ECO/Scan/*", recursive = TRUE, force = TRUE)

#--------------------------------Notas----------------------------------------------

"
1. Cuando tengo una seed de 20 hay un error, en esta matriz hay negativos en la diagonal de las interacciones.
2. Las poblaciones llegan a alcanzar negativos, esto no debe ser posible
3. No me permite tener una tolerancia maxima menor a e+1, lo mismo sucede con rtol
4. El numero de maxsteps depende mucho de las poblaciones iniciales
"

"
NUEVAS:
Para las alphas una distribucion uniforme runif de 0 a 1.
READR PARA TABLAS 

- 20 generaciones
- Especies constantes
- Por cada red que genere, calcular cuanto son positivos y cuantos son negativos. Tambien lo mismo con la matriz de interaciones

En cuantas generaciones se convierten en negativos.
Hacer una correlación entre el numero de simulaciones negativas

COMPLETADO
- Hacer una funcion para las graficas de las especies, tomando en cuenta las especies que queremos.
- Tambien graficar las redes.
- Guardar en 2 tablas (parametros y output)
- En la ultima generacion contar cuantas terminan negativo.
- Cuantas especies tuvieron al menos un negativo.
- Eliminar todo lo de las carpetas para probar codigo
- Plot con el numero de interacciones negativas (y) y el numero de especies que terminan en negativas (x)
- Para fuera de la diagonal en la matriz de interacciones una distribucion normal estandar.
- Para la diagonal usar una version normal y otra con negativos (runif)

"