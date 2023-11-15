# Diferential equations solvers
# Limpiar memoria
rm(list=ls()) # Clear the memory

#---------------Definir ecuaciones----

GLV = function(t,state,params){

  with(as.list(c(state, params)),{ # Access the names of the variables in the function Pobl
    
    # Extraer numero de especies
    n <- length(state)
    
    # Vector para guardar cambio de especies
    dN <- numeric(length = n)
    
    # Calculate the rate of change for each species
    for (i in 1:n) {
      result <- (params$r[i] * state[i]) + (state[i] * sum(params$alpha[i, ] * state))
      #result <- (params$r * state) + (state * sum(params$alpha %*% state))
      dN[i] <- result
    }
    return(list(dN)) # Return
  })
}

#Guardar resultados
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

#Generar datos
generate <- function(N,seeds,C){
  
  # Sample from the vector 'a' 1 element.
  S_p <- sample(seeds, 1)
  S_i <- sample(seeds, 1)
  S_g <- sample(seeds, 1)
  
  # Vector de especies POBLACIONES INCIALES
  Pobl <- vector("numeric", length = N)
  set.seed(S_p)
  for (i in 1:N) {
    Pobl[i] <- runif(1, min = 0.1, max = 1)
  }
  
  # Hacer tabla de INTERACCIONES
  set.seed(S_i)
  vec <- rbinom(n=N*N,size=1,p=1-C)*runif(n=N*N)
  inter <- matrix(vec, nrow = N, ncol = N) 
  diag(inter) <- rnorm(N, mean = 0, sd = 1)
  
  # Vector de crecimiento de especies
  Grow <- vector("numeric", length = N)
  set.seed(S_g)
  for (i in 1:N){
    Grow[i] <- runif(1, min = 0.001, max = 1)
  }
  
  seed = c(S_p, S_i, S_g)
  result <- list(inter, Grow, Pobl, seed)
  return(result)
}


#Population graphs
Pobl_plot <- function(output, specs){
  
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

#Graficar redes nodos
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

# Calcular negativos output
Scan_neg <- function(output, params, ID, N, C){
  
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
  scn1 <- "./Scan/Scan_res.tsv" #Scan path
 
  if (file.exists(scn1)!=TRUE) { # False-> NO EXISTE
     system(paste("touch", scn1))
  }
  
  # Guardar NUMERO ESPECIES
  N_species <- N
  
  # Guardar PROBABILIDAD DE 0
  Prob_0 <- C
  
  # Guardar INTERACCIONES=0
  interac_0 <- sum(params$alpha==0)
    
  # Generar data frame
  Vec_scan <- data.frame(ID, least_1neg, end_n, neg_int, pos_int, N_species, Prob_0, interac_0)
  
  #Si el archivo esta vacio añadir nombres de columnas
  if (file.info(scn1)$size == 0) {
    write.table(Vec_scan, file = scn1, sep = "\t", col.names = TRUE, row.names = FALSE, append = TRUE)
    
  } else {
      write.table(Vec_scan, file = scn1, sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
  }
}

# Plot interactions-negative end
Scan_plot <- function() {
  
  scn <- "./Scan/Scan_res.tsv"
  data <- read.table(scn, header = TRUE, sep = "\t")
  
  library(ggplot2) 
  ggplot(data, aes(x=neg_int, y=end_n)) +
    geom_point() +                 # Add points
    labs(x = "# Negative interactions",  y = "# species end negative")   
  
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
  
  # Obtener tabla de interacciones
  tmp <- tsv_lines[(Interactions_line + 1):Grows_line-1]
  Interacs <- read.table(text = tmp, sep = "\t", skip=1, header=TRUE, row.names = 1)
  
  # Obtener tabla de grows rates
  tmp <- tsv_lines[(Grows_line + 1):Population_line-1]
  Grows <- read.table(text = tmp, sep = "\t", skip=1, header=TRUE, row.names = 1)
  
  # Obtener tablas de poblaciones iniciales
  tmp <- tsv_lines[(Population_line + 1):Seed_line-1]
  Pobl <- read.table(text = tmp, sep = "\t", skip=1, header=TRUE, row.names = 1)
  
  result <- list(Interacs, Grows, Pobl)
  return(result)
  
}

#---------------------------Generar simulacion----------------------------------

# Generar las semillas posibles
setwd("~/Documents/LAB_ECO") # Set Working Directory
library (random)
seeds <- randomNumbers(n=300, min=1, max=6000)

ITERS <- 10 # Introducir numero de simulaciones
N <- 20 # Introducir el numero de especies
C <- 0.45 # Probabilidad de interacción=0
times <- seq(0, 100, by = 1) # Numero de generaciones

for (i in 1:ITERS) {
  
  # Generar datos y extraerlos
  res <- generate(N,seeds,C)
  V_inter <- unlist(res[1])
  params <- list(
    r = unlist(res[2]), # Grow rates
    alpha = matrix(V_inter, nrow = N, ncol = N) # Interaction
  )
  Pobl <- unlist(res[3])
  Semilla <- unlist(res[4])
  
  
  # Evaluar la euacion
  library(deSolve)
  output <- ode(y = Pobl, times = times, func = GLV, parms = params,  method = "lsoda", atol = 1e+1, rtol = 1e+1,
                maxsteps = 200)
  
  # Guardar resultados
  ID <- save(output, params, Pobl, Semilla)
  
  # Contar negativos
  Scan_neg(output, params, ID, N, C)
  
}

####
# Plot interacciones negativas (x) especies que terminan en netivas (y)
Scan_plot ()   

# Graficar poblaciones
esp <- c("time", "1","3")
Pgraph(output, esp)

# Graficar red de interaccion
red(params)

# Obtener los parametros nuevamente
# [1] contiene la matriz de interacciones
# [2] contiene los crecimientos
# [3] contiene las poblaciones iniciales
Read_params(ID)

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

# Guardar resultados ANTES QUE SCAN
save(output, params, Pobl, Semilla)

# Contar negativos
scan <- scan_neg(output, params, ID)

# Graficar redes
red(params)

# Grafica poblacion
esp <- c("time", "1","3")
Pgraph(output, esp)

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
Para los growth rate modificarlos para que sean >0 y <1, con runif u otra
Intentar la diagonal con 1 
Intentar con 20 especies
- Para la diagonal intentar con 1.
- ki = 1
- Hacer operaciones en vez de vector como matriz

- 20 generaciones
- Especies constantes
- Por cada red que genere, calcular cuanto son positivos y cuantos son negativos. Tambien lo mismo con la matriz de interaciones

- 20 especies a 100 generaciones, ve estabilidad 

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

- Añadir un parametro C que simule el numero de interacciones promedio que sean diferentes a 0, fuera de la diagonal. 
C es un parametro que yo añado
- FUNCION READ PARA LEER PARAMETROS INICIALES
- Para las poblaciones iniciales generarlas de una uniforme de .1 a 1

- Añadir numero de especies totales en el TSV de scan
- corregir parametro C y como funciona

"