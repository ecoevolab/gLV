# Laboratorio Sur
# Ecuaciones de Lotka-Volterra generalizado




GLV_mine <- function (n_species, # NUmber of species
                      A = NULL, # Interaction matrix
                      x0 = NULL, # Initial abundances
                      growth_rates = NULL, 
                      t_external_events = NULL, 
                      t_external_durations = NULL,  
                      error_variance = 0, 
                      norm = FALSE, 
                      metacommunity_probability = NULL,
                      t_end = 1000, ...)  
  {
  
  if (is.null(metacommunity_probability)) {
    metacommunity_probability <- .rdirichlet(1, alpha = rep(1, n_species))
  }
  
  # Generate dirichlet random deviates
  metacommunity_probability <- metacommunity_probability/sum(metacommunity_probability)
  
  # Generate simulation times
  t_dyn <- .simulationTimes(t_end = t_end, ...)
  tEvent <- .eventTimes(t_events = t_external_events, 
                        t_duration = t_external_durations, 
                        t_end = t_end, ... = ...)
  
  parameters <- list(growth_rates = growth_rates, A = A, n_species = n_species, 
                     tEvent = tEvent)
  
  out <- ode(y = x0, times = t_dyn$t_sys, func = glvModel, 
             parms = parameters, events = list(func = perturb, time = t_dyn$t_sys), 
             maxsteps = 10^9, method = "ode45")
  
  out_matrix <- out[, 2:ncol(out)]
  out_matrix <- out_matrix[t_dyn$t_index, ]
  
  if (error_variance > 0) {
    measurement_error <- rnorm(n = length(t_dyn$t_index) * 
                                 n_species, mean = 0, sd = sqrt(error_variance))
    measurement_error <- matrix(measurement_error, nrow = length(t_dyn$t_index))
    out_matrix <- out_matrix + measurement_error
  }
  if (norm) {
    out_matrix <- out_matrix/rowSums(out_matrix)
  }
  colnames(out_matrix) <- names_species
  out_matrix <- cbind(out_matrix, time = t_dyn$t_sys[t_dyn$t_index])
  
  TreeSE <- TreeSummarizedExperiment(assays = list(counts = t(out_matrix[,1:n_species])), 
                                     colData = DataFrame(time = out_matrix[,"time"]), 
                                     metadata = list(migration_p = migration_p, error_variance = error_variance)
                                     )
  return(TreeSE)
}








#--------------------------------------------------------------------Algo-----------------------------------
# Load interactions table
setwd("~/Documents/LAB ECO/Codes")
inter <- read.csv(file = './Interacciones.csv', row.names = 1, header = TRUE)
N <- ncol(inter) # Count number of species

###############################
# Make interactions table OPTIONAL
N <- 3
inter <- matrix(nrow = N, ncol = N)

# Primer numero en [i,0] es fila, segundo numero es columna
for (r in 1:N) { # column
  for (c in 1:N){ # row
    inter [r,c] <- rnorm(1)
  }
}

# Asignar nombres
# Create an empty vector of strings for ROW NAMES
names_vector <- character(0)


for (i in 1:N){ # column
  # Change name column
  name <- paste("Especie_", i , sep = "")
  
  # Add more values to the vector using c()
  names_vector <- c(names_vector, name)
}

colnames(inter) <- names_vector
rownames(inter) <- names_vector

#############

# Vector de especies
Pobl <- vector("numeric", length = N)

for (i in 1:N){
  Pobl[i] <- round(abs(rnorm(1, 10, 5)))
}

# Vector de crecimiento de especies
Grow <- vector("numeric", length = N)

for (i in 1:N){
  Grow[i] <- round(abs(rnorm(1, 1, 1)))
  # Grow[i] <- 0 #PRUEBA
 
}
# Grow[i] <- 1 #prueba
# Pobl[i] <- 0 #prueba

# Comenzar a clacular las fluctuaciones de especies
Final <- Pobl
iters <- 1000


for (t in 1:iters){
  
  for (i in 1:N) {
    
    # Calculo aparte
    Sigma = 0
    for (j in 1:N){
      
      Mij <- inter[i,j] # Obtener interaccion de especie i con j
      Sum <- Mij*Pobl[j]  # Calculo a sumar
      Sigma <- Sigma + Sum # Añadir a la variable sigma
    }
    
    Final[i] <- (Grow[i]*Pobl[i]) + (Final[i]*Sigma)
  }
}

# Integradores numericos, Ordinary diferential equiations (desolve)
# Usar vectores para las interraciones.
# Toda la diagonal de negativos en las interracciones, y fuera de la diagonal algunos
# SEMILLAS
# Guardar en un csv la semillas y las fluctuaciones

# SIGUIENTE PASO: Implementar la extincion
# FOTO Y RESEÑAAAAAAAAAAAAAAAAAAAAA

