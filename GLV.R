# Laboratorio Sur
# Ecuaciones de Lotka-Volterra generalizado


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

