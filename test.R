
# Test using pracma::ode23s

# Antes de utilizar este archivo correr las funciones de "DE_LV.R"
#-----------------------------------------------------------------------------------------------------------------------

# Generar las semillas posibles
setwd("~/Documents/LAB_ECO") # Set Working Directory
library (random)
seeds <- randomNumbers(n=300, min=1, max=100)

N <- 20 # Number of species
C0 <- 0.45 # Prob. interaction =0
CN <- 0.2 # Prob. interaction <0
times <- seq(0, 200, by = 1) # Generations

# Generar datos y extraerlos
res <- generate(N,seeds,C0, CN) # Results

# V_inter <- unlist(res[[1]]) 
params <- list(
  r = unlist(res[2]), # Grow rates
  alpha = matrix(unlist(res[[1]]) , nrow = N, ncol = N) # Interaction
)
Pobl <- unlist(res[3])
Semilla <- unlist(res[4])

# Utilizar nueva función para ver funcionamiento miaSim
# 


#-----------------------PRACMA------------------------------------------------------------------------------------------

# y <- Poblaciones iniciales
Prac_GLV <- function(t,y, ...) {
  (r * y ) + ( (alpha %*% y) * y ) 
  
}

library(pracma)
alpha <- as.matrix(params$alpha)
r <- params$r
y0 <- as.matrix(Pobl)
t0 <-0; tf <-1000; 
t <- seq(from=1, to=100, by= 1)
t_points = linspace(0, 100, 80);

output3 <- pracma::ode23s(Prac_GLV, t0, tf, y0 ,  rtol = 1e+1, atol= 1e+1, ...=c(alpha, r, t_points))

# Organizar en un dataframe
times <- output3[["t"]]
out <- output3[["y"]]
MAT_out <- as.matrix(out, nrows=length(times))
rownames(MAT_out) <- times # Los nombres de los renglones son los tiempos

ID <- save(MAT_out, params, Pobl, Semilla)
Scan_neg(MAT_out, params, ID, N, C)

# Graficar

# Population Graph
library(reshape2)
library(ggplot2)
library(tidyverse)

result_matrix <- cbind(MAT_out, times) #Add column of time
result_matrix <- as.data.frame(result_matrix)
df_long <- result_matrix %>% gather(key = "species", value = "population", -times)
tmp <- ggplot(df_long, aes(x = times, y = population, color = species)) + geom_line(size = 1.5) + 
  labs(title = "Population Over Time", x = "Time", y = "Population", color = "Species") +
  ylim(c(-100,200)) +
  xlim(c(0, 150)) +
  theme_minimal() 

# jpeg(file="~/Documents/LAB_ECO/Poster/image1.jpeg")
# tmp
# dev.off()

#------------------------------------------------------miaSim-----------------------------------------------------------

# Para leer parametros primero correr la funcion "R_params"
# Read paremeters
setwd("~/Documents/LAB_ECO")
ID <- "826ba2"
R_params <- Read_params(ID)

r <- unlist(R_params[2]) # Grows
alpha <- as.matrix(R_params[[1]]) # Interactions
Pobl <- unlist(R_params[3]) # Poblation

#---------------------
library(miaSim)
library(miaViz)
glvmodel <- simulateGLV(n_species = 20, 
                        A = alpha, # interaction matrix
                        x0 = Pobl, # Initial abundances
                        growth_rates = r, # Growth rates
                        t_start = 0, 
                        t_store = 10, 
                        t_end=10, 
                        migration_p = 0,
                        stochastic = FALSE, # Ignorar ruido
                        norm = TRUE) # FALSE=conteo, TRUE=proporciones

out <- glvmodel@assays@data@listData[["counts"]]

miaViz::plotSeries(glvmodel, "time")



