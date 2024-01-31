
# Test using pracma::ode23s

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



#-----------------------PRACMA------------------------------------------------------------------------------------------

# y <- Poblaciones iniciales
Prac_GLV <- function(t,y, ...) {
  (r * y ) + ( (alpha %*% y) * y ) 
  
}

library(pracma)
alpha <- as.matrix(params$alpha)
r <- params$r
y0 <- as.matrix(Pobl)
t0 <-0; tf <-100; 
t <- seq(from=1, to=100, by= 1)
t_points = linspace(0, 100, 80);

output3 <- pracma::ode23s(Prac_GLV, t0, tf, y0, t_points, rtol = 1e+1, atol= 1e+1, hmax=300, ...=c(alpha, r,))

# Organizar en un dataframe
times <- output3[["t"]]
out <- output3[["y"]]
MAT_out <- as.matrix(out, nrows=length(times))
rownames(MAT_out) <- times # Los nombres de los renglones son los tiempos

