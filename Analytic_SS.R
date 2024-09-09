
wd <- "~/Documents/LAB_ECO"
setwd(wd)
source("./gLV/Functions.R")

#--------------------------Load all seeds-------------------------------#
library(data.table)
seeds <- as.matrix(fread("./Seeds.tsv", sep = "\t"))

#-----------------------------Generate data----------------------------------#
N <- 2 # Number of species
C0 <- 0.45 # Prob. interaction =0
CN <- 0.2 # Prob. interaction <0
tol <- 0.05 # Tolerance
V_diag <- -0.5

res <- generate(N =2, 
                seeds, 
                C0 <- 0.45, # Prob. interaction =0, 
                CN <- 0.2, # Prob. interaction <0, 
                V_diag <- -0.5)

params <- list(
  alpha = res$Interacs, # Interactions
  r = res$Growth, # Growth rates
  Pobl = res$Population, # Population
  Semilla = res$Seeds # Seeds
)

#----------------------------Simulate ------------------------------------#
glvmodel <- miaSim::simulateGLV(
  n_species = N, 
  A= params$alpha, 
  x0 = params$Pobl, 
  growth_rates = params$r, 
  #t_start = 0, 
  #t_store = times, 
  t_end = 700,
  migration_p = 0, 
  stochastic = FALSE, 
  norm = FALSE
)

output <- glvmodel@assays@data@listData[["counts"]]
cat("Primer renglon columnas 1:4:", output[1, 1:4], "\n")

#-------------------------- Save results---------------------------------#
ID <- Ms_save(N, C0, CN, V_diag, output, params)
cat("El numero de ID es:", ID, "\n")


#-------------------------find steady state------------------------------#
source("./gLV/Steady_fun.R")
ID = "e33382"
tmp <- All_SS_save(ID, tol, params, wd)



#------------------------Testing---------------------------------------#

A_inv <- solve(as.matrix(params$alpha) ) # Inverso de la matriz
detA <- det(params$alpha)
test = A_inv %*% params$alpha
P <- -A_inv %*% params$r
cat("Steady state populations: \n" , P)

col_index <- which(abs(output[1, ] - P[1]) < 0.05)[1]

col_index <- which(abs(output[2, ] - P[2]) < 0.05)[1]
