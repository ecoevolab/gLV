

#--------------------------------------------------------Testing-----------------------------------------------------

# wd <- "~/Documents/LAB_ECO/"
wd <- "/mnt/atgc-d3/sur/users/mrivera"
setwd(wd)

#-----------------------------------Date of simulation-----------------#
cat("\n" , "\n", strrep("#", 25))
current_date <- Sys.Date()
formatted_date <- format(current_date, "%d %B %Y")
cat("Simulation started at:", formatted_date , strrep("#", 25), "\n")

#----------------------------------Simultaion data--------------------#
n_sim = 5 # Number of Simulations
counter <- 1 
N = 5  # Number of species
times = 50 # Generations

# source("./gLV/Functions.R")
source("./Functions.R")

for (sim in 1:n_sim) {
  # Run the simulation function
  CPr_sim(N,  # Number of species
          C0 = 0.45, # Prob. interaction =0 
          CN = 0.2, # Prob. interaction <0
          times, # Generations 
          tol = 0.05, # Tolerance 
          individual = FALSE, 
          V_diag = -0.5, # Diagonal interaction
          wd =  "~/Documents/LAB_ECO/")
  
  cat("Number of generations", times , "\n")
  cat("Species number", N , "\n")
  cat("Simulation number", sim, "DONE\n\n")
  cat("-----------------------------------------", "\n")
  
  if (counter >= 10) {  # Check if `times` has been updated 3 times
    N <- (N + sample(1:200, 1))    # Add X random species 
    counter <- 1  # Reset `times` update counter 
  } else {
    times <- times * 10  # Add 100 generations
    counter <- counter + 1  # Increment `times` update counter
  }
}


