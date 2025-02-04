
# ----------------------------------------Function----------------------------------

CPr_sim <- function(N, C0, CN, times,seeds_n, tol, individual) { 
    
  CPr_time1 <- system.time({
    seeds <- random::randomNumbers(n=seeds_n, min=1, max=100)
  })
    
    #-----------------------------Generate data----------------------------------#
  CPr_time2 <- system.time({
      res <- generate(N,seeds,C0, CN) # Results
      
      # V_inter <- unlist(res[[1]]) 
      params <- list(
        r = unlist(res[2]), # Grow rates
        alpha = matrix(unlist(res[[1]]) , nrow = N, ncol = N) # Interaction
      )
      Pobl <- unlist(res[3])
      Semilla <- unlist(res[4])
  })
  
    #----------------------------------------Simulate----------------------------#
  CPr_time3 <- system.time({
    #library(miaViz)
    interacs <- params$alpha
    glvmodel <- miaSim::simulateGLV(n_species = N, 
                            A = params$alpha, # interaction matrix
                            x0 = Pobl, # Initial abundances
                            growth_rates = params$r, # Growth rates
                            t_start = 0, 
                            t_store = times, 
                            t_end=times, 
                            migration_p = 0,
                            stochastic = FALSE, # Ignorar ruido
                            norm = TRUE) # FALSE=conteo, TRUE=proporciones
    
    output <- glvmodel@assays@data@listData[["counts"]]
  })
    
    #------------------Save Simulation results-----------------------------------#
  CPr_time4 <- system.time({
    ID <- save(output,params,Pobl,Semilla)
  })
  
  #---------------------Search for steady state-----------------------------------#
  CPr_time5 <- system.time({
    out_path <- paste("./Outputs/O_", ID , ".tsv", sep = "") # output
    result <- st_search(ID,out_path,tol,individual)
  })
  
  #-----------------------------------Save Code profiling times------------------#
  CPr_df <- data.frame(
    ID = ID ,
    Specs = N , 
    Generations = times,
    Seed_num = seeds_n,
    'Seeds_time-s' = round(CPr_time1[3], 6) , 
    'GenDat_time-s' = round(CPr_time2[3], 6) ,
    'Simulation_time-s' = round(CPr_time3[3], 6) ,
    'Save_time-s' = round(CPr_time4[3], 6),
    'SSS_time-s' = round(CPr_time5[3], 6),
    Individual = individual
  )
  
  library(readr)
  CPr_path <- paste("./Scan/CPr_time", ".tsv", sep = "") # output
  exist <- file.exists(CPr_path) # Flag
  
  if (!exist) { # File doesnt exist
    
    system(paste("touch", CPr_path)) # Make file
    write.table(CPr_df, file = CPr_path, sep = "\t", row.names = FALSE, col.names = TRUE)  # Save
  } else {
    
    CPr_table <- read.delim(CPr_path, sep = "\t", header = TRUE) # Read table
    Join_CPr <- rbind(CPr_table, CPr_df) # Join tables
    write.table(Join_CPr, file = CPr_path, sep = "\t", row.names = FALSE, col.names = TRUE) # Save
  }
  
  

}


#--------------------------------------Generate simulations--------------------------------

N <- 50 # Number of species
C0 <- 0.45 # Prob. interaction =0
CN <- 0.2 # Prob. interaction <0
times <- 3000 # Generations
counter <- 0
STOP <- FALSE
  
while (!STOP) {
  
  # Run the simulation function
  CPr_sim(N, C0, CN, times)
  
  if (counter < 2) {  # Check if `times` has been updated less than 2 times
    times <- times + 100  # Add 100 generations
    counter <- counter + 1  # Increment `times` update counter
  } else {
    N <- N + 50  # Add 50 species
    counter <- 0  # Reset `times` update counter
  }
  
  # Stop simulation when there are at least 500 species and 10,000 generations
  if (N==500) {
    STOP <- TRUE
  }
  
}
 
  
#------------------------------------------------------Testing----------------------------------------------------------
N <- 50 # Number of species
C0 <- 0.45 # Prob. interaction =0
CN <- 0.2 # Prob. interaction <0
times <- 3000 # Generations
seeds_n <- 300
tol <- 5
individual  <- TRUE
setwd("~/Documents/LAB_ECO") # Set Working Directory

CPr_sim(N,C0,CN,times,seeds_n,tol,individual)
  


unlink("~/Documents/LAB_ECO/Scan/CPr_time.tsv", recursive = TRUE, force = TRUE) # Remove all data from Scan 








########################################################################################################################

    #----------------------------------------------Code profiling "Rprof"---------------------#
    Rprof("./Scan/Cpr_code.out")
    
    testing(N, C0, CN, times)
    
    Rprof(NULL)
    
    summary <- summaryRprof("./Scan/Cpr_code.out")
    print(summary)

    #-----------------------------------------Code profiling "profvis"------------------------#
    library(profvis)
    
    prof_data <- profvis({
      seeds <- random::randomNumbers(n=300, min=1, max=100)
      
      #-----------------------------Generate data----------------------------------#
      res <- generate(N,seeds,C0, CN) # Results
      
      # V_inter <- unlist(res[[1]]) 
      params <- list(
        r = unlist(res[2]), # Grow rates
        alpha = matrix(unlist(res[[1]]) , nrow = N, ncol = N) # Interaction
      )
      Pobl <- unlist(res[3])
      Semilla <- unlist(res[4])
      
      #----------------------------------------Simulate----------------------------#
      #library(miaViz)
      interacs <- params$alpha
      glvmodel <- miaSim::simulateGLV(n_species = N, 
                                      A = params$alpha, # interaction matrix
                                      x0 = Pobl, # Initial abundances
                                      growth_rates = params$r, # Growth rates
                                      t_start = 0, 
                                      t_store = times, 
                                      t_end=times, 
                                      migration_p = 0,
                                      stochastic = FALSE, # Ignorar ruido
                                      norm = TRUE) # FALSE=conteo, TRUE=proporciones
      
      out <- glvmodel@assays@data@listData[["counts"]]
      
      #------------------Save Simulation results-----------------------------------#
      ID <- save(out,params,Pobl,Semilla)
      
    })
    
    prof_info <- prof_data$data
    head(prof_info)
