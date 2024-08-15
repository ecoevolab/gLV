extinct <- function(ID, t_ext, S_ext, t_gens) {
  
        # Alert
        if (t_ext==0){
          print(paste("The time of extinction (", t_ext, ")", "hast to be greater than 0"))
          return()
        }
        
        #-----------------------Read parameters--------------#
        R_params <- Read_params(ID)
        r <- unlist(R_params[2]) # Grows
        alpha <- as.matrix(R_params[[1]]) # Interactions
        N <- ncol(alpha)
        
        # Read output
        file <- paste("./Outputs/O_",ID,".tsv", sep = "")
        output <- read.csv(file, sep = "\t")
        
        # Alert 2 
        if (t_ext>ncol(output)){
          print(paste("The time of extinction (", t_ext, ")", "is not valid it is greater than the number of times simulated(", ncol(output) ,")"))
          return()
        }
        
        #----------------------Change parameters-----------#
        # Columns-> times 
        # Rows -> species
        if (t_ext-1>0){ 
          Pobl <- output[,t_ext-1] # The extinction time is >1
        } else {
          Pobl <- output[,t_ext] # The extinction time is 1
        }
        
        Pobl[S_ext] <- 0 # Asign 0 to the specie in the population
      
        #---------------------New simulation---------------#
      library(miaSim)
      glvmodel <- simulateGLV(n_species = N, 
                              A = alpha, # interaction matrix
                              x0 = Pobl, # Initial abundances
                              growth_rates = r, # Growth rates
                              t_start = 0, 
                              t_store = t_gens, 
                              t_end=t_gens, 
                              migration_p = 0,
                              stochastic = FALSE, # Ignorar ruido
                              norm = TRUE) # FALSE=conteo, TRUE=proporciones
      
      out <- glvmodel@assays@data@listData[["counts"]]
      
      #--------------------Save new simulation-------------#
      f_out <- paste("./Extinction/E_", ID ,"_S0", S_ext, "-T0", t_ext, ".tsv", sep = "") # output
      
      # Revisar si un archivo con ese ID existe
      exist <- file.exists(f_out) # Bandera
      if (exist==TRUE) { 
        print(paste("A file with the same conditions already exist, exiting..."))
        return()
      }
      
      # Generar los archivos
      system(paste("touch", f_out))
      
      # Guardar la tabla OUTPUT
      write.table(out, file = f_out, sep = "\t", row.names = FALSE, col.names = TRUE)
}


#----------------------------------------Extinction-----------------------------#
t_ext=3
S_ext=5
t_gens=25
extinct(ID, t_ext, S_ext, t_gens)