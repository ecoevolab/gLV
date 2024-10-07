

run_extinction <- function(nspecs_extinct, output, tolerance, uniqueID, wd) {
  
  
  # Additional parameters
  N_species <- nrow(output)
  
  # Generate original parameters
  params <- params_seed_reader(uniqueID, wd)
  
  # Get on which generation does all the species of the system reach Steady State
  result <- all_prop_SS(output, tolerance)
  
  # Get which species will go extinct and their extinction times
  ext_specs <- round(runif(n = nspecs_extinct, min = 1, max = N_species))
  ext_times <- round(runif(n = nspecs_extinct, min = result - 1, max = ncol(output)))
  
  # ext_specs = c(1,6,5,9,8)
  # ext_times = c(5,6,9,9,4)
  # Combine into a data frame and order by extinction times
  df <- data.frame(ext_specs, ext_times)
  df <- df[order(df$ext_times), ]
  
  # Get unique extinction times and initialize empty data frame for storing results
  unique_times <- unique(df$ext_times)
  diff_times <- diff(unique_times)
  empty_df <- data.frame(matrix(ncol = 0, nrow = N_species))
  
  for (e in 1:length(unique_times) ) {
    
    # Get rows where the generation column has the value t
    df_filtered <- df[df$ext_times == unique_times[e], ]
    
    # Vector of species that extinct at t time
    vec <- df_filtered$ext_specs
    
    cat("Generation ", unique_times[e] , "... specie", vec, "\n")
    
    # Get new population size FIRST TIME
    if (e == 1) {
      population <- output[, unique_times[e]] 
      population[vec] <- 0 # Assign 0 if the specie will extinct
      
      # Modify the population parameters
      params$Population <- population
      
      # Run simulation for the difference in times as t
      test_df <- run_simulation(N_species, params = params, times = 200 )
      
      # Join simulation output to df
      empty_df <- cbind(empty_df, test_df)
      
    } else {
      # Get new population size T TIMES
      if (1 < nspecs_extinct && nspecs_extinct > 2) {
        population <- empty_df[, ncol(empty_df)] 
        population[vec] <- 0 # Assign 0 if the specie will extinct
        
        # Modify the population parameters
        params$Population <- population
        
        # Run simulation for the difference in times as t
        test_df <- run_simulation(N_species, params = params, times = diff_times[e] )
        
        # Join simulation output to df
        empty_df <- cbind(empty_df, test_df)
        
      } 

      if ( 1 < nspecs_extinct && nspecs_extinct == 2 && e == length(unique_times)) {
        population <- empty_df[, ncol(empty_df)]
        population[vec] <- 0 # Assign 0 if the specie will extinct

        # Modify the population parameters
        params$Population <- population

        # Run simulation for the difference in times as t
        test_df <- run_simulation(N_species, params = params, times = 200 )

        #  Join simulation output to df
        empty_df <- cbind(empty_df, test_df)

      }
    }
    #print(e)
  }
  
  colnames(empty_df) <- seq(from = unique_times[1], to = unique_times[1] + ncol(empty_df) - 1)
  empty_df <- as.matrix(empty_df)
  
  #----------------------------------Save extinction----------------------------------------#
  # Extinction path
  ext_path <- file.path(wd, "Extinctions", paste0("E_", uniqueID, ".tsv") )
  
  # Save output data
  utils::write.table(empty_df, file = ext_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  
  
  return(empty_df)
  
  
}











