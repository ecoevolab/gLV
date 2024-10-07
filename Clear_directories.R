
#------------------------------Limpiar carpetas---------------------------------

# Outputs
unlink("~/Documents/LAB_ECO/Simulations/Outputs/*", recursive = TRUE, force = TRUE)
# Parameters
unlink("~/Documents/LAB_ECO/Simulations/Parameters/*", recursive = TRUE, force = TRUE)
# Scan
unlink("~/Documents/LAB_ECO/Simulations/Scan/*", recursive = TRUE, force = TRUE)
# Differences
unlink("~/Documents/LAB_ECO/Simulations/Differences/*", recursive = TRUE, force = TRUE)

unlink("~/Documents/LAB_ECO/testing/*", recursive = TRUE, force = TRUE) # test

#------------------------------------Extra code---------------------------------
# Generar las semillas posibles
setwd("~/Documents/LAB_ECO") # Set Working Directory


# Outdated graph
out <- glvmodel@assays@data@listData[["counts"]]
miaViz::plotSeries(glvmodel, "time")

#------------------------------
# Append new entry or initialize new data
if (file.exists(RDS_path)) {
  
  # Load existing data
  old_list <- readRDS(RDS_path)  
  
  # Look if the ID is already on the list and check if the new search is with a different tolerance
  all_ids <- sapply(tmp, function(x) x$ID)
  id_found <- uniqueID %in% all_ids
  
  # If the ID is present append new data to it
  if (id_found) {
    index <- which(uniqueID == all_ids)
    old_list[[index +  1]]
    #old_data <- old_list$data[[index]]
  } else {
    # Append new entry
    updated_entry <- list(old_list, nested_entry)
    
    # Save appended entry
    saveRDS(updated_entry, file = RDS_path) 
    cat("ID not found on RDS file, adding ID...")
  }
  
} else {
  
  # file doesn't exist, so we save it 
  saveRDS(nested_entry, file = RDS_path) 
  cat("RDS file generated, it didnt exist")
}

  
            
