
#------------------------------Limpiar carpetas---------------------------------

# Outputs
unlink("~/Documents/LAB_ECO/Outputs/*", recursive = TRUE, force = TRUE)
# Parameters
unlink("~/Documents/LAB_ECO/Parameters/*", recursive = TRUE, force = TRUE)
# Scan
unlink("~/Documents/LAB_ECO/Scan/*", recursive = TRUE, force = TRUE)

unlink("~/Documents/LAB_ECO/testing/*", recursive = TRUE, force = TRUE) # test

#------------------------------------Extra code---------------------------------
# Generar las semillas posibles
setwd("~/Documents/LAB_ECO") # Set Working Directory


# Outdated graph
out <- glvmodel@assays@data@listData[["counts"]]
miaViz::plotSeries(glvmodel, "time")

####
# Compare
negative_count <- sum(interacs == 0)
count_intrf = sum(interacs > 0)
count_test = sum(interacs < 0)

  
            
