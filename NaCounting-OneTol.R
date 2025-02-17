cat(
  paste0(rep("=", 20), collapse = ""), "  Running code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)

library(dplyr)
library(data.table) 

# Define the root path
root_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/D13M02Y25-outs"

# Función para calcular la cantidad de NAs en un archivo
calculate_nas <- function(file_path) {
  sum(is.na(data.table::fread(file_path)))
}

# Obtener la lista de archivos
files <- list.files(root_path, full.names = TRUE)

# Calcular NAs y extraer IDs en un solo paso
na_counts_list <- lapply(files, function(file) {
  data.frame(
    TSV_ID = sub(".*_(.*?)\\.tsv$", "\\1", basename(file)),  # Extraer ID
    NA_Counts = calculate_nas(file),
    stringsAsFactors = FALSE
  )
})

# Unir todas las filas en un solo data frame
final_df <- bind_rows(na_counts_list)

# Optionally, save the final data frame to a CSV
save_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/D13M02Y25/NAcount-D13M02Y25.tsv"
data.table::fwrite(final_df, file = save_path, sep = "\t")

cat(
  paste0(rep("=", 20), collapse = ""), "  Ending code at: ", 
  format(Sys.time(), "%B %d, %Y %I:%M:%S %p"), " ", 
  paste0(rep("=", 20), collapse = ""), 
  "\n"
)


