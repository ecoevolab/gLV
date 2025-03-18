


# Load libraries
library(dplyr) 
library(tidyr)
require(tidyverse, lib.loc = "/mnt/atgc-d3/sur/modules/pkgs/tidyverse_mrc")
library(purrr)

#======= Prepare data ======

# Load parameters data
params_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/ExpM279-D06M03.tsv")
head(params_table)

# Load Na counting table
na_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Results/ExpM279-D06M03/Nas-counting.tsv")
head(na_table)

# Filtering step
na_ids <- na_table$id[na_table$Total.NAs != 0]

df_table <- params_table %>%
  filter(id %in% na_ids)

# Source gLV parameter generation function
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/Method2-normal/generate-gLV-normal.R")
print(regenerate_Dnormal)

#======= Define gLV model ======
# Define the equation
glv_model <- function(t, x0, params) {
  r <- params$mu         # Growth rate vector
  A <- params$M          # Interaction matrix
  
  # Compute dx/dt for each species
  dx <- x0 * (r + A %*% x0)
  list(dx)
}

#======= Example step ======
# Regenerate parameters and apply Dorman-Prince method
tmp <- df_table[1:5,]
apply(tmp, 1, function(x){
  
  # Regenerate parameters
  params <- regenerate_Dnormal(x)
  cat("Parameters generated for simulation: ", params$id, "\n")
  
  # Apply Dorman Prince method
})