#' Declare libraries
library(tidyverse, lib.loc = "/mnt/atgc-d3/sur/modules/pkgs/tidyverse_mrc")
library(tidyr)

#' Generate grid and filter to keep rows where the sum of `p_neg + p_noint <= 1`.
Sims <- expand_grid(
  n_species = seq(5, 100, by = 5),
  p_noint = seq(0, 1, by = 0.05),
  p_neg = seq(0, 1, by = 0.05)
) %>%
  filter(p_noint + p_neg <= 1) %>%
  mutate(id = ids::random_id(n = length(n_species), bytes = 4)) %>%
  mutate(Pop_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE))) %>%
  mutate(Growth_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE))) %>%
  mutate(A_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE)))

#' The next lines are for testing if all combinations made meet the criteria of being <= 1.
test <- Sims %>%
  mutate(valid = if_else(p_noint + p_neg <= 1, TRUE, FALSE))
unique(test$valid)

#' Verify if ids are unique and in case they are, save the parameters.
if (nrow(Sims) == length(unique(Sims$id) ) ) {
  Params_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Data-D25M02.tsv"
  
  # Create directory
  dir.create(dirname(Params_path), recursive = TRUE)
  
  # Save Parameters as TSV
  data.table::fwrite(Sims, Params_path, sep = "\t")
  
}




