#' Declare libraries
library(tidyverse, lib.loc = "/mnt/atgc-d3/sur/modules/pkgs/tidyverse_mrc")
library(tidyr)

#' Generate grid and filter to keep rows where the sum of `p_neg + p_noint <= 1`.
params.to.sim <- expand_grid(
  n_species = seq(5, 100, by = 5),
  p_noint = seq(0, 1, by = 0.05),
  p_neg = seq(0, 1, by = 0.05)
) %>%
  filter(!(p_noint == 1 & p_neg != 0)) %>%
  mutate(id = ids::random_id(n = length(n_species), bytes = 4)) %>%
  mutate(Pop_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE))) %>%
  mutate(Growth_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE))) %>%
  mutate(A_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE)))

#' Verify if ids are unique and in case they are, save the parameters.
if (nrow(params.to.sim) == length(unique(params.to.sim$id) ) ) {
  Params_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Exp03-D25M02.tsv"
  
  # Save Parameters as TSV
  data.table::fwrite(params.to.sim, Params_path, sep = "\t")
}




