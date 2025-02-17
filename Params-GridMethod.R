


#------------------------Prepare grid cell for simulations--------#
library(tidyverse, lib.loc = "/mnt/atgc-d3/sur/modules/pkgs/tidyverse_mrc")

Sims <- expand_grid(n_species = seq(from = 5, to = 100, by = 5),
                    p_noint = seq(from = .1, to = 1, by = 0.05),
                    p_neg = seq(from = .1, to = 1, by = 0.05)) %>%
  mutate(id = ids::random_id(n = length(n_species), bytes = 4)) %>%
  mutate(Pop_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE))) %>%
  mutate(Growth_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE))) %>%
  mutate(A_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE)))


if (nrow(Sims) == length(unique(Sims$id) ) ) {
  Params_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/D13M02Y25.tsv"
  
  # Create directory
  dir.create(dirname(Params_path), recursive = TRUE)
  
  # Save Parameters as TSV
  data.table::fwrite(Sims, Params_path, sep = "\t")
  
}




