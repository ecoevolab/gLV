#' # Exp04 using Normal Distribution for interactions
#'
#' This script generates a parameter grid, applies filters, verifies unique IDs, and saves the output if valid.

#' ## Attach Packages
suppressPackageStartupMessages(library(dplyr))

#' Declare libraries
suppressPackageStartupMessages(library(tidyverse, lib.loc = "/mnt/atgc-d3/sur/modules/pkgs/tidyverse_mrc"))
suppressPackageStartupMessages(library(tidyr))

#' ## Generate Parameters for simulation
#'
#' Create a grid of parameters.
#'
#' The grid includes:
#' - Species counts (`n_species`)
#' - Normal distribution mean (`Norm_mu`)
#' - Normal distribution standard deviation (`sigma`)
#' - No interaction probability (`p_noint`)

#+ generate_grid, eval=FALSE
# Generate grid
params.to.sim <- expand_grid(
  n_species = rep(x = c(10,100), times = 10),
  Norm_mu = 0,
  sigma = seq(from = 0.005, to = 1, by = 0.05),
  p_noint = seq(0, 1, by = 0.05)
) %>%
  mutate(id = ids::random_id(n = length(n_species), bytes = 4)) %>%
  mutate(x0_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE))) %>%
  mutate(mu_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE))) %>%
  mutate(A_seed = as.vector(sample(1:1e6, length(n_species), replace = FALSE)))

#' Verify Unique IDs and Save. If all generated IDs are unique, the parameters are saved as a TSV file.

#+ verify_and_save, eval=FALSE
if (nrow(params.to.sim) == length(unique(params.to.sim$id))) {
  Params_path <- "/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Exp04-NormalD.tsv"
  
  # Save Parameters as TSV
  data.table::fwrite(params.to.sim, Params_path, sep = "\t")
  message("Parameters successfully saved to: ", Params_path)
} else {
  warning("Duplicate IDs detected. Parameters not saved.")
}
