# ==== Load libraries ====
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
library(tidyverse, lib.loc = "/mnt/atgc-d3/sur/modules/pkgs/tidyverse_mrc")

# ==== Generate parameters ====
param.generate.fun <- function() {
  n_species <- rep(c(10, 30, 100), times = 10)
  
  grid <- expand_grid(
    n_species = n_species,
    Norm_mu = 0,
    sigma = seq(from = 1, to = 20, by = 0.5),
    p_noint = seq(from = 0, to = 1, by = 0.05)
  )
  
  grid %>%
    mutate(
      id = ids::random_id(nrow(grid), bytes = 4),
      x0_seed = sample(1:1e6, nrow(grid), replace = FALSE),
      mu_seed = sample(1:1e6, nrow(grid), replace = FALSE),
      A_seed = sample(1:1e6, nrow(grid), replace = FALSE),
      diagonal = -0.5
    )
}

params.to.sim <- param.generate.fun()
head(params.to.sim)


# ==== Save parameters ====
generate_id <- function() {
  day <- format(Sys.Date(), "%d")
  month <- format(Sys.Date(), "%m")
  
  char_string <- paste0(sample(c(LETTERS, 0:9), 4, replace = TRUE), collapse = "")
  paste0("Exp", char_string, "-D", day, "M", month)
}

# Set working directory
wd <- "/mnt/atgc-d3/sur/users/mrivera/glv-research"
exp.name <- generate_id()

# Paths
params.path <- file.path(wd, "Data", paste0(exp.name, ".tsv"))
results.path <- file.path(wd, "Results", exp.name)

# Ensure unique parameters
save_params <- function(params.to.sim, params.path, results.path) {
  if (!file.exists(params.path)) {
    data.table::fwrite(params.to.sim, params.path, sep = "\t")
    dir.create(results.path, recursive = TRUE, showWarnings = FALSE)
    message("Parameters successfully saved to: ", params.path)
    return(TRUE)
  } else {
    message("File already exists. Verify path...")
    return(FALSE)
  }
}

if (nrow(params.to.sim) != length(unique(params.to.sim$id))) {
  repeat {
    params.to.sim <- param.generate.fun()
    if (nrow(params.to.sim) == length(unique(params.to.sim$id))) {
      if (save_params(params.to.sim, params.path, results.path)) break
    } else {
      warning("Duplicate IDs detected. Parameters not saved.")
    }
  }
} else {
  save_params(params.to.sim, params.path, results.path)
}
