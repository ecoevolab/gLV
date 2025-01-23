# (C) Copyright 2025 Sur Herrera Paredes

# This file is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This file is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this file  If not, see <http://www.gnu.org/licenses/>.

library(tidyverse)
library(miaSim)
library(miaViz)

generate_params <- function(n_species = 20,
                            p_noint = 0.5,
                            p_neg = 0.5){
  
  # Simulate proportion of interactions
  n_int <- n_species ^ 2 - n_species
  n_zero <- round(n_int * p_noint)
  n_neg <- round((n_int - n_zero) * p_neg)
  signs <- rep(c(-1, 0, 1), times = c(n_neg, n_zero, n_int - n_neg - n_zero))
  
  
  # Randomize parameters
  x0 <- runif(n = n_species, min = 0, max = 1)
  mu <- runif(n = n_species, min = 0, max = 1)
  signs <- sample(signs, size = n_int, replace = FALSE)
  
  # Interaction matrix
  M <- matrix(NA, ncol = n_species, nrow = n_species)
  ii_diag <- diag(n_species) == 1
  M[ ii_diag ] <- rep(-0.5, times = n_species)
  M[ ! ii_diag  ] <- runif(n = n_int, min = 0, max = 0.5) * signs

  # M <- randomA(n_species = n_species,
  #              diagonal = -0.5,
  #              connectance = 1 - p_noint,
  #              scale_off_diagonal = 1,
  #              mutualism = 1,
  #              commensalism = 1,
  #              parasitism = 1,
  #              amensalism = 1, 
  #              competition = 1,
  #              interactions = runif(n = n_int, min = 0, max = 0.5) * signs,
  #              symmetric = FALSE
  # )
  
  return(list(x0 = x0, mu = mu, M = M))
}


sim_glv <- function(params = params, n_t = n_t){
  
  # Check that matrrix is square
  if(nrow(params$M) != ncol(params$M))
    stop("ERROR", call. = TRUE)
  
  # Use miaSim to simulate standard gLV
  msim <- simulateGLV(n_species = nrow(params$M), 
                      A = params$M,
                      x0 = params$x0,
                      growth_rates = params$mu,
                      sigma_migration = 0,
                      epoch_p = 0,
                      t_external_events = NULL,
                      t_external_durations = NULL,
                      stochastic = FALSE,
                      migration_p = 0,
                      error_variance = 0,
                      norm = FALSE,
                      t_end = n_t)
  
  # Check if NA's or Inf's
  if(any(colSums(is.na(assay(msim))) != 0))
    stop("ERROR: NAs", call. = TRUE)
  if(any(colSums(is.infinite(assay(msim))) != 0))
    stop("ERROR: Infinites", call. = TRUE)
  
  return(msim)
}

measure_start_end_change <- function(msim){
  sim <- (assay(msim))
  t_end <- dim(sim)[2]
  
  sim.prcomp <- prcomp(sim, center = TRUE, scale. = TRUE)
  
  Res <- tibble(shannon_start = vegan::diversity(sim[,1]),
                shannon_end = vegan::diversity(sim[,t_end]),
                richness_start = sum(sim[,1] > 0),
                richness_end = sum(sim[,t_end] > 0),
                bray_change = vegan::vegdist(t(sim[,c(1,t_end)]), method = "bray"),
                euclidean_change = vegan::vegdist(t(sim[,c(1,t_end)]), method = "euclidean"),
                manhattan_change = vegan::vegdist(t(sim[,c(1,t_end)]), method = "manhattan"),
                pc1_change = abs(sim.prcomp$rotation[1,"PC1"] - sim.prcomp$rotation[t_end,"PC1"]),
                pc1_propvar = summary(sim.prcomp)$sdev[1]^2 / sum(summary(sim.prcomp)$sdev^2),
                pc_totvar = sum(summary(sim.prcomp)$sdev^2)
  )
  
  
  return(Res)
}



n_species <- 20
p_noint <- 0.2
p_neg <- 0.8
n_t <- 1000
seed <- 7654
id <- "aaaaaaaa" # NEED TO CHANGE THIS!!!!!

set.seed(seed)
params <- generate_params(n_species = n_species,
                          p_noint = p_noint, 
                          p_neg = p_neg)

msim <- sim_glv(params = params, n_t = n_t)


# colSums(is.na(assay(msim)))
# colSums(is.infinite(assay(msim)))
# summary(t(assay(msim)))
# plotAbundanceDensity(msim)
# plotSeries(msim, x = "time") 
# warnings()

# Identify surviving species
x_t <- assay(msim)[,n_t]
surv_specs <- which(x_t > 0)

Sims <- tibble(id = id,
               n_species = n_species,
               p_noint = p_noint,
               p_neg = p_neg,
               n_t = n_t,
               seed = seed,
               extinct_species = NA,
               spec_abun = NA,
               spec_freq = NA,
               params = list(params),
               sim = list(assay(msim))) %>%
  bind_cols(measure_start_end_change(msim))

# Simulate each species extinction
for(spec in surv_specs){
  # spec <- surv_specs[1]
  cat("Extinguishing species ", spec, "\n")
  
  # Simulate extinction
  x_e <- x_t
  x_e <- x_e[-spec] 
  
  mu_e <- params$mu
  mu_e <- mu_e[-spec]
  
  M_e <- params$M
  M_e <- M_e[-spec, -spec]
  
  # New M params
  n_species_e <- nrow(M_e)
  ii_diag <- diag(n_species_e) == 1
  n_int <- n_species_e ^ 2 - n_species_e
  p_noint_e <- sum(M_e[!ii_diag] == 0) / n_int
  p_neg_e <- sum((M_e[!ii_diag] != 0) & (M_e[!ii_diag] < 0)) / (n_int - sum(M_e[!ii_diag] == 0))

  # Update parameters
  params_e <- list(x0 = x_e, mu = mu_e, M = M_e)
  
  # Simulate time after extinction
  msim_e <- sim_glv(params = params_e, n_t = n_t)
  
  res <- tibble(id = id,
                n_species = n_species_e,
                p_noint = p_noint_e,
                p_neg = p_neg_e,
                n_t = n_t,
                seed = seed,
                extinct_species = spec,
                spec_abun = x_t[spec],
                spec_freq = x_t[spec] / sum(x_t),
                params = list(params_e),
                sim = list(assay(msim_e))) %>%
    bind_cols(measure_start_end_change(msim_e))

  
  Sims <- bind_rows(Sims, res)
  res <- NULL
}

save(Sims, file = paste0(id, ".rdat") )

Sims %>%
  select(-params, -sim) %>%
  write_tsv(file = paste0(id, ".tsv"))


Sims %>%
  filter(!is.na(extinct_species)) %>%
  ggplot(aes(x = spec_abun, y = shannon_end - shannon_start)) +
  geom_point() +
  theme_classic()
