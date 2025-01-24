#!/usr/bin/env Rscript
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

library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Simulate extinctions with gLV" ))
  
  # Positional arguments
  p <- add_argument(p, "sims",
                    help = paste("Table of simulations. Needs columns: p_neg,",
                                 "id, seed, n_species, and p_noint"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--outdir",
                    help = paste("Output directory"),
                    type = "character",
                    default = "output/")
  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  args$rdats <- file.path(args$outdir, "rdats")
  args$tsvs <- file.path(args$outdir, "tsvs")
  
  return(args)
}

args <- process_arguments()
args <- list()
args$sims <- "/Users/sur/lab/exp/2025/today/simulations_batch1.tsv"
args$outdir <- "output"
args$rdats <- file.path(args$outdir, "rdats")
args$tsvs <- file.path(args$outdir, "tsvs")
print(args)

library(tidyverse)
library(miaSim)
library(miaViz)


#' Generate parameters for gLV simulation
#' 
#' Uses overall parameter constraints to generate starting condition
#' and specific sets of parameters
#' 
#' Sets diagonal elements to -0.5, and interaction values capped at an absolute
#' value of 0.5. Based on miaSim
#'
#' @param n_species Number of species to simulate
#' @param p_noint Probability of no interaction, i.e. probability that a
#' non-diagonal element in the interaction matrix is zero. 
#' @param p_neg Probability of a negative effect, i.e. probability that
#' a non-diagonal non-zero element of the interaction matrix is negative
#'
#' @return A list with starting conditions (x0), growth rates (mu), and
#' interaction matrix (M) for gLV simulation
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
  
  return(list(x0 = x0, mu = mu, M = M))
}

#' Simulate gLV
#' 
#' Wrapper for miaSim. Takes output from enerate_params and runs
#' gLV simulation. No stochasticity, measurement error, perturbation or
#' immigration is considered.
#' 
#' Checks for numerical issues after the simulation
#'
#' @param params A list with starting conditions (x0), growth rates (mu), and
#' interaction matrix (M) for gLV simulation
#' @param n_t Number of timepoints to simulate
#'
#' @return A TreeSummarizedExperiment class object produced by 
#' miaSim::simulateGLV
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

#' Compares start and end
#' 
#' Makes some quick calculation between the start and end timepoints of
#' a gLV simulation produced with miaSim::simulateGLV
#'
#' @param msim A TreeSummarizedExperiment class object produced by 
#' miaSim::simulateGLV
#'
#' @return A tibble with various diversity and change estimates
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

#' Simulate all extinction
#' 
#' For a given simulation and parameters, simulate the extinction of all
#' surviving species for the same time as the original simulation
#'
#' @param msim A TreeSummarizedExperiment class object produced by 
#' miaSim::simulateGLV
#' @param params  list with starting conditions (x0), growth rates (mu), and
#' interaction matrix (M) for gLV simulation
#'
#' @return A tibble with the results of all the extinction simulations
simulate_all_extinctions <- function(msim, params){
  # Identify surviving species
  x_t <- assay(msim)[,n_t]
  surv_specs <- which(x_t > 0)
  
  # Simulate each species extinction
  Ext <- NULL
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
    
    
    Ext <- bind_rows(Ext, res)
    res <- NULL
  }
  
  return(Ext)
}



if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
  
  if(!dir.exists(args$rdats)){
    dir.create(args$rdats)
  }
  
  if(!dir.exists(args$tsvs)){
    dir.create(args$tsvs)
  }
}

date()

hyper <- read_tsv(args$sims)
read_tsv(args$sims) %>%
  select(id, n_species, p_noint, p_neg, seed) %>%
  pmap(.f = function(id, n_species, p_noint, p_neg, seed){
    
    i <- 2
    n_species <- hyper$n_species[i]
    p_noint <- hyper$p_noint[i]
    p_neg <- hyper$p_neg[i]
    seed <- hyper$seed[i]
    id <- hyper$id[i]
    
    message(paste0("Simulation ", id))
    n_t <- 1000
    
    # Generate parameters and simulate
    set.seed(seed)
    params <- generate_params(n_species = n_species,
                              p_noint = p_noint, 
                              p_neg = p_neg)
    
    # print(params)
    msim <- sim_glv(params = params, n_t = n_t)
    
    # Combine results in tibble
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
    
    # Simulate extinctions
    Ext <- simulate_all_extinctions(msim = msim, params = params)
    Sims <- bind_rows(Sims, Ext)
    
    # Save full results and table
    save(Sims, file = file.path(args$rdats, paste0(id, ".rdat") ))
    Sims %>%
      select(-params, -sim) %>%
      write_tsv(file = file.path(args$tsvs, paste0(id, ".tsv")) )
  })

date()
sessionInfo()
