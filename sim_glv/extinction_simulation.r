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
  p <- add_argument(p, "--ks_boost",
                    help = paste("Boost factor for keystone speice. Zero does",
                                 "nothing."),
                    type = "numeric",
                    default = 0)
  

  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  args$rdats <- file.path(args$outdir, "rdats")
  args$tsvs <- file.path(args$outdir, "tsvs")
  
  return(args)
}

# args <- process_arguments()
setwd("/home/sur/lab/exp/2026/today2/")
args <- list()
args$sims <- "/home/sur/lab/exp/2026/today2/glv_simulation_parameters_controls.tsv"
args$outdir <- "sims"
args$ks_boost <- 10
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
#' @param ks_boost Indicates whether to and how much to boost the interaction parameters
#' for 1 'keystone' species. If zero nothing is done. If greater than zero,
#' then the effects of a random species are multiplied by the specified number.
#'
#' @return A list with starting conditions (x0), growth rates (mu),
#' interaction matrix (M) for gLV simulation, keystone species index (ks_ii),
#' and kesytonge boost factor (ks_boost)
#' 
#' @export
generate_params <- function(n_species = 20,
                            p_noint = 0.5,
                            p_neg = 0.5,
                            ks_boost = 0){
  
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
  M[ ! ii_diag  ] <- runif(n = n_int, min = 0, max = 0.5) * signs
  
  # Use keystone boost
  if(ks_boost > 0){
    ks_ii <- sample(n_species, size = 1)
    M[ , ks_ii ] <- ks_boost *  M[ , ks_ii ] 
  }else{
    ks_ii <- NA
  }
  
  # Set diagonal
  M[ ii_diag ] <- rep(-0.5, times = n_species)
  
  return(list(x0 = x0, mu = mu, M = M, ks_ii = ks_ii, ks_boost = ks_boost))
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
#' @param timeout Number of seconds to wait for simulation to end before
#' stopping it
#'
#' @return A tibble where column sim has the results of the simulation
#' 
#' @export
sim_glv <- function(params = params, n_t = n_t, timeout = 600){

  # Check that matrix is square
  if(nrow(params$M) != ncol(params$M))
    stop("ERROR", call. = TRUE)
  
  # Use miaSim to simulate standard gLV
  # Added timeout for dealing with rare instance where simulation
  # keeps going forever. An issue is that we cannot distinguish failure
  # by timeout from other types of failure.
  msim <- tryCatch(
    R.utils::withTimeout(msim <- simulateGLV(n_species = nrow(params$M), 
                                             names_species = names(params$x0),
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
                                             t_end = n_t),
                         timeout = timeout), 
    error = function(e){
      cat(">>Simulation failed...skipping")
      matrix(NA)})
  sim <- assay(msim)
  
  # We need to define extinctions, miasim sets abundances to zero below
  # 10^-8 to calculate derivative. Given that differential equation
  # is defined as:
  # dx <- x0 * (growth_rates + A %*% x0)
  # This means that once abundance goes below threshold it cannot change
  # anymore, being effectively extinct.
  sim[ sim < 10^-8 ] <- 0
  
  # Check if NA's or Inf's
  if( (any(colSums(is.na(sim)) != 0)) || (any(colSums(is.infinite(sim)) != 0))){
    sim <- NA
  }
  
  # M params
  n_species <- nrow(params$M)
  ii_diag <- diag(n_species) == 1
  n_int <- n_species ^ 2 - n_species
  p_noint <- sum(params$M[!ii_diag] == 0) / n_int
  p_neg <- sum((params$M[!ii_diag] != 0) & (params$M[!ii_diag] < 0)) / (n_int - sum(params$M[!ii_diag] == 0))
  
  # Sim
  
  
  # Combine results in tibble
  Sim <- tibble(n_species = n_species,
                p_noint = p_noint,
                p_neg = p_neg,
                n_t = n_t,
                ks_ii = params$ks_ii,
                extinct_species = NA,
                spec_abun = NA,
                spec_freq = NA,
                params = list(params),
                sim = list(sim)) %>%
    bind_cols(measure_start_end_change(sim))
  
  return(Sim)
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
#' 
#' @export
measure_start_end_change <- function(msim){
  
  if(is(msim, "TreeSummarizedExperiment")){
    sim <- assay(msim)
  }else if(is(msim, "matrix")){
    # Nothing specific
    sim <- msim
  }else if(is(msim, "logical") && is.na(msim)){
    return(NULL)
  }else{
    stop("ERROR measure_change", call. = TRUE)
  }
  
  # Identify number of timepoints
  t_end <- dim(sim)[2]
  
  # Calculate relative abundance
  sim <- apply(sim, 2, function(x) {x / sum(x)})
  
  # sim.prcomp <- prcomp(sim, center = TRUE, scale. = TRUE)
  
  Res <- tibble(shannon_start = vegan::diversity(sim[,1]),
                shannon_end = vegan::diversity(sim[,t_end]),
                richness_start = sum(sim[,1] > 0),
                richness_end = sum(sim[,t_end] > 0),
                bray_change = as.numeric(vegan::vegdist(t(sim[,c(1,t_end)]), method = "bray"))
                # euclidean_change = as.numeric(vegan::vegdist(t(sim[,c(1,t_end)]), method = "euclidean")),
                # manhattan_change = as.numeric(vegan::vegdist(t(sim[,c(1,t_end)]), method = "manhattan")),
                # pc1_change = abs(sim.prcomp$rotation[1,"PC1"] - sim.prcomp$rotation[t_end,"PC1"]),
                # pc1_propvar = summary(sim.prcomp)$sdev[1]^2 / sum(summary(sim.prcomp)$sdev^2),
                # pc_totvar = sum(summary(sim.prcomp)$sdev^2)
  )
  
  
  return(Res)
}

#' Simulate all extinctions
#' 
#' For a given simulation and parameters, simulate the extinction of all
#' surviving species for the same time as the original simulation
#'
#' @param sim A matrix with species as rows and timepoints as columns with
#' the results of a simulation
#' @param params  list with starting conditions (x0), growth rates (mu), and
#' interaction matrix (M) for gLV simulation
#' @param timeout Number of seconds to wait before stopping simulation 
#'
#' @return A tibble with the results of all the extinction simulations
#' 
#' @export
simulate_all_extinctions <- function(sim, params, timeout = 600){
  
  # sim <- Sims$sim[[1]]
  # timeout <- 600
  
  # Check if simulations are needed
  if(is(sim, "logical") && is.na(sim)){
    return(NULL)
  }
  
  n_t <- dim(sim)[2]
  
  # Identify surviving species
  x_t <- sim[,n_t]
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
    
    # Remove extra species
    # Ensure simulation works when only one species is left after
    # extinction
    ii <- x_e > 0
    x_e <- x_e[ii]
    mu_e <- mu_e[ii]
    M_e <- M_e[ii,ii, drop = FALSE]
    
    # New M params
    # n_species_e <- nrow(M_e)
    # ii_diag <- diag(n_species_e) == 1
    # n_int <- n_species_e ^ 2 - n_species_e
    # p_noint_e <- sum(M_e[!ii_diag] == 0) / n_int
    # p_neg_e <- sum((M_e[!ii_diag] != 0) & (M_e[!ii_diag] < 0)) / (n_int - sum(M_e[!ii_diag] == 0))
    
    # Update parameters
    params_e <- list(x0 = x_e, mu = mu_e, M = M_e)
    
    # Simulate time after extinction
    sim_e <- sim_glv(params = params_e, n_t = n_t, timeout = timeout)
    sim_e$extinct_species <- as.numeric(spec)
    sim_e$spec_abun <- as.numeric(x_t[spec])
    sim_e$spec_freq <- as.numeric(sim_e$spec_abun / sum(x_t))
    
    Ext <- bind_rows(Ext, sim_e)
    sim_e <- NULL
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
hyper
Tab <- read_tsv(args$sims) %>%
  select(simid=id, n_species, p_noint, p_neg, seed) %>% # Renamming unnecesary
  pmap_dfr(.f = function(simid, n_species, p_noint, p_neg, seed){
    
    i <- 1
    n_species <- hyper$n_species[i]
    p_noint <- hyper$p_noint[i]
    p_neg <- hyper$p_neg[i]
    seed <- hyper$seed[i]
    simid <- hyper$id[i]
    
    message(paste0("===== Simulation ", simid," ====="))
    n_t <- 1000
    
    # Generate parameters and simulate
    set.seed(seed)
    params <- generate_params(n_species = n_species,
                              p_noint = p_noint, 
                              p_neg = p_neg,
                              ks_boost = args$ks_boost)
    
    # print(params)
    Sims <- sim_glv(params = params, n_t = n_t, timeout = 600)
    Sims$id <- simid
    
    # Sims
    
    if( !is.logical(Sims$sim[[1]]) ){
      if(Sims$richness_end > 1){
        # Simulate extinctions
        Ext <- simulate_all_extinctions(sim = Sims$sim[[1]], params = params,
                                        timeout = 600)
        if(!is.null(Ext)){
          Ext <- Ext %>%
            arrange(extinct_species) %>%
            mutate(ext_ii = 1:length(extinct_species)) %>%
            mutate(id = simid) 
        }
        Sims <- bind_rows(Sims, Ext)
        
        # Extract & write networks
        x_t <- Sims$sim[[1]][,n_t]
        surv_specs <- which(x_t > 0)
        M_surv <- params$M[surv_specs, surv_specs]
        write_tsv(M_surv %>% as_tibble, paste0(simid, "_net.tsv"), col_names = FALSE)
        
        # Extract metrics of species importance that we want to predict later
        Ext %>%
          transmute(prop_extinction = (richness_start - richness_end) / (length(surv_specs) - 1),
                    bray_change,
                    ks_bray = bray_change * (1 - spec_freq),
                    extinct_species) %>%
          mutate(ks_label = 1*(extinct_species == params$ks_ii)) %>%
          select(-extinct_species) %>%
          write_tsv(paste0(simid, "_target_obs.tsv"))
        
        
        # Extract node level features
        # help(package = "igraph")
        M_surv_graph <- igraph::graph_from_adjacency_matrix(t(M_surv), mode = "directed", weighted = TRUE, diag = FALSE)
        tibble(deg_in = igraph::degree(M_surv_graph, mode = "in"),
               deg_out = igraph::degree(M_surv_graph, mode = "out"),
               strength_in = igraph::strength(M_surv_graph, mode = "in"),
               strength_out = igraph::strength(M_surv_graph, mode = "out"),
               
               betweenness = igraph::betweenness(M_surv_graph, weights = abs(igraph::E(M_surv_graph)$weight)),
               closeness = igraph::closeness(M_surv_graph, , weights = abs(igraph::E(M_surv_graph)$weight)),
               harmonic_centrality = igraph::harmonic_centrality(M_surv_graph, weights = abs(igraph::E(M_surv_graph)$weight)),
               eigen_centrality = igraph::eigen_centrality(M_surv_graph)$vector,
               page_rank = igraph::page_rank(M_surv_graph, weights = abs(igraph::E(M_surv_graph)$weight))$vector, # Can't use negative values
               authority_score = igraph::authority_score(M_surv_graph)$vector,
               hub_score = igraph::hub_score(M_surv_graph)$vector,
               
               transitivity = igraph::transitivity(M_surv_graph, type = "local"),
               local_efficiency = igraph::local_efficiency(M_surv_graph, weights = abs(igraph::E(M_surv_graph)$weight)),
               knn = igraph::knn(M_surv_graph)$knn,                        # average neighbor degree,
               eccentricity = igraph::eccentricity(M_surv_graph, weights = abs(igraph::E(M_surv_graph)$weight)),
               # walktrap = igraph::membership(igraph::cluster_walktrap(M_surv_graph, weights = abs(igraph::E(M_surv_graph)$weight)))
               triangles = igraph::count_triangles(M_surv_graph)
              
               )
        

        
        

        
        
        

        
        
      }else{
        message("\t>Only one species survived, no additional extinctions.")
      }
    } else{
      message("\t>Main simulation failed. Extinctions cannot be tested.")
    }
    
    Sims
    
    # Save full results and table
    save(Sims, file = file.path(args$rdats, paste0(simid, ".rdat") ))
    Tab <- Sims %>%
      select(-params, -sim) 
    Tab %>%
      write_tsv(file = file.path(args$tsvs, paste0(simid, ".tsv")) )
    
    return(Tab)
  })

write_tsv(Tab, file.path(args$outdir, "summary.tsv"))

date()
sessionInfo()

