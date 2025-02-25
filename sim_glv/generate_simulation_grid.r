# (C) Copyright 2025 Sur Herrera Paredes
# 
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

#' Prepare combinations of parameters for simulations. This is not
#' a script to be re-run

#' Here we prepare a grid of simulations with gLV

Sims <- expand_grid(n_species = rep(c(20, 100), 10),
                    p_noint = seq(from = .1, to = .9, by = 0.1),
                    p_neg = seq(from = .1, to = .9, by = 0.1)) %>%
  mutate(id = ids::random_id(n = length(n_species), bytes = 10)) %>%
  mutate(seed = as.vector(random::randomNumbers(n = length(n_species),
                                      min = 1,
                                      max = 1e6,
                                      col = 1)))
Sims
write_tsv(Sims, "glv_simulation_parameters.tsv")

#' Split in chunks of same number of `n_species` and same `p_noint`
#' This might be useful if we need to parallelize in cluster.
Sims %>%
  group_by(n_species, p_noint) %>%
  mutate(group_id = cur_group_id()) %>%
  group_walk(~ write_tsv(.x %>% 
                           bind_cols(n_species = .y$n_species,
                                     p_noint = .y$p_noint), 
                         paste0("simulations_batch", 
                                unique(.x$group_id), 
                                ".tsv")))

