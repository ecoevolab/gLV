# src-sims

## Description

This directory contains the functions and scripts to generate **control** and **training** datasets
of in silico microbial communities using the **generalized Lotka-Volterra (gLV)** model.

ODEs are solved with `ode45` from the [`deSolve`](https://cran.r-project.org/package=deSolve) package,
and simulations are parallelized across cores via `mclapply` from the 
[`parallel`](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) package.

## Structure
```
src-sims/
└── FUN/
    ├── generate-boosted-params.R
    ├── generate-cascade-params.R
    ├── generate-controls-params.R
    ├── solver-gLV.R
    ├── forge-symls.R
    └── extinctions-fun.R
└── controls_runner.R

```

### `FUN/`

Parameter generation, ODE solving, and post-processing utilities for gLV simulations.

| Script | Description |
|---|---|
| `generate-boosted-params.R` | Generates gLV parameters by boosting the full column of a keystone species (excluding its diagonal) |
| `generate-cascade-params.R` | Generates gLV parameters with cascading interaction boosts (e.g. A→B→C, where A↔B and B↔C are boosted by factor *k*) |
| `generate-controls-params.R` | Generates gLV parameters where one column is one order of magnitude weaker than the rest — no boosting applied |
| `solver-gLV.R` | ODE solver for the gLV equations |
| `forge-symls.R` | Merges per-core output directories into a single dataset, preventing race conditions during parallelization |
| `extinctions-fun.R` | Performs systematic species extinctions and computes each species' removal impact on the community |


