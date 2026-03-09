# gLV In Silico Community Simulator

> ⚠️ This project is under active development. 

This repository provides tools for generating generalized Lotka-Volterra (gLV) 
in silico communities through simulation, and training Graph Neural 
Networks (GNNs) to model their dynamics.

## Overview

The generalized Lotka-Volterra model describes the growth and interactions of 
species in a community. This pipeline simulates synthetic communities, 
structures the resulting data as tensors, and trains GNN models to learn from 
those ecological dynamics.

## Repository Structure

| Directory | Description |
|---|---|
| [`src-sims/`](src-sims/README.md) | Scripts and functions for simulating control and training datasets |
| [`src-GNN/`](src-GNN/README.md) | Tensor data generation, GNN model definitions, and training/testing scripts |
| `Markdowns/` | Notes and records of experiments performed |

