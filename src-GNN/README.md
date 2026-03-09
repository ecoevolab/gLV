# src-GNN

## Description

This directory contains the scripts to convert extinction summaries into tensors for **PyTorch**,
and to define, train, and evaluate **Graph Neural Network (GNN)** models on the resulting data.

## Structure
```
src-GNN/
├── Batching-tensors.py
├── filter_nodes_tensors.py
├── all_features_runner.py
└── dummy_runner.py
```

## Scripts

| Script | Description |
|---|---|
| `Batching-tensors.py` | Converts extinction summaries into tensors for all network nodes. Parallelized for performance |
| `filter_nodes_tensors.py` | Generates tensors excluding species (nodes) whose relative abundance falls below a set threshold. Parallelized for performance |
| `all_features_runner.py` | Declares, trains, and evaluates the GNN using all node features. Outputs: predicted vs. true values, node with highest vs. predicted maximum value, and loss over time |
| `dummy_runner.py` | Same pipeline as `all_features_runner.py` but uses dummy data as node features, for testing/debugging purposes |
