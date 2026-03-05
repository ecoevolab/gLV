"""
Filtered nodes Data Generation Script for GNN model.
=================================================
Purpose:
    Generate tensors required for the GNN model.
    Species with low relative abundance at perturbation time are filtered out prior to tensor construction;
    feature dimensionality is preserved across the remaining nodes.

Data Format:
    - X : Node feature matrix. Each row corresponds to its respective node statistics in the network.
    - Y : Target matrix. Each row contains the metric(s) of interest for the corresponding node.
    - edge_index : Graph connectivity in COO format, encoding which species are connected to which.
    - edge_attr  : Edge weights representing the strength of interaction.

Dependencies:
    torch==2.8, pandas==2.3.3, numpy==2.0.2

Author: Manuel Rivera
Date:   4-March-2026
"""

#-------------------------------------------------------------
# Section: Import data
# Load data
import pandas as pd
import numpy as np
import torch 
import os
from datetime import datetime

# Section: Generate-paths
experiment_dir = '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone'
interactions_dir = os.path.join(experiment_dir, "Interactions")
targets_dir = os.path.join(experiment_dir, "ExtSummaries")
topologies_dir = os.path.join(experiment_dir, "Topologies")
data_path = os.path.join(experiment_dir, "simulation-params.tsv")

# Generate ID for training.
timeID = datetime.now().strftime("Y%YM%mD%d")

#  Load-data
data_rows = pd.read_csv(data_path, sep="\t", usecols=['id','key'])

#-------------------------------------------------------------
# SECTION: Generate function to load data for one simulation.
from torch_geometric.data import Data
import pyarrow.feather as feather

def load_single_data(row, interactions_dir, targets_dir, topologies_dir):
    #------------------------------
    # Section: Load adjacency matrix 
    id = row['id']
    keystone = row['key'] - 1 # Convert to python 0 index
    #-------------------------------
    # Section: Read target features
    tgt_path = os.path.join(targets_dir, f"ExtSummary_{id}.feather")
    cols = ['rel_pop_initial', 'keystoneness']
    tgt_table =  feather.read_table(tgt_path, columns=cols).to_pandas()
    # Filter extinct species
    alive_vector = tgt_table['rel_pop_initial'] > 1e-06     # vector
    alive_indices = np.where(alive_vector)[0]                           # indices
    # Subset indices from table
    y_values = tgt_table['keystoneness'].iloc[alive_indices]
    y_tensor = torch.from_numpy(y_values.to_numpy(dtype=np.float32))  
    #------------------------------
    # Section: Load adjacency matrix 
    A_path = os.path.join(interactions_dir, f"A_{id}.feather")
    A = feather.read_table(A_path).to_pandas()
    # Subset interactions of alive population
    A_subset = A.iloc[alive_indices, alive_indices].to_numpy()
    # Vector of edge weights
    row_idx, col_idx = np.nonzero(A_subset) # who with whom
    edge_weights = A_subset[row_idx, col_idx] # edge weights 
    # Convert to torch tensors efficiently
    edge_index_tensor = torch.from_numpy(np.vstack([row_idx, col_idx]).astype(np.int64))
    edge_weights_tensor = torch.from_numpy(edge_weights)   
    #------------------------------ 
    # Section: Load node features 
    topology_path = os.path.join(topologies_dir, f"Topology_{id}.feather")
    network_stats = feather.read_feather(topology_path).iloc[alive_indices]
    network_stats = network_stats.drop(columns=["total_degree"])
    x_tensor = torch.from_numpy(network_stats.to_numpy(dtype=np.float32))  
    # Add dummy data, vector of ones
    n = x_tensor.shape[0]
    x_tensor = torch.cat([x_tensor, torch.ones(n, 1, dtype=torch.float32)], dim=1)
    #------------------------------
    # Clean up large intermediate
    del A, tgt_table, A_subset
    #------------------------------
    # Create Data object
    data = Data(
        x=x_tensor,
        edge_weights=edge_weights_tensor,
        edge_index=edge_index_tensor,
        y=y_tensor
    )
    return data

# Sanity check: load one data sample
row = data_rows.iloc[0]
eg = load_single_data(row, interactions_dir, targets_dir, topologies_dir)

#-------------------------------------------------------------
# Section: Parallelize function
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from tqdm import tqdm  # optional, for progress bar

def load_dataset_parallel(df, interactions_dir, targets_dir, topologies_dir, max_workers=8):
    # Fix arguments
    func = partial(load_single_data, interactions_dir=interactions_dir, targets_dir=targets_dir, topologies_dir=topologies_dir)
    # Generate list of rows
    rows = [row for _, row in df.iterrows()]
    # Preallocate results
    results = [None] * len(rows)  
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_idx = {executor.submit(func, row): i for i, row in enumerate(rows)}
        for future in tqdm(as_completed(future_to_idx), total=len(rows)):
            idx = future_to_idx[future]
            try:
                results[idx] = future.result()
            except Exception as e:
                print(f"Row {idx} failed: {e}")
                results[idx] = None
    return results


#-------------------------------------------------------------
# Section: Create batches
import time

def batching(df, A_dir, targets_dir, topology_dir, path, batch_size=250, num_workers=6, prefix='TrainBatch'):
    #----------------------------
    # Section: Calculate batch size
    os.makedirs(path, exist_ok=True) # Generate directory if it does not exist
    start = time.time()
    rows = df.shape[0]
    num_batches = rows // batch_size
    #------------------------
    # Section: Process rows
    for i in tqdm(range(num_batches), desc="Processing batches"):
        batch_rows = df[i*batch_size : (i+1)*batch_size]  # slice the DataFrame
        batch_data = load_dataset_parallel(batch_rows, A_dir, targets_dir, topology_dir, num_workers)  # different name
        name = os.path.join(path, f'{prefix}_{i}.pt')
        torch.save(batch_data, name)
        del batch_data  # free memory
        print(f"Saved batch {i}/{num_batches-1}: {name}")
    #------------------------
    # Handle leftover rows that don't fill a complete batch
    remainder = rows % batch_size
    if remainder > 0:
        batch_rows = df[num_batches*batch_size:]
        batch_data = load_dataset_parallel(batch_rows, A_dir, targets_dir, topology_dir, num_workers)
        name = os.path.join(path, f'{prefix}_{num_batches}.pt')
        torch.save(batch_data, name)
        del batch_data
        print(f"Saved remainder batch: {name}")
    #--------------
    elapsed = time.time() - start
    print(f">> Batch size: {batch_size} | Num batches: {num_batches} | Elapsed: {elapsed:.2f}s")

#-----------------------
# Run parallelization
tensors_path = '/mnt/data/sur/users/mrivera/Cuda-tensors/Boosted_filtered'
batching(df = data_rows,
         A_dir = interactions_dir, 
         targets_dir = targets_dir,
         topology_dir = topologies_dir, 
         path = tensors_path, 
         batch_size=100, num_workers=6, prefix='TrainBatch_')



#===============================================
# Create Zip of batches
import shutil

# Zip entire directory
name = os.path.basename(tensors_path)
zip_dir = f'/mnt/data/sur/users/mrivera/Cuda-tensors/{name}'
shutil.make_archive(zip_dir, 'zip', tensors_path)

