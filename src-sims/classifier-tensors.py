# This code is done for creating the batces of tensors required for training the GNN. 
# This code is used to generate data required for the classifier (classes)
# It is done in a parallelized way, and the batches are saved in a directory as .pt files. 
# Finally, the entire directory is zipped.

#-------------------------------------------------------------
# Section: Import data
# Load data
import pandas as pd
import numpy as np
import torch 
import os
from datetime import datetime

# Section: Generate-paths
experiment_dir = '/mnt/data/sur/users/mrivera/Controls/Cascade_keystone'
interactions_dir = os.path.join(experiment_dir, "Interactions")
targets_dir = os.path.join(experiment_dir, "ExtSummaries")
topologies_dir = os.path.join(experiment_dir, "Topologies")
data_path = os.path.join(experiment_dir, "simulation-params.tsv")

# Generate ID for training.
timeID = datetime.now().strftime("Y%YM%mD%d")

#  Load-data
data = pd.read_csv(data_path, sep="\t", usecols=['id', 'key'])

#-------------------------------------------------------------
# SECTION: Generate function to load data for one simulation.
from torch_geometric.data import Data
import pyarrow.feather as feather

def load_single_data(row, interactions_dir, topologies_dir):
    #------------------------------
    # Section: Load adjacency matrix 
    id = row['id']
    keystone = row['key'] - 1 # Convert to python 0 index
    A_path = os.path.join(interactions_dir, f"A_{id}.feather")
    A = pd.read_feather(A_path).to_numpy(dtype=np.float32)
    # Vector of edge weights
    row_idx, col_idx = np.nonzero(A) # who with whom
    edge_weights = A[row_idx, col_idx] # edge weights 
    # Convert to torch tensors efficiently
    edge_index = torch.from_numpy(np.vstack([row_idx, col_idx]).astype(np.int64))
    edge_weights = torch.from_numpy(edge_weights)
    #------------------------------
    # Section: Load target features 
    n = A.shape[0]
    y_tensor = torch.zeros(n,1) 
    y_tensor[keystone] = 1 
    #------------------------------ 
    # Section: Load node features 
    # Topological features
    topology_path = os.path.join(topologies_dir, f"Topology_{id}.feather")
    network_stats = feather.read_feather(topology_path)
    network_stats = network_stats.drop(columns=["total_degree"])
    x_tmp = torch.from_numpy(network_stats.to_numpy(dtype=np.float32))  
    # Add dummy data, vector of ones
    dummy_tensor = torch.ones(n, 1, dtype=torch.float32)
    x_tensor = torch.cat([x_tmp, dummy_tensor], dim=1)
    #------------------------------
    # Clean up large intermediate
    del A
    #------------------------------
    # Create Data object
    data = Data(
        x=x_tensor,
        edge_weights=edge_weights,
        edge_index=edge_index,
        y=y_tensor
    )
    return data

row = data.iloc[0]
x = load_single_data(row, interactions_dir, topologies_dir)


#-------------------------------------------------------------
# Section: Parallelize function
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from tqdm import tqdm  # optional, for progress bar

def load_dataset_parallel(df, interactions_dir, topologies_dir, max_workers=8):
    # Fix arguments
    func = partial(load_single_data, interactions_dir=interactions_dir, topologies_dir=topologies_dir)
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

def batching(df, A_dir, topology_dir, path, batch_size=250, num_workers=6, prefix='TrainBatch'):
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
        batch_data = load_dataset_parallel(batch_rows, A_dir, topology_dir, num_workers)  # different name
        name = os.path.join(path, f'{prefix}_{i}.pt')
        torch.save(batch_data, name)
        del batch_data  # free memory
        print(f"Saved batch {i}/{num_batches-1}: {name}")
    #------------------------
    # Handle leftover rows that don't fill a complete batch
    remainder = rows % batch_size
    if remainder > 0:
        batch_rows = df[num_batches*batch_size:]
        batch_data = load_dataset_parallel(batch_rows, A_dir, topology_dir, num_workers)
        name = os.path.join(path, f'{prefix}_{num_batches}.pt')
        torch.save(batch_data, name)
        del batch_data
        print(f"Saved remainder batch: {name}")
    #--------------
    elapsed = time.time() - start
    print(f">> Batch size: {batch_size} | Num batches: {num_batches} | Elapsed: {elapsed:.2f}s")

#-----------------------
# Run parallelization
tensors_path = '/mnt/data/sur/users/mrivera/Cuda-tensors/Cascade_classes'
batching(df = data,
         A_dir = interactions_dir, 
         topology_dir = topologies_dir, 
         path = tensors_path, 
         batch_size=100, num_workers=6, prefix='TrainBatch_')

#===============================================
# Create Zip of batches
import shutil

# Zip entire directory
name = os.path.basename(experiment_dir)
zip_dir = f'/mnt/data/sur/users/mrivera/Cuda-tensors/{name}'
shutil.make_archive(zip_dir, 'zip', tensors_path)