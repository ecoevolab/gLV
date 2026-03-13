"""
Tensors Generation Script for GNN model.
=================================================
Purpose:
    Generate tensors required for the GNN model.
    Species with low relative abundance at perturbation time are filtered out already.

Data Format:
    - X : Node feature matrix. Each row corresponds to its respective node statistics in the network.
    - Y : Target matrix. Each row contains the metric(s) of interest for the corresponding node.
    - edge_index : Graph connectivity in COO format, encoding which species are connected to which.
    - edge_attr  : Edge weights representing the strength of interaction.

Dependencies:
    torch==2.8, pandas==2.3.3, numpy==2.0.2

Author: Manuel Rivera
Date:   13-March-2026
"""

#-------------------------------------------------------------
# Section: Import data
# Load data
import torch 
import pandas as pd
import numpy as np
import os
from datetime import datetime

# Section: Generate-paths
experiment_dir = '/mnt/data/sur/users/mrivera/Controls/Boosted_keystone'
filtered_inters = os.path.join(experiment_dir, "Filter_Interactions")
filtered_exts = os.path.join(experiment_dir, "Filter_ExtSummaries")
filtered_stats = os.path.join(experiment_dir, "Filter_Topologies")

# Generate ID for training.
timeID = datetime.now().strftime("Y%YM%mD%d")

#  Load-data
data_path = os.path.join(experiment_dir, "simulation-params.tsv")
data_rows = pd.read_csv(data_path, sep="\t", usecols=['id','key'])

import pyarrow.feather as feather
import pyarrow.compute as pc
success_path = os.path.join(experiment_dir, "failed_filters.feather")
success_ids = feather.read_table(success_path)
success_ids = success_ids.filter(pc.field('n') > 1)['id'].to_pylist()


#-------------------------------------------------------------
# SECTION: Generate function to load data for one simulation.
from torch_geometric.data import Data
import pyarrow.feather as feather

def load_single_data(id, filtered_inters, filtered_exts, filtered_stats):
    #-------------------------------
    # Section: Read target features
    tgt_path = os.path.join(filtered_exts, f"ExtSummary_{id}.feather")
    cols = ['keystoneness']
    tgt_table =  feather.read_table(tgt_path, columns=cols).to_pandas()
    y_tensor = torch.from_numpy(tgt_table.to_numpy(dtype=np.float32)) # Convert to tensor
    #------------------------------
    # Section: Load adjacency matrix 
    A_path = os.path.join(filtered_inters, f"A_{id}.feather")
    A = feather.read_table(A_path).to_pandas()
    # Vector of edge weights
    row_idx, col_idx = np.nonzero(A) # who with whom
    edge_weights = A.to_numpy()[row_idx, col_idx]# edge weights 
    # Convert to torch tensors efficiently
    edge_index_tensor = torch.from_numpy(np.vstack([row_idx, col_idx]).astype(np.int64))
    edge_weights_tensor = torch.from_numpy(edge_weights).float()   
    #------------------------------ 
    # Section: Load node features 
    features_path = os.path.join(filtered_stats, f"topology_{id}.feather")
    network_stats = feather.read_feather(features_path)
    x_tensor = torch.from_numpy(network_stats.to_numpy(dtype=np.float32))  
    # Add dummy data, vector of ones
    n = x_tensor.shape[0]
    x_tensor = torch.cat([x_tensor, torch.ones(n, 1, dtype=torch.float32)], dim=1)
    #------------------------------
    # Clean up large intermediate
    del A, tgt_table
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
id = success_ids[0]
eg = load_single_data(id, filtered_inters, filtered_exts, filtered_stats)

#-------------------------------------------------------------
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from tqdm import tqdm

def load_all_data(ids, filtered_inters, filtered_exts, filtered_stats, max_workers=8):
    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(load_single_data, id, filtered_inters, filtered_exts, filtered_stats): id for id in ids
        }
        for future in tqdm(as_completed(futures), total=len(ids)):
            id = futures[future]
            try:
                results[id] = future.result()
            except Exception as e:
                print(f"Failed for id {id}: {e}")
    return results

#-------------------------------------------------------------
# Section: Create batches
import time

def batching(list_ids, A_dir, targets_dir, topology_dir, path, batch_size=250, num_workers=6, prefix='TrainBatch'):
    #----------------------------
    # Section: Calculate batch size
    os.makedirs(path, exist_ok=True) # Generate directory if it does not exist
    start = time.time()
    n_ids = len(list_ids)
    num_batches = n_ids // batch_size
    #------------------------
    # Section: Process ids
    for i in tqdm(range(num_batches), desc="Processing batches"):
        batch_ids = list_ids[i*batch_size : (i+1)*batch_size]  # slice the DataFrame
        batch_data = load_all_data(batch_ids, A_dir, targets_dir, topology_dir, num_workers)  # different name
        name = os.path.join(path, f'{prefix}_{i}.pt')
        torch.save(batch_data, name)
        del batch_data  # free memory
        print(f"Saved batch {i}/{num_batches-1}: {name}")
    #------------------------
    # Handle leftover rows that don't fill a complete batch
    remainder = n_ids % batch_size
    if remainder > 0:
        batch_ids = list_ids[num_batches*batch_size:]
        batch_data = load_all_data(batch_ids, A_dir, targets_dir, topology_dir, num_workers) 
        name = os.path.join(path, f'{prefix}_{num_batches}.pt')
        torch.save(batch_data, name)
        del batch_data
        print(f"Saved remainder batch: {name}")
    #--------------
    elapsed = time.time() - start
    print(f">> Batch size: {batch_size} | Num batches: {num_batches} | Elapsed: {elapsed:.2f}s")

#-----------------------
# Run parallelization
tensors_path = '/mnt/data/sur/users/mrivera/Cuda-tensors/Boost_Surv_Only'
batching(list_ids = success_ids,
         A_dir = filtered_inters, 
         targets_dir = filtered_exts,
         topology_dir = filtered_stats, 
         path = tensors_path, 
         batch_size=100, num_workers=6, prefix='TrainBatch')



#===============================================
# Create Zip of batches
import shutil

# Zip entire directory
name = os.path.basename(tensors_path)
zip_dir = f'/mnt/data/sur/users/mrivera/Cuda-tensors/{name}'
shutil.make_archive(zip_dir, 'zip', tensors_path)

