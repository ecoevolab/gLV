# This code is done for creating the batces of tensors required for training the GNN. 
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
experiment_dir = "/mnt/data/sur/users/mrivera/Controls/NegCtrl-V1"
interactions_dir = os.path.join(experiment_dir, "Interactions")
targets_dir = os.path.join(experiment_dir, "ExtSummaries")
topologies_dir = os.path.join(experiment_dir, "Topologies")
data_path = os.path.join(experiment_dir, "simulation-params.tsv")

# Generate ID for training.
timeID = datetime.now().strftime("Y%YM%mD%d")

#  Load-data
data_ids = pd.read_csv(data_path, sep="\t", usecols=['id'])['id']

#-------------------------------------------------------------
# SECTION: Generate function to load data for one simulation.
from torch_geometric.data import Data
import pyarrow.feather as feather
id = data_ids.iloc[0]

def load_single_data(id, interactions_dir, targets_dir, topologies_dir):
    #------------------------------
    # Section: Load adjacency matrix 
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
    tgt_path = os.path.join(targets_dir, f"ExtSummary_{id}.feather")
    cols = ['prop_extinctions', 'dissimilarity_bc', 'keystoneness']
    tgt_table =  feather.read_table(tgt_path, columns=cols)
    y_tensor = torch.from_numpy(tgt_table.to_pandas().to_numpy(dtype=np.float32))   
    # Convert to 0 nan if any
    y_tensor[torch.isnan(y_tensor)] = 0.0
    #------------------------------ 
    # Section: Load node features 
    # Topological features
    topology_path = os.path.join(topologies_dir, f"Topology_{id}.feather")
    network_stats = feather.read_feather(topology_path)
    network_stats = network_stats.drop(columns=["total_degree"])
    x_tmp = torch.from_numpy(network_stats.to_numpy(dtype=np.float32))  
    # Add dummy data, vector of ones
    n = A.shape[0]
    dummy_tensor = torch.ones(n, 1, dtype=torch.float32)
    x_tensor = torch.cat([x_tmp, dummy_tensor], dim=1)
    #------------------------------
    # Clean up large intermediate
    del A, tgt_table
    #------------------------------
    # Create Data object
    data = Data(
        x=x_tensor,
        edge_weights=edge_weights,
        edge_index=edge_index,
        y=y_tensor
    )
    return data

# Sanity check: load one data sample
load_single_data(data_ids.iloc[0], interactions_dir, targets_dir, topologies_dir)

#-------------------------------------------------------------
# Section: Create batches of data
# We create batches of data to speed up training.
# Batches with more than 500 samples can cause memory issues, so we will use a batch size of 100.
import random

# Load all data samples (for demo, we use only first 100 samples)
indices = list(range(1, len(data_ids)))  # Indices 1-100
random.shuffle(indices)  # Uses Python's random module (already seeded)

# We will use all data for training and some for validation.

# Now select first 80 for training, rest for validation
# indx = round(len(indices) * .8)
# train_indices = indices[:indx]            # First 80 shuffled indices
# val_indices = indices[indx:]              # Last 20 shuffled indices

#-------------------------------------------------------------
# Section: Parallelize function
from concurrent.futures import ThreadPoolExecutor
import time
from tqdm import tqdm

# Generate data and separate in batches
tensors_path = f'/mnt/data/sur/users/mrivera/Cuda-tensors/{os.path.basename(experiment_dir)}'
os.makedirs(tensors_path, exist_ok=True)

def generate_data_parallel(idx, A_dir, tgt_dir, num_workers=4):  # idx is a list
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        data_list = list(executor.map(
            load_single_data,           
            idx,                       # List of IDs to iterate over
            [A_dir]*len(idx),          # Repeat A_dir for each ID
            [tgt_dir]*len(idx),         # Repeat tgt_dir for each ID
            [topologies_dir]*len(idx)  # Repeat topologies_dir for each ID
        ))
    return data_list

#-------------------------------------------------------------
# Section: Batching function
def batching(all_ids, indexes, A_dir, tgt_dir, path, batch_size=250, num_workers=6, prefix='TrainBatch'):
    start = time.time()
    batching_ids = [all_ids[i] for i in indexes]
    num_batches = (len(batching_ids) + batch_size - 1) // batch_size
    for i in tqdm(range(num_batches), desc="Processing batches"):
        batch_ids = batching_ids[i*batch_size : (i+1)*batch_size]
        data = generate_data_parallel(batch_ids, A_dir, tgt_dir, num_workers)
        name = os.path.join(path, f'{prefix}_{i}.pt')
        torch.save(data, name)  # Actually save the data
        del data  # Free memory immediately
        print(f"Saved batch {i}/{num_batches-1}: {name}")
    elapsed = time.time() - start
    print(f">> The batch size is of {batch_size}, the number of batches was {num_batches}, while the elapsed time is of: {elapsed:.2f}")
    return None

batching(all_ids = data_ids, 
         indexes = indices, A_dir = interactions_dir, 
         tgt_dir = targets_dir, path = tensors_path, 
         batch_size=100, num_workers=6, prefix='TrainBatch')

# batching(all_ids = data_ids, indexes = val_indices, A_dir = A_dir, tgt_dir = tgt_dir, path = path, batch_size=100, num_workers=6, prefix='ValBatch')

#===============================================
# Create Zip of batches
import shutil

# Zip entire directory
name = os.path.basename(experiment_dir)
zip_dir = f'/mnt/data/sur/users/mrivera/Cuda-tensors/{name}'
shutil.make_archive(zip_dir, 'zip', tensors_path)

