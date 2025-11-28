
# Load data
import pandas as pd
import numpy as np
import torch 
import os
from datetime import datetime

# Section: Generate-paths
exp_dir = "/mnt/data/sur/users/mrivera/Controls/exp_20251125"
A_dir = os.path.join(exp_dir, "A-mat")
tgt_dir = os.path.join(exp_dir, "GNN-targets")
odes_path = os.path.join(exp_dir, "raw-ODEs")
data_path = os.path.join(exp_dir, "Simulation-parameters.tsv")

# Generate ID for training.
timeID = datetime.now().strftime("Y%YM%mD%d")

#  Load-data
data_ids = pd.read_csv(data_path, sep="\t", usecols=['id'])['id']

#===============================================
# Loading functions
# SECTION: Load-function
from torch_geometric.data import Data
import pyarrow.feather as feather

def load_single_data(id, A_dir, tgt_dir):
    # Load adjacency matrix 
    A_path = os.path.join(A_dir, f"A_{id}.feather")
    A = pd.read_feather(A_path).to_numpy(dtype=np.float32)
    # Vector of edge weights
    row_idx, col_idx = np.nonzero(A)
    edge_weights = A[row_idx, col_idx]
    # Convert to torch tensors efficiently
    edge_index = torch.from_numpy(np.vstack([row_idx, col_idx]).astype(np.int64))
    edge_weights = torch.from_numpy(edge_weights)
    # Section: Load target features 
    tgt_path = os.path.join(tgt_dir, f"tgt_{id}.feather")
    tgt_table =  feather.read_table(tgt_path, columns=['K_s'])
    # columns_to_load = ['new_ext', 'BC_diss', 'K_s']  # all columns except K_s
    #tgt_table =  feather.read_table(tgt_path, columns=columns_to_load)
    y_tensor = torch.from_numpy(tgt_table.to_pandas().to_numpy(dtype=np.float32))   
    # Node features - simple ones vector
    n = A.shape[0]
    x_tensor = torch.ones(n, 1, dtype=torch.float32)
    # Clean up large intermediate
    del A, tgt_table
    # Create Data object
    data = Data(
        x=x_tensor,
        edge_weights=edge_weights,
        edge_index=edge_index,
        y=y_tensor
    )
    return data

# Testing line
id = data_ids.iloc[0]
# load_single_data(data_ids.iloc[0], A_dir, tgt_dir)
#===============================================
# SECTION: Divide-data
import random

# Load all data samples (for demo, we use only first 100 samples)
indices = list(range(1, len(data_ids)))  # Indices 1-100
random.shuffle(indices)  # Uses Python's random module (already seeded)

# Now select first 80 for training, rest for validation
indx = round(len(indices) * .8)
train_indices = indices[:indx]            # First 80 shuffled indices
val_indices = indices[indx:]              # Last 20 shuffled indices

#===============================================
# Parallelize function
from concurrent.futures import ThreadPoolExecutor
import time
from tqdm import tqdm

# Generate data and separate in batches
path = f'/mnt/data/sur/users/mrivera/Cuda-tensors/{os.path.basename(exp_dir)}'
os.makedirs(path, exist_ok=True)

def generate_data_parallel(idx, A_dir, tgt_dir, num_workers=4):  # idx is a list
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        data_list = list(executor.map(
            load_single_data,           
            idx,                       # List of IDs to iterate over
            [A_dir]*len(idx),          # Repeat A_dir for each ID
            [tgt_dir]*len(idx)         # Repeat tgt_dir for each ID
        ))
    return data_list


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

batching(all_ids = data_ids, indexes = train_indices, A_dir = A_dir, tgt_dir = tgt_dir, path = path, batch_size=100, num_workers=6, prefix='TrainBatch')
batching(all_ids = data_ids, indexes = val_indices, A_dir = A_dir, tgt_dir = tgt_dir, path = path, batch_size=100, num_workers=6, prefix='ValBatch')

#===============================================
# Create Zip of batches
import shutil

# Zip entire directory
name = os.path.basename(exp_dir)
shutil.make_archive(f'/mnt/data/sur/users/mrivera/Cuda-tensors/{name}', 'zip', 
                    '/mnt/data/sur/users/mrivera/Cuda-tensors/exp_20251125')

