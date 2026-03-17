"""
Hybrid data generation for GNN.
=================================================
Purpose:
    Generate tensors required for the GNN model.
    Species with low relative abundance at perturbation time are filtered out already.
    Survival nodes statistics are taken out from the whole network.

Important:
    Node features include the network statistics calculated using the whole network.
    Targets are the extinctions performed in the subcommunities of the model.

Data Format:
    - X : Node feature matrix. Each row corresponds to its respective node statistics in the network.
    - Y : Target matrix. Each row contains the metric(s) of interest for the corresponding node.
    - edge_index : Graph connectivity in COO format, encoding which species are connected to which.
    - edge_attr  : Edge weights representing the strength of interaction.

Dependencies:
    torch==2.8, pandas==2.3.3, numpy==2.0.2

Author: Manuel Rivera
Date:   16-March-2026
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
experiment_dir = '/mnt/data/sur/users/mrivera/Controls/KBoost_dataset_v2'
networks_dir = os.path.join(experiment_dir, "Interactions")
filtered_tgt_dir = os.path.join(experiment_dir, "Filtered_ExtSummaries")
features_dir = os.path.join(experiment_dir, "Topologies")
outs_dir = os.path.join(experiment_dir, "RawOutputs")

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


# What files are good?
import re

pass_ids = [re.search(r"ExtSummary_(.+)\.feather", f).group(1) for f in os.listdir(filtered_tgt_dir)]
id = pass_ids[0]

"""
We did not save the indices of which species should be filtered, so we will calculate them 
from the last column of the output. The last column corresponds to t=1000, the time at 
which the community is expected to be in quasi-stable state. NumPy's np.where() already 
returns 0-based Python indices, so no conversion is needed.

Node features include the network statistics calculated using the whole network.
Targets are the extinctions performed in the subcommunities of the model.
"""
# SECTION: Generate function to load data for one simulation.
from torch_geometric.data import Data
import pyarrow.feather as feather

def load_single_data(id, output_dir, target_dir, networks_dir, features_dir):
    #-------------------------------
    # Section: Load which species to filter
    output_path = os.path.join(output_dir, f"RawOutput_{id}.feather")
    output = feather.read_table(output_path).to_pandas().iloc[:, 20].to_numpy()
    relative = output/sum(output)
    to_filter = np.where(relative > 1e-06)[0]
    #-------------------------------
    # Section: Read target features
    tgt_path = os.path.join(target_dir, f"ExtSummary_{id}.feather")
    tgt_table =  feather.read_table(tgt_path, columns = ['keystoneness']).to_pandas()
    y_tensor = torch.from_numpy(tgt_table.to_numpy(dtype=np.float32)) # Convert to tensor
    #------------------------------
    # Section: Load adjacency matrix 
    A_path = os.path.join(networks_dir, f"A_{id}.feather")
    A = feather.read_table(A_path).to_pandas().to_numpy()
    A_filter = A[np.ix_(to_filter, to_filter)]
    # Vector of edge weights
    row_idx, col_idx = np.nonzero(A_filter) # who with whom
    edge_weights = A_filter[row_idx, col_idx]# edge weights 
    # Convert to torch tensors efficiently
    edge_index_tensor = torch.from_numpy(np.vstack([row_idx, col_idx]).astype(np.int64))
    edge_weights_tensor = torch.from_numpy(edge_weights).float()   
    #------------------------------ 
    # Section: Load node features 
    features_path = os.path.join(features_dir, f"Topology_{id}.feather")
    network_stats = feather.read_feather(features_path).iloc[to_filter]
    x_tensor = torch.from_numpy(network_stats.to_numpy(dtype=np.float32))  
    # Add dummy data, vector of ones
    n = x_tensor.shape[0]
    x_tensor = torch.cat([x_tensor, torch.ones(n, 1, dtype=torch.float32)], dim=1)
    #------------------------------
    # Clean up large intermediate
    del A, A_filter, tgt_table
    #------------------------------
    # Create Data object
    data = Data(
        x=x_tensor,
        edge_weights=edge_weights_tensor,
        edge_index=edge_index_tensor,
        y=y_tensor
    )
    return data

load_single_data(id = id, output_dir = outs_dir, target_dir = filtered_tgt_dir, networks_dir = networks_dir, features_dir = features_dir)

"""
Once we have our function to load data, we parallelize it.
"""
# Section: Parallelize function
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from tqdm import tqdm  # optional, for progress bar

def load_dataset_parallel(ids_list, outs_dir, target_dir, networks_dir, features_dir, max_workers=8):
    # Fix arguments
    func = partial(load_single_data, output_dir = outs_dir, target_dir = target_dir, networks_dir = networks_dir, features_dir = features_dir)
    # Generate list of ids
    ids_list = [id for id in pass_ids]
    # Preallocate results
    results = [None] * len(ids_list)  
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_idx = {executor.submit(func, id): i for i, id in enumerate(ids_list)}
        for future in tqdm(as_completed(future_to_idx), total=len(ids_list)):
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

def batching(ids_list, outputs_dir, targets_dir, networks_dir, features_dir, save_path, batch_size=250, num_workers=6, prefix='TrainBatch'):
    #----------------------------
    # Section: Calculate batch size
    os.makedirs(save_path, exist_ok=True) # Generate directory if it does not exist
    start = time.time()
    n_ids = len(ids_list)
    num_batches = n_ids // batch_size
    #------------------------
    # Section: Process rows
    for i in tqdm(range(num_batches), desc="Processing batches"):
        batch_ids = ids_list[i*batch_size : (i+1)*batch_size]  
        batch_data = load_dataset_parallel(batch_ids, outputs_dir, targets_dir, networks_dir, features_dir, num_workers)  
        name = os.path.join(save_path, f'{prefix}_{i}.pt')
        torch.save(batch_data, name)
        del batch_data  # free memory
        print(f"Saved batch {i}/{num_batches-1}: {name}")
    #------------------------
    # Handle leftover rows that don't fill a complete batch
    remainder = n_ids % batch_size
    if remainder > 0:
        batch_ids = ids_list[num_batches*batch_size:]
        batch_data = load_dataset_parallel(batch_ids, outputs_dir, targets_dir, networks_dir, features_dir, num_workers)
        name = os.path.join(save_path, f'{prefix}_{num_batches}.pt')
        torch.save(batch_data, name)
        del batch_data
        print(f"Saved remainder batch: {name}")
    #--------------
    elapsed = time.time() - start
    print(f">> Batch size: {batch_size} | Num batches: {num_batches} | Elapsed: {elapsed:.2f}s")



# We divide our data into 80/20 for cross validation.
import random
import numpy as np 

# Shuffle indixes
pass_ids = np.array(pass_ids) 
indexes = np.arange(0, len(pass_ids), 1)
np.random.shuffle(indexes) 

split = int(0.8 * len(indexes))
train_ids = pass_ids[indexes[:split]]
validation_ids = pass_ids[indexes[split:]]

#-----------------------
# Run parallelization
tensors_dir = '/mnt/data/sur/users/mrivera/Cuda-tensors/'
name = 'KBoost_v2_hybrid'
tensors_path = os.path.join(tensors_dir, name)
batching(ids_list = train_ids, 
        outputs_dir = outs_dir, 
        targets_dir = filtered_tgt_dir, 
        networks_dir = networks_dir, 
        features_dir =features_dir, 
        save_path = tensors_path, 
        batch_size=250, 
        num_workers=6,
        prefix='TrainBatch')

batching(ids_list = validation_ids, 
        outputs_dir = outs_dir, 
        targets_dir = filtered_tgt_dir, 
        networks_dir = networks_dir, 
        features_dir =features_dir, 
        save_path = tensors_path, 
        batch_size=100, 
        num_workers=6,
        prefix='ValBatch')

#===============================================
# Create Zip of batches
import shutil

# Zip entire directory
name = os.path.basename(tensors_path)
zip_dir = f'/mnt/data/sur/users/mrivera/Cuda-tensors/{name}'
shutil.make_archive(zip_dir, 'zip', tensors_path)

