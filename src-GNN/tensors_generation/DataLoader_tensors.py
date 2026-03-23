"""
Filtered nodes Data Generation Script for GNN model.
=================================================
Purpose:
    Generate tensors required for the GNN model.
    Species with low relative abundance at perturbation time are filtered out prior to tensor construction;
    feature dimensionality is preserved across the remaining nodes.

Important:
    Node features include the network statistics calculated using the whole network.
    Targets are the extinctions of the survival nodes performed in the communities.

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
import torch 
import pandas as pd
import numpy as np
import os
from datetime import datetime
import pyarrow.feather as feather
import pyarrow.compute as pc

# Section: Generate-paths
experiment_dir = '/mnt/data/sur/users/mrivera/clean_controls/91074c4e25b4'
networks_dir = os.path.join(experiment_dir, "Interactions")
targets_dir = os.path.join(experiment_dir, "Full_ExtSummaries")
features_dir = os.path.join(experiment_dir, "Topologies")
outs_dir = os.path.join(experiment_dir, "RawOutputs")

# Generate ID for training.
timeID = datetime.now().strftime("Y%YM%mD%d")

#  Load-data
data_path = os.path.join(experiment_dir, "simulation_summary.feather")
ids_list = feather.read_table(data_path).filter(pc.field('ext_performed') == True)['id'].to_numpy()

"""
We did not save the indices of which species should be filtered, so we will calculate them 
from the last column of the output. The last column corresponds to t=1000, the time at 
which the community is expected to be in quasi-stable state. NumPy's np.where() already 
returns 0-based Python indices, so no conversion is needed.

Node features and node statistics are calculated using the whole network.
"""
# SECTION: Generate function to load data for one simulation.
from torch_geometric.data import Data
import pyarrow.feather as feather

def load_single_data(id, output_dir, target_dir, networks_dir, features_dir):
    #-------------------------------
    # Section: Load which species to filter
    output_path = os.path.join(output_dir, f"RawOutput_{id}.feather")
    output = feather.read_table(output_path).to_pandas().iloc[:, 20].to_numpy()
    relative = output / output.sum()
    to_filter = np.where(relative > 1e-06)[0]
    #-------------------------------
    # Section: Read target features
    tgt_path = os.path.join(target_dir, f"ExtSummary_{id}.feather")
    tgt_table = feather.read_table(tgt_path, columns=['keystoneness']).to_pandas()
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
    edge_index_tensor = torch.from_numpy(np.stack([row_idx, col_idx]).astype(np.int64))
    edge_weights_tensor = torch.from_numpy(edge_weights).float()   
    #------------------------------ 
    # Section: Load node features 
    features_path = os.path.join(features_dir, f"Topology_{id}.feather")
    network_stats = feather.read_table(features_path).to_pandas().iloc[to_filter]
    x_tensor = torch.from_numpy(network_stats.to_numpy(dtype=np.float32))  
    # Add dummy data, vector of ones
    n = x_tensor.shape[0]
    x_tensor = torch.cat([x_tensor, torch.ones(n, 1, dtype=torch.float32)], dim=1)
    #------------------------------
    # Clean up large intermediate
    del output, relative, row_idx, col_idx, edge_weights, network_stats
    #------------------------------
    # Create Data object
    data = Data(
        x=x_tensor,
        edge_weights=edge_weights_tensor,
        edge_index=edge_index_tensor,
        y=y_tensor,
        num_nodes=n
    )
    return data

# Sanity check: load one data sample
id = ids_list[0]
load_single_data(id = id, output_dir = outs_dir, target_dir = targets_dir, networks_dir = networks_dir, features_dir = features_dir)

#-------------------------------------------------------------
# Section: Batch data
#-------------------------------------------------------------
import time, gc
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from tqdm import tqdm
import traceback
import torch.multiprocessing
torch.multiprocessing.set_sharing_strategy('file_system')

def batching(ids_list, outputs_dir, targets_dir, networks_dir, features_dir, save_path, 
             batch_size=250, num_workers=6, prefix='TrainBatch'):
    # Create directory to save tensors
    os.makedirs(save_path, exist_ok=True)
    start = time.time()
    func = partial(load_single_data, output_dir=outputs_dir, target_dir=targets_dir,networks_dir=networks_dir, features_dir=features_dir)
    #---------------------------------
    # Generate batches
    for i, batch_start in enumerate(range(0, len(ids_list), batch_size)):
        batch_ids = ids_list[batch_start : batch_start + batch_size]
        results = [None] * len(batch_ids)
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            future_to_idx = {executor.submit(func, id): j for j, id in enumerate(batch_ids)}
            for future in tqdm(as_completed(future_to_idx), total=len(batch_ids), desc=f"Batch {i}"):
                j = future_to_idx[future]
                try:
                    results[j] = future.result()
                except Exception as e:
                    print(f"  ID {batch_ids[j]} failed: {type(e).__name__}: {e}")
                    traceback.print_exc()
        # Saving path
        name = os.path.join(save_path, f'{prefix}_{i}.pt')
        torch.save(results, name)
        del results
        gc.collect()
        print(f">> Saved {name} \n")
    #------------------------------------
    elapsed = time.time() - start
    num_batches = -(-len(ids_list) // batch_size)  # ceiling division
    print(f">> Batch size: {batch_size} | Num batches: {num_batches} | Elapsed: {elapsed:.2f}s")

"""
We will divide our data in 80/20 for cross validation. 
Generate batches of data, save them and compress it to reduce storage.
"""
# We divide our data into 80/20 for cross validation.
import numpy as np 

# Shuffle indixes
indexes = np.arange(0, len(ids_list), 1)
np.random.shuffle(indexes) 

split = int(0.8 * len(indexes))
train_ids = ids_list[indexes[:split]]
validation_ids = ids_list[indexes[split:]]
#----------------------------------------------
# Run parallelization
#----------------------------------------------
tensors_dir = '/mnt/data/sur/users/mrivera/Cuda-tensors'
name = os.path.basename(experiment_dir)
tensors_path = os.path.join(tensors_dir, name)

batching_fn = partial(batching,
    outputs_dir   = outs_dir,
    targets_dir   = targets_dir,
    networks_dir  = networks_dir,
    features_dir  = features_dir,
    save_path     = tensors_path,
    num_workers   = 6
)

if __name__ == '__main__':
    batching_fn(ids_list=train_ids, batch_size=250, prefix='TrainBatch')
    batching_fn(ids_list=validation_ids, batch_size=250, prefix='ValBatch')

#===============================================
# Create Zip of batches
import shutil

# Zip entire directory
name = os.path.basename(tensors_path)
zip_dir = f'/mnt/data/sur/users/mrivera/Cuda-tensors/{name}'
shutil.make_archive(zip_dir, 'zip', tensors_path)

