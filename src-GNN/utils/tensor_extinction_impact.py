"""
Generate tensors for sub and full community impact
=================================================
Description:
    -> Generates graph tensors for GNN training and evaluation.
    -> Only surviving species are retained — defined as those with a relative 
        abundance >= 1e-06 at perturbation time. Node features, targets, and 
        the adjacency matrix are all filtered consistently against this mask.
    -> Two tensor sets are produced (sub, full), sharing the same node features,
        adjacency matrix and edge weights. They differ only in how the target 
        tensor (keystoneness) is computed:
            - sub:  computed over surviving species only
            - full: computed over the full community
Data Format:
    - X : Node feature matrix. Each row corresponds to its respective node statistics in the network.
    - Y : Target matrix. Each row contains the metric(s) of interest for the corresponding node.
    - edge_index : Graph connectivity in who to whom format, encoding which species are connected to which.
    - edge_attr  : Edge weights representing the strength of those interaction.

Author: Manuel Rivera
Date:   3-May-2026
"""

# Imports
import numpy as np
import os
import pyarrow.feather as feather

import torch
from torch_geometric.data import Data
import pyarrow.feather as feather

import time, gc
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from functools import partial
from tqdm import tqdm
import traceback
import torch.multiprocessing
torch.multiprocessing.set_sharing_strategy('file_system')

import zipfile
import io

# Section: Generate-paths
experiment_dir = '/mnt/data/sur/users/mrivera/Data/clean_controls/Kboost_proof_0c7b0f'
A_dir = os.path.join(experiment_dir, "Interactions")
sub_dir = os.path.join(experiment_dir, "Sub_ExtSummaries")
full_dir = os.path.join(experiment_dir, "Full_ExtSummaries")
x_dir = os.path.join(experiment_dir, "Topologies")
outs_dir = os.path.join(experiment_dir, "RawOutputs")

#  Load-data
data_path = os.path.join(experiment_dir, "simulation_summary.feather")
ids_list = feather.read_table(data_path)['id'].to_numpy()

# SECTION: Generate function to load data for one simulation
def generate_data(id, output_dir, y_dir, A_dir, x_dir, sub=False):
    # Section: Load filter mask
    output_col = feather.read_table( os.path.join(output_dir, f"RawOutput_{id}.feather") ).column(20).to_numpy()
    rel = output_col/output_col.sum()
    to_filter = np.where(rel > 1e-6)[0]
    # Secion: Load target features
    if sub:
        # Important note: targets are already filtered with sub-community impact
        y = feather.read_table(f"{y_dir}/ExtSummary_{id}.feather", columns=['keystoneness'])
        y_tensor = torch.from_numpy(y.column('keystoneness').to_numpy().astype(np.float32)).unsqueeze(1)
    else: 
        # Important note: targets are unfiltered with full impact
        y = feather.read_table(f"{y_dir}/ExtSummary_{id}.feather", columns=['keystoneness'])
        y_tensor = torch.from_numpy(y.column('keystoneness').to_numpy()[to_filter].astype(np.float32) ).unsqueeze(1)
    # Section: Load adjacency matrix
    A_filter = feather.read_table(f"{A_dir}/A_{id}.feather").to_pandas().to_numpy()[np.ix_(to_filter, to_filter)]
    row_idx, col_idx = np.nonzero(A_filter)
    edge_index_tensor = torch.from_numpy( np.stack([row_idx, col_idx]).astype(np.int64) )
    edge_weights_tensor = torch.from_numpy( A_filter[row_idx, col_idx].astype(np.float32) )
    # Section: Load node features
    x_np = feather.read_table(f"{x_dir}/Topology_{id}.feather").to_pandas().iloc[to_filter].to_numpy(dtype=np.float32)
    dummy = torch.ones(x_np.shape[0], 1, dtype=torch.float32)
    x_tensor = torch.from_numpy( np.hstack([x_np, dummy]) )
    return Data(
        x=x_tensor,
        edge_index=edge_index_tensor,
        edge_weights=edge_weights_tensor,
        y=y_tensor,
    )

# Sanity check: load one data sample
id = ids_list[0]
output_dir = outs_dir

data_full = generate_data(id = id, output_dir = outs_dir, y_dir = full_dir, A_dir=A_dir, x_dir =x_dir, sub=False)
for key, value in data_full:
    print(f"{key}: {value.dtype}")

data_sub = generate_data(id = id, output_dir = outs_dir, y_dir = sub_dir, A_dir=A_dir, x_dir =x_dir, sub=True)
for key, value in data_sub:
    print(f"{key}: {value.dtype}")


#-------------------------------------------------------------
# Section: Batch data
#-------------------------------------------------------------
"""
    Data will be stored in a zip file, with each batch as a separate .pt file inside. 
    Multiprocessing is useed to speed up the loading and processing of individual samples, and tqdm for progress tracking.
"""
def batching(ids_list, partial_fun, save_path, batch_size=250, num_workers=6, prefix='TrainBatch', overwrite=False):
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    zip_mode = 'w' if (overwrite or not os.path.exists(save_path)) else 'a'
    
    def serialize_and_write(zf, batch_name, results):
        """Runs in a separate thread so workers don't idle during zip writes."""
        buffer = io.BytesIO()
        torch.save(results, buffer)
        zf.writestr(batch_name, buffer.getvalue())
        print(f">> Written {batch_name} into zip")
    start = time.time()
    batches = [ids_list[i:i + batch_size] for i in range(0, len(ids_list), batch_size)]
    with zipfile.ZipFile(save_path, zip_mode, compression=zipfile.ZIP_DEFLATED) as zf, \
         ProcessPoolExecutor(max_workers=num_workers) as executor, \
         ThreadPoolExecutor(max_workers=1) as writer:  # 1 writer thread
        write_future = None
        for i, batch_ids in enumerate(batches):
            # Load batch (ordered, simpler than submit/as_completed)
            results = list(tqdm(
                executor.map(partial_fun, batch_ids, chunksize=max(1, len(batch_ids) // num_workers)),
                total=len(batch_ids),
                desc=f"Batch {i+1}/{len(batches)}"
            ))
            # Wait for previous write to finish before submitting next
            if write_future:
                write_future.result()
            # Submit write in background while next batch loads
            write_future = writer.submit(serialize_and_write, zf, f'{prefix}_{i}.pt', results)
            del results
        # Ensure last write completes
        if write_future:
            write_future.result()
    elapsed = time.time() - start
    print(f">> Batch size: {batch_size} | Num batches: {len(batches)} | Elapsed: {elapsed:.2f}s")
    print(f'>> Completed batching. Saved at: {save_path}')

# Generate sub-community impact tensors
func = partial(generate_data, output_dir = outs_dir, y_dir = sub_dir, A_dir=A_dir, x_dir =x_dir, sub=True)
tensors_dir = '/mnt/data/sur/users/mrivera/Tensors/clean_ctrls'
tensors_path = f'{tensors_dir}/gnn_proof_sub.zip'
if os.path.exists(tensors_path):
    os.remove(tensors_path)

batching(
    ids_list      = ids_list,        # list of IDs to process
    partial_fun   = func,
    save_path     = tensors_path,
    batch_size    = 250,
    num_workers   = 6,
    prefix        = 'TrainBatch',
    overwrite     = False
)

# Generate full-community impact tensors
func = partial(generate_data, output_dir = outs_dir, y_dir = full_dir, A_dir=A_dir, x_dir =x_dir, sub=False)
tensors_dir = '/mnt/data/sur/users/mrivera/Tensors/clean_ctrls'
tensors_path = f'{tensors_dir}/gnn_proof_full.zip'
if os.path.exists(tensors_path):
    os.remove(tensors_path)
    
batching(
    ids_list      = ids_list,        # list of IDs to process
    partial_fun   = func,
    save_path     = tensors_path,
    batch_size    = 250,
    num_workers   = 6,
    prefix        = 'TrainBatch',
    overwrite     = False
)