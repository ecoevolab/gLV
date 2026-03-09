
from pathlib import Path

#----------------------------------------------------------
# Section: Mount-cluster
import os
import subprocess

# Paths
remote = "/mnt/data/sur/users/mrivera"
mount_p = "/home/mriveraceron/fenix_mount"

# Mount only if not already mounted
if os.path.ismount(mount_p):
    print(">> Cluster is already mounted")
else:
    subprocess.run([
        'sshfs',
        '-o', 'ro',   # <-- read-only option
        f'mrivera@fenix.lavis.unam.mx:{remote}',
        mount_p
    ], capture_output=True, text=True)
    print(">> Cluster mounted")

# Unmount command (if needed)
# subprocess.run(["fusermount", "-u", mount_p], check=True)
#---------------------------------------
# Section: Set seeds

def set_seed(seed=42):
    # Set ALL seeds for full reproducibility
    torch.manual_seed(seed)                 # Seed CPU 
    torch.cuda.manual_seed(seed)            # Seed GPU
    np.random.seed(seed)                    # Seed numpy
    random.seed(seed)                       # Seed python random
    torch.backends.cudnn.deterministic = True   # Ensure deterministic behavior
    torch.backends.cudnn.benchmark = False      

#---------------------------------------
import pandas as pd
import numpy as np
import random
import torch 
from datetime import datetime

# Section: Generate-paths
# Target-path
exp = "c748247a-8dc2"
# A_dir = os.path.join(mount_p, f"Experiments/{exp}/A-mat")
# tgt_dir = os.path.join(mount_p, f"Experiments/{exp}/Replica2/GNN-targets")
# data_path = os.path.join(mount_p, f"Data/{exp}.tsv")

exp_dir = "/home/mriveraceron/fenix_mount/Train-sims/4379fd40-9f0a"
A_dir = os.path.join(exp_dir, "A-mat")
tgt_dir = os.path.join(exp_dir, "GNN-targets")
odes_path = os.path.join(exp_dir, "raw-ODEs")
data_path = os.path.join(exp_dir, "parameters-sims.tsv")

# Generate quint for train_id
seednum=42
set_seed(seednum)  # Ensure reproducibility

# Generate ID for training.
timeID = datetime.now().strftime("Y%YM%mD%d")

#  Load-data
data = pd.read_csv(data_path, sep="\t")             # Load data
data_ids = data['id']                               # Extract ids

#------------------------------------------
# SECTION: Load-function
from torch_geometric.data import Data
import pyarrow.feather as feather

def load_single_data(id, A_dir, tgt_path):
    # Load adjacency matrix 
    A_path = os.path.join(A_dir, f"A_{id}.feather")
    A = pd.read_feather(A_path).to_numpy(dtype=np.float32)
    # Vector of edge weights
    row_idx, col_idx = np.nonzero(A)
    edge_weights = A[row_idx, col_idx]
    # Convert to torch tensors efficiently
    edge_index = torch.from_numpy(np.vstack([row_idx, col_idx]).astype(np.int32))
    edge_weights = torch.from_numpy(edge_weights)
    # Load target features 
    tgt_path = os.path.join(tgt_dir, f"tgt_{id}.feather")
    tgt_table =  feather.read_table(tgt_path, columns=['K_s'])
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

#----------------------------------------------------------
# Section: Parallelizing-data-loading

from concurrent.futures import ThreadPoolExecutor

def generate_data_parallel(idx, A_dir, tgt_dir, num_workers=4):  # idx is a list
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        data_list = list(executor.map(
            load_single_data,           
            idx,                       # List of IDs to iterate over
            [A_dir]*len(idx),          # Repeat A_dir for each ID
            [tgt_dir]*len(idx)         # Repeat tgt_dir for each ID
        ))
    return data_list


#----------------------------------------------------------
# SECTION: Define-GNN

import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv

class simple_gnn_gcn(nn.Module):
    def __init__(self, num_node_features=1, hidden_channels=64,  num_predictions=1):
        super().__init__()
        self.conv1 = GraphConv(num_node_features, hidden_channels)
        self.conv2 = GraphConv(hidden_channels, num_predictions)
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_weights
        x = self.conv1(x, edge_index, edge_weight)
        x = F.relu(x)
        x = self.conv2(x, edge_index, edge_weight)
        x = torch.sigmoid(x)  # Outputs between 0-1
        return x  # [num_nodes]


#----------------------------------------------------------
from tqdm import tqdm
import time

def train_batches(tidx, batch_size = 1000, epochs=500):
    num_batches = (len(tidx) + batch_size - 1) // batch_size
    start = time.time()
    model.train()
    # Empty lists for predictions, targets, loss at each epoch
    x_train, y_train, loss_epochs  = [], [], []
    for i in tqdm(range(0, len(tidx), batch_size), total=num_batches, desc="Loading batches"):
        data = [load_single_data(data_ids[id], A_dir, tgt_dir) for id in tidx[i:i + batch_size]]             # First 80 after shuffling
        Dloaded = DataLoader(data, batch_size=round(len(data)/10), shuffle=True)
        for epoch in range(epochs):
            total_loss = 0
            for data in Dloaded:
                optimizer.zero_grad()
                data = data.to(device)
                out = model(data)
                loss = loss_fn(out, data.y)
                loss.backward()
                optimizer.step()
                total_loss += loss.item()   # Accumulate loss
                if epoch==(epochs-1):
                    x_train.append(out.cpu().detach().numpy()) 
                    y_train.append(data.y.cpu().detach().numpy())
            loss_epochs.append(total_loss)
            if epoch % 200 == 0:
                print(f"Epoch {epoch}: Loss = {total_loss:.4f}")
    elapsed_time = time.time() - start
    print(f"Elapsed time for batching and training with batch: {elapsed_time:.2f}")

# Testing purpose 
train_batches(train_indices[:1000], batch_size = 100, epochs=500)
len(train_indices)
tidx = train_indices[:1000]
batch_size, epochs = 100, 500 

#-----------------------------------------------------------
num_batches = (len(x) + batch_size - 1) // batch_size

for i in tqdm(range(0, len(x), batch_size), total=num_batches, desc="Loading batches"):
    batch_ids = [data_ids[idx] for idx in x[i:i + batch_size]]
   
    data = generate_data_parallel(batch_ids, A_dir, tgt_dir, num_workers=8)            # First 80 after shuffling
    par_time = time.time() - start
    print(f">> The not parallelized time is of: {not_par_time:.2f}, while the parallel time is of: {par_time:.2f}")


data = generate_data_parallel(data_ids[idx], A_dir, tgt_dir) for idx in train_indices[i:i + batch_size]


#----------------------------------------------------------
# SECTION: Divide-data
from torch_geometric.loader import DataLoader
import random

# Load all data samples (for demo, we use only first 100 samples)
indices = list(range(1, len(data_ids)))  # Indices 1-100
random.shuffle(indices)  # Uses Python's random module (already seeded)

# Now select first 80 for training, rest for validation
indx = round(len(indices) * .8)
train_indices = indices[:indx]            # First 80 shuffled indices
val_indices = indices[indx:]              # Last 20 shuffled indices

# test
batch_size = 1000
num_batches = (len(indices) + batch_size - 1) // batch_size
id = data_ids[1]
train_data = generate_data_parallel(train_indices, data_ids, A_dir, tgt_dir, num_workers=4)



# Estimate memory per sample
mem_per_sample = sys.getsizeof(sample_data) / 100 / (1024**2)  # MB
total_mem_estimate = mem_per_sample * 10000  # MB for full dataset

print(f"Estimated memory for 10k samples: {total_mem_estimate:.2f} MB")


test = data_ids[train_indices[1]]
train_data = [load_single_data(data_ids[idx], A_dir, tgt_dir) for idx in train_indices]             # First 80 after shuffling
val_data = [load_single_data(data_ids[idx], A_dir, tgt_dir) for idx in val_indices]                 # Remaining 20 for validation

#----------------------------------------------------
# Section: Plotting
import matplotlib.pyplot as plt

def loss_plotter(loss_epochs = None, epochs = None, path = None ):
    # After collecting your data
    y = np.round(loss_epochs, 3)
    x = list(range(0, epochs))
    # Create scatter plot
    plt.figure(figsize=(8, 8))
    plt.plot(x, y, alpha=0.5)
    # Add perfect prediction line (y=x)
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title('Loss over epochs')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    ymin = min(y) 
    plt.ylim(ymin, max(y))
    plt.xlim(0, max(x))
    plt.savefig(path, dpi=150, bbox_inches='tight')

def preds_plotter(preds = None, tgts = None, path = None ):
    # After collecting your data
    preds = np.concatenate(preds)  # predictions
    tgts = np.concatenate(tgts)  # targets
    # Create scatter plot
    plt.figure(figsize=(8, 8))
    plt.scatter(preds, tgts, alpha=0.5)
   # Add perfect prediction line (y=x)
    plt.plot([0,  np.max(tgts)], [0,  np.max(tgts)], 'r--', label='Perfect prediction')
    plt.xlabel('Predictions')
    plt.ylabel('True Values')
    plt.title('Predictions vs True Values')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(0, max(tgts))
    plt.xlim(0, max(tgts))
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches='tight')
#--------------------------------------------------------
# Section: Example-run
import torch.optim as optim
import time

model = simple_gnn_gcn(num_node_features=1, hidden_channels=72, num_predictions=1).to('cuda' if torch.cuda.is_available() else 'cpu')
loss_fn = nn.MSELoss()                                                # Loss function for regression
optimizer = optim.Adam(model.parameters(), lr=0.01) 
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

start_time = time.time()
TR_load = DataLoader(train_data, batch_size=round(len(train_data)/10), shuffle=True)
model.train()
x_train, y_train, loss_epochs  = [], [], []

epochs=500
for epoch in range(epochs):
    total_loss = 0
    num_samples = 0
    for data in TR_load:
        optimizer.zero_grad()
        data = data.to(device)
        out = model(data)
        loss = loss_fn(out, data.y)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()   # Accumulate loss
        num_samples += data.num_graphs 
        if epoch==(epochs-1):
            x_train.append(out.cpu().detach().numpy()) 
            y_train.append(data.y.cpu().detach().numpy())
    # avg_loss = total_loss / num_samples
    loss_epochs.append(total_loss)
    if epoch % 50 == 0:
        print(f"Epoch {epoch}: Loss = {total_loss:.4f}")

elapsed_time = time.time() - start_time
print(f"Training completed in {elapsed_time:.2f} seconds")

# Version 1 is without sigmoid function
# Version 2 is with it
# Version 3 is with 1000 epochs
preds_plotter(preds = x_train, tgts = y_train, path = '/home/mriveraceron/glv-research/plots/GraphConv-pred_tgt-train.png')
loss_plotter(loss_epochs, epochs = 500, path = '/home/mriveraceron/glv-research/plots/GraphConv-LossEpochs-train.png')


#----------------------------------------------------------
# SECTION: Validation-loop
model.eval()  # Set to evaluation mode
total_loss = 0
num_samples = 0

start_time = time.time()
Val_load = DataLoader(val_data, batch_size=round(len(val_data)/10))
x_val, y_val = [], []

with torch.no_grad():  # Disable gradient computation
    for data in Val_load:
        data = data.to(device)
        out = model(data)
        loss = loss_fn(out, data.y)
        x_val.append(out.cpu().numpy())
        y_val.append(data.y.cpu().numpy())
        total_loss += loss.item() 
        num_samples += data.num_graphs
        
avg_loss = total_loss / num_samples
print(f"Validation Loss = {avg_loss:.4f}")
elapsed_time = time.time() - start_time
print(f"Validation completed in {elapsed_time:.2f} seconds")

preds_plotter(preds = x_train, tgts = y_train, path = '/home/mriveraceron/glv-research/plots/GraphConv-pred_tgt-val.png')
