
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
from torch_geometric.data import Data
from proquint import uint2quint, quint2uint

# Section: Generate-paths
# Target-path
exp = "c748247a-8dc2"
A_dir = os.path.join(mount_p, f"Experiments/{exp}/A-mat")
tgt_dir = os.path.join(mount_p, f"Experiments/{exp}/Replica2/GNN-targets")
data_path = os.path.join(mount_p, f"Data/{exp}.tsv")
odes_path = os.path.join(mount_p, f"Experiments/{exp}/raw-ODEs/raw-ODEs-{exp}.tsv")

# Generate quint for train_id
seednum=42
set_seed(seednum)  # Ensure reproducibility
num = np.random.randint(0, 2**32 - 1)               # Random seed for train_id      
train_id = uint2quint(num).split('-')[0]            # Generate quint

#  Load-data
data = pd.read_csv(data_path, sep="\t")             # Load data
data_ids = data['id']                               # Extract ids

 
# SECTION: Function
def load_single_data(exp_id, A_dir, tgt_path):
    # Get weights (x)
    A_path = os.path.join(A_dir, f"A_{exp_id}.feather")
    A = pd.read_feather(A_path).to_numpy(dtype=np.float32)
    A_tensor = torch.from_numpy(A.round(4))
    # Adjacency-matrix in COO format ([2, num_edges])
    edge_index = torch.nonzero(A_tensor, as_tuple=False).t().contiguous()
    edge_weights = edge_weights = A_tensor[A_tensor != 0]  # Shape: (|E|,) - only actual edges
    # Target-features                                                
    tgt_path = os.path.join(tgt_dir, f"tgt_{exp_id}.feather")
    tgt_numpy = pd.read_feather(tgt_path).filter(regex='^K_s').to_numpy(dtype=np.float32)
    tgt_tensor = torch.from_numpy(tgt_numpy)
    # Node-features
    n = A_tensor.shape[0]
    x_tensor = torch.ones(n,1)  
    # Create Data object
    data = Data(
        x= x_tensor,
        edge_weights = edge_weights,
        edge_index = edge_index,
        y=tgt_tensor
    )
    return data

#----------------------------------------------------------
# SECTION: Define-GNN

import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv

class simple_gnn_gcn(nn.Module):
    def __init__(self, num_node_features=1, hidden_channels=16,  num_predictions=1):
        super().__init__()
        self.conv1 = GraphConv(num_node_features, hidden_channels)
        self.conv2 = GraphConv(hidden_channels, num_predictions)
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_weights
        x = self.conv1(x, edge_index, edge_weight)
        x = F.relu(x)
        x = self.conv2(x, edge_index, edge_weight)
        return x  # [num_nodes]

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

train_data = [load_single_data(data_ids[idx], A_dir, tgt_dir) for idx in train_indices]             # First 80 after shuffling
val_data = [load_single_data(data_ids[idx], A_dir, tgt_dir) for idx in val_indices]                 # Remaining 20 for validation

#----------------------------------------------------
# Section: Plotting
import matplotlib.pyplot as plt

def plotter(x_axis = None, y_axis = None, path = None ):
    # After collecting your data
    x_axis = np.concatenate(x_axis)  # predictions
    y_axis = np.concatenate(y_axis)  # targets
    # Create scatter plot
    plt.figure(figsize=(8, 8))
    plt.scatter(x_axis, y_axis, alpha=0.5)
    # Add perfect prediction line (y=x)
    min_val = min(x_axis.min(), y_axis.min())
    max_val = max(x_axis.max(), y_axis.max())
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', label='Perfect prediction')
    plt.xlabel('Predictions')
    plt.ylabel('True Values')
    plt.title('Predictions vs True Values')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.axis('equal')  # Equal aspect ratio
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches='tight')

#--------------------------------------------------------
# Section: Example-run
import torch.optim as optim
import time

model = simple_gnn_gcn(num_node_features=1, hidden_channels=16, num_predictions=1).to('cuda' if torch.cuda.is_available() else 'cpu')
loss_fn = nn.MSELoss()                                                # Loss function for regression
optimizer = optim.Adam(model.parameters(), lr=0.01) 
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

start_time = time.time()
TR_load = DataLoader(train_data, batch_size=round(len(train_data)/10), shuffle=True)
model.train()
x_train = []
y_train = []
for epoch in range(100):
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
        if epoch==99:
            x_train.append(out.cpu().detach().numpy()) 
            y_train.append(data.y.cpu().detach().numpy())
    avg_loss = total_loss / num_samples
    print(f"Epoch {epoch}: Loss = {avg_loss:.4f}")

elapsed_time = time.time() - start_time
print(f"Training completed in {elapsed_time:.2f} seconds")

path =  '/home/mriveraceron/glv-research/plots/train-pred_tgt.png'
plotter(x_axis = x_train, y_axis = y_train, path = '/home/mriveraceron/glv-research/plots/GraphConv-V1-training.png')
#----------------------------------------------------------
# SECTION: Validation-loop
model.eval()  # Set to evaluation mode
total_loss = 0
num_samples = 0

start_time = time.time()
Val_load = DataLoader(val_data, batch_size=round(len(val_data)/10))
x_val = []
y_val = []
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

plotter(x_axis = x_val, y_axis = y_val, path = '/home/mriveraceron/glv-research/plots/GraphConv-V1-validation.png')
