import pandas as pd

import os
import subprocess

import torch 
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader  
from torch_geometric.nn import GCNConv
from torch_geometric.nn import GraphConv

from sklearn.metrics import mean_squared_error, r2_score
import torch.nn.functional as F
import torch.optim as optim
import torch.nn as nn
import rpy2.robjects as ro

import numpy as np
from datetime import datetime

import time
from pathlib import Path
from proquint import uint2quint, quint2uint

#=======================  Data loading =======================
# Section: Mount-cluster
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
# Section: Generate-paths
# Target-path
exp = "c748247a-8dc2"
A_dir = os.path.join(mount_p, f"Experiments/{exp}/A-mat")
tgt_dir = os.path.join(mount_p, f"Experiments/{exp}/Replica2/GNN-targets")
data_path = os.path.join(mount_p, f"Data/{exp}.tsv")
odes_path = os.path.join(mount_p, f"Experiments/{exp}/raw-ODEs/raw-ODEs-{exp}.tsv")
num = np.random.randint(0, 2**32 - 1)
train_id = uint2quint(num).split('-')[0]
data = pd.read_csv(data_path, sep="\t")
data_ids = data['id']
# ids_20_species = params[params['n_species'] == 20]['id']
# ids_100_species = params[params['n_species'] == 100]['id']

# SECTION: Function
def load_single_data(exp_id, A_dir, tgt_path):
    # Get weights (x)
    A_path = os.path.join(A_dir, f"A_{exp_id}.feather")
    A = pd.read_feather(A_path).to_numpy(dtype=np.float32)
    A_tensor = torch.from_numpy(A.round(4))
    # Adjacency-matrix in COO format ([2, num_edges])
    edge_index = torch.nonzero(A_tensor, as_tuple=False).t().contiguous()
    edge_weights = edge_weights = A_tensor[A_tensor != 0]  # Shape: (|E|,) - only actual edges
    # Section: Target-features                                                
    # Load Targets
    # Order:
    #   new_ext = new_ext,            # new-extinctions
    #   BC_diss = bray_curtis,        
    #   K_s = K_s,                    # Keystoness
    #   ext_ts = ext_ts               # Time to stability after perturbation
    tgt_path = os.path.join(tgt_dir, f"tgt_{exp_id}.feather")
    tgt_numpy = pd.read_feather(tgt_path).filter(regex='^K_s').to_numpy(dtype=np.float32)
    tgt_tensor = torch.from_numpy(tgt_numpy)
    # Section: Node-features
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
        return x.view(-1)  # [num_nodes]
    
#=======================  Training  =======================
dlist = [load_single_data(eid, A_dir, tgt_dir) for eid in data_ids[1:101]]  # Load first 100 samples for demo

def training_loop(n_epochs, train_data, model, criterion, optimizer, device, exp, train_id, save_dir='/home/mriveraceron/glv-research/GNN-params'):
    train_lines = []
    model.to(device)
    # Initialize tracking variables
    best_epoch = 1
    best_loss = float('inf')
    for epoch in range(1, n_epochs + 1):
        start_time = time.time()
        # Training phase
        model.train()
        total_train_loss = 0.0
        batch_count = 0
        for data in train_data:
            data = data.to(device, non_blocking=True)  # Async GPU transfer
            # Forward pass
            optimizer.zero_grad()
            out = model(data)
            loss = criterion(out, data.y)
            # Backward pass
            loss.backward()
            optimizer.step()
            total_train_loss += loss.item()
            batch_count += 1
        # Calculate average loss per batch
        avg_train_loss = total_train_loss / batch_count if batch_count > 0 else 0.0
        # Track best epoch
        if epoch ==1:
            best_loss = avg_train_loss
            best_epoch = epoch
        if avg_train_loss < best_loss:
            best_epoch = epoch
            best_loss = avg_train_loss
        # Calculate epoch duration
        epoch_duration = time.time() - start_time
        # Logging
        epoch_line = f">> Epoch {epoch:3d}, Average total training Loss: {avg_train_loss:.6f}"
        duration_line = f">> Epoch duration: {epoch_duration:.2f} seconds"   
        print(epoch_line)
        print(duration_line)
        train_lines.extend([epoch_line + '\n', duration_line + '\n'])
    # Final best epoch summary
    best_line = f">> Best epoch: {best_epoch} with loss: {best_loss:.6f}"
    print(best_line)
    train_lines.append(best_line + '\n')
    # Save model checkpoint
    save_path = Path(save_dir) / f'{train_id}-Exp_{exp}.pt'
    save_path.parent.mkdir(parents=True, exist_ok=True)  # Create directory if needed
    torch.save({
        'epoch': n_epochs,
        'best_epoch': best_epoch,
        'best_loss': best_loss,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'final_loss': avg_train_loss
    }, save_path)
    return train_lines, save_path

nlayer = 3
neurons = 16
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = simple_gnn_gcn(num_node_features=1, hidden_channels=neurons, num_predictions=1)
criterion = nn.MSELoss()                                                # Loss function for regression
optimizer = optim.Adam(model.parameters(), lr=0.01)                     # Optimizer
train_data = dlist[1:80]             # Data for training

now = time.time()
ltrain, save_path = training_loop(100, train_data, model, criterion, optimizer, device, exp, train_id)
train_duration = time.time() - now 

#======================================================
# SECTION: Validation
val_data = dlist[81:]            
def val_loop(val_data, save_path):
    val_lines = []
    model.load_state_dict(torch.load(save_path,weights_only=True)['model_state_dict'])
    model.eval()
    total_val_loss = 0
    with torch.no_grad():
        for data in val_data:
            data = data.to(device)
            out = model(data)
            loss = criterion(out, data.y)
            total_val_loss += loss.item()
    avg_val_loss = total_val_loss / len(val_data)
    val_lines.append(f'>> Validation Loss is {avg_val_loss:.6f}\n')
    return val_lines

# Get validation lines
now = time.time()
lval = val_loop(val_data, save_path )
val_duration = time.time() - now 

# SECTION: Create-log-TXT
# Create a unique log directory for this training run
log_dir = f"/home/mriveraceron/glv-research/GNN-Logs"
log_txt =  os.path.join(log_dir, f"{train_id}-Exp_{exp}.txt")
print(f">> Log file will be saved to: {log_txt}")

if os.path.exists(log_txt):  # check if file exists
    os.remove(log_txt)       # delete it

now = datetime.now()
formatted_time = now.strftime("%Y-%m-%d %H:%M:%S.%f")
with open(log_txt, 'w') as file:
    file.write("-" * 40 + "\n")
    file.write("Starting training\n")
    file.write(f"Timestamp: {formatted_time}\n")
    file.write(f"Experiment ID: {exp}\n")
    file.write(f"The number of layer used is: {nlayer}\n")
    file.write(f"The number of neurons used is: {neurons}\n")
    file.write(f"The time elapsed for training was: {train_duration:.2f} seconds\n")
    file.write("-" * 40 + "\n")
    file.writelines(ltrain)
    file.write("-" * 40 + "\n")
    file.write("Starting validation\n")
    file.write(f"The time elapsed for validation was: {val_duration}\n")
    file.writelines(lval)
    file.write("-" * 40 + "\n")
   
    

