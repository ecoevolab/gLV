

from pathlib import Path

#=======================  Data loading =======================
import os
import subprocess

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
#=======================  Set seeds  =======================
# Section: Set seeds
def set_seed(seed=42):
    # Set ALL seeds for full reproducibility
    torch.manual_seed(seed)                 # Seed CPU 
    torch.cuda.manual_seed(seed)            # Seed GPU
    np.random.seed(seed)                    # Seed numpy
    random.seed(seed)                       # Seed python random
    torch.backends.cudnn.deterministic = True   # Ensure deterministic behavior
    torch.backends.cudnn.benchmark = False      

#==============================================
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
    ## Order:
    ##   new_ext = new_ext,            # new-extinctions
    ##  BC_diss = bray_curtis,        
    ##   K_s = K_s,                    # Keystoness
    ##   ext_ts = ext_ts               # Time to stability after perturbation
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
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv

# SECTION: Define-GNN
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

#----------------------------------------------------------
# SECTION: Divide-data
from torch_geometric.loader import DataLoader
import random

# Load all data samples (for demo, we use only first 100 samples)
indices = list(range(1, 100))  # Indices 1-100
random.shuffle(indices)  # Uses Python's random module (already seeded)

# Now select first 80 for training, rest for validation
indx = round(len(indices) * .8)
train_indices = indices[:indx]            # First 80 shuffled indices
val_indices = indices[indx:]              # Last 20 shuffled indices

train_data = [load_single_data(data_ids[idx], A_dir, tgt_dir) for idx in train_indices]             # First 80 after shuffling
val_data = [load_single_data(data_ids[idx], A_dir, tgt_dir) for idx in val_indices]                 # Remaining 20 for validation

#----------------------------------------------------------
# SECTION: Training-loop
import torch.optim as optim

src_file = "/home/mriveraceron/glv-research/gLV/GNN/src-TrVal/training_src.py"
exec(open(src_file).read())

nlayer = 3
neurons = 16
model = simple_gnn_gcn(num_node_features=1, hidden_channels=neurons, num_predictions=1)
criterion = nn.MSELoss()                                                # Loss function for regression
optimizer = optim.Adam(model.parameters(), lr=0.01) 
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
batch_size=16

final_save_path, best_model_path = optimized_training_loop(100, train_data, model, criterion, optimizer, device, exp, train_id, 
                          params_dir='/home/mriveraceron/glv-research/GNN-params', 
                          logs_dir = '/home/mriveraceron/glv-research/GNN-Logs',
                          seed=seednum, batch_size=batch_size, use_mixed_precision=True, patience=20, 
                          save_every=50, gradient_clip_val=1.0, nlayer= nlayer, neurons= neurons)


















# FIX ME
#======================================================
# SECTION: Validation-loop          
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
    file.write("-" * 40 + "\n")
    file.writelines(ltrain)
    file.write("-" * 40 + "\n")
    file.write("Starting validation\n")
    file.write(f"The time elapsed for validation was: {val_duration}\n")
    file.writelines(lval)
    file.write("-" * 40 + "\n")
   
    

