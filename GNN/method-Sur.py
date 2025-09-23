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


data = pd.read_csv(data_path, sep="\t")
data_ids = data['id']
sim_id = data_ids[1]
# ids_20_species = params[params['n_species'] == 20]['id']
# ids_100_species = params[params['n_species'] == 100]['id']

# SECTION: Function
# Functions-to-retrieve-data
def load_single_data(exp_id, A_dir, tgt_path):
    # Get weights (x)
    A_path = os.path.join(A_dir, f"A_{exp_id}.feather")
    A = pd.read_feather(A_path).to_numpy(dtype=np.float32)
    A_tensor = torch.from_numpy(A.round(4))
    # Adjacency-matrix in COO format ([2, num_edges])
    edges = torch.nonzero(A_tensor, as_tuple=False).t().contiguous()
    # Section: Target-features                                                
    # Load Targets
    # Order:
    #   spec = i,                     # specie-extinct
    #   new_ext = new_ext,            # new-extinctions
    #   BC_diss = bray_curtis,        
    #   K_s = K_s,                    # Keystoness
    #   ext_ts = ext_ts               # Time to stability after perturbation
    tgt_path = os.path.join(tgt_dir, f"tgt_{exp_id}.feather")
    tgt_numpy = pd.read_feather(tgt_path).filter(regex='^K_s').to_numpy(dtype=np.float32)
    tgt_tensor = torch.from_numpy(tgt_numpy)
    # Create Data object
    data = Data(
        x= A_tensor,
        edge_index = edges,
        y=tgt_tensor
        # nsp=A_tensor.shape[0]
    )
    return data


class simple_gnn_gcn(nn.Module):
    def __init__(self, num_node_features=20, hidden_channels=16,  num_predictions=1):
        super().__init__()
        self.conv1 = GraphConv(num_node_features, hidden_channels)
        self.conv2 = GraphConv(hidden_channels, num_predictions)
    def forward(self, data):
        x, edge_index = data.x, data.edge_index 
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index)
        return x  # [num_nodes]
    
#=======================  Training  =======================
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = simple_gnn_gcn(num_node_features=20, hidden_channels=16, num_predictions=1)
criterion = nn.MSELoss()                                                # Loss function for regression
optimizer = optim.Adam(model.parameters(), lr=0.01)                     # Optimizer

exp_id = data_ids[1]
total_train_loss = 0
d = load_single_data(exp_id, A_dir, tgt_dir).to(device)
dlist = [load_single_data(eid, A_dir, tgt_dir) for eid in data_ids[1:101]]  # Load first 100 samples for demo



def training_loop(n_epochs, train_data, model, criterion, optimizer, device,exp):
    train_lines = []
    model.to(device)
    for epoch in range(1, n_epochs + 1):
        # ===== TRAINING PHASE =====
        model.train()                                        # training-mode
        total_train_loss = 0
        for data in train_data:
            data = data.to(device)                          # Copy-data-to-gpu
            optimizer.zero_grad()                           # Clear-gradients
            out = model(data)                               # output
            loss = criterion(out, data.y)                   # MSE-loss
            loss.backward()                                 # Backpropagation
            optimizer.step()                                # Update-weights 
            total_train_loss += loss.item()
        line = '>> Epoch %d, Training Loss %f \n' % (epoch, float(total_train_loss))
        print(line)
        train_lines.append(line)
        # Look for best epoch
        if epoch == 1:
            best_epoch = epoch
            best_loss = float(total_train_loss)
        else:
            if float(total_train_loss) < best_loss:
                best_epoch = epoch
                best_loss = float(total_train_loss)
    line = '>> The best epoch was %d with loss %f \n' % (best_epoch, best_loss)
    print(line)
    train_lines.append(line)
    torch.save({
    'epoch': n_epochs,
    'model_state_dict': model.state_dict(),
    'optimizer_state_dict': optimizer.state_dict(),
    'loss': loss
    }, f'/home/mriveraceron/glv-research/train-logs/GCN-params/WeightsEp{n_epochs}-Exp-{exp}.pt')
    return train_lines

train_data = dlist[1:80]             # Data for training
ltrain = training_loop(100, train_data, model, criterion, optimizer, device, exp)


# SECTION: Create-log-TXT
# Create a unique log directory for this training run
log_dir = f"/home/mriveraceron/glv-research/train-logs/GCN-logs"
log_txt =  os.path.join(log_dir, f"{exp}.txt")

if os.path.exists(log_txt):  # check if file exists
    os.remove(log_txt)       # delete it

current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

with open(log_txt, 'w') as file:
    file.write("Starting training\n")
    file.write(f"{current_time}\n")
    file.write(f"Experiment ID: {exp}\n")
    file.write("-" * 40 + "\n")
    file.writelines(ltrain)


#======================================================
# SECTION: Validation
val_data = dlist[81:]            
def val_loop(n_epochs, val_data):
    val_lines = []
    model.load_state_dict(torch.load(f'/home/mriveraceron/glv-research/train-logs/GCN-params/WeightsEp{n_epochs}-Exp-{exp}.pt',
                                         weights_only=True)['model_state_dict'])
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
lval = val_loop(100, val_data)

# Open file in append mode and add validation results
with open(log_txt, 'a') as file:
    file.write("-" * 40 + "\n")
    file.write("Starting validation\n")
    file.write(f"{current_time}\n")
    file.write("-" * 40 + "\n")
    file.writelines(lval)
#======================================================
# SECTION: Load-checkpoint
# checkpoint = torch.load(f'checkpoint_epoch_{n_epochs}.pt')
# model.load_state_dict(checkpoint['model_state_dict'])
# optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
# epoch = checkpoint['epoch']
# loss = checkpoint['loss']
# model.eval()  # set model to evaluation mode
