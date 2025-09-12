import pandas as pd
import os

import torch 
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader  
from torch_geometric.nn import GCNConv
from sklearn.metrics import mean_squared_error, r2_score
import torch.nn.functional as F
import torch.optim as optim
import torch.nn as nn
import rpy2.robjects as ro
import numpy as np
from datetime import datetime
import subprocess


#=======================  Data loading =======================
# Section: Mount-cluster
remote = "/mnt/data/sur/users/mrivera"
mount_p = "/home/mriveraceron/fenix_mount"

# Mount-personal-dir
if os.path.ismount(mount_p):
    print(">> Cluster is already mounted") 
else:
    subprocess.run(['sshfs', f'mrivera@fenix.lavis.unam.mx:{remote}', mount_p],  capture_output=True, text=True)


# Section: Generate-paths
# Target-path
A_dir = os.path.join(mount_p, "Experiments/c748247a-8dc2/A-mat")
tgt_dir = os.path.join(mount_p, "Experiments/c748247a-8dc2/Replica2/GNN-targets")
data_path = os.path.join(mount_p, "Data/c748247a-8dc2.tsv")
odes_dir = os.path.join(mount_p, "Experiments/c748247a-8dc2/raw-ODEs")

data = pd.read_csv(data_path, sep="\t")
ids = data['id']
# ids_20_species = params[params['n_species'] == 20]['id']
# ids_100_species = params[params['n_species'] == 100]['id']

# SECTION: Function
# Functions-to-retrieve-data
def load_single_data(id, A_dir, odes_dir, tgt_dir):
    # Section: Edges weights
    A_path = os.path.join(A_dir, f"A_{id}.feather")
    A_mat = pd.read_feather(A_path).to_numpy(dtype=np.float32)
    A_tensor = torch.from_numpy(A_mat.round(4))
    edge_weights = A_tensor.flatten()  
    # Weights should be  a one-dimensional edge weight tensor
    # Section: Adjacency-matrix 
    edges = torch.nonzero(A_tensor, as_tuple=False).t().contiguous()
    # Section: Input-features
    feat_path = os.path.join(odes_dir, f"O_{id}.feather")
    feat_np = pd.read_feather(feat_path).iloc[:, 999].to_numpy(dtype=np.float32)
    tmp = torch.from_numpy(feat_np) 
    freq_tensor = tmp / tmp.sum()       # Relative-freq
    freq_tensor = freq_tensor.unsqueeze(1)
    # Section: Target-features                                                
    # Load Targets
    # Order:
    #   spec = i,                     # specie-extinct
    #   new_ext = new_ext,            # new-extinctions
    #   BC_diss = bray_curtis,        
    #   K_s = K_s,                    # Keystoness
    #   ext_ts = ext_ts               # Time to stability after perturbation
    tgt_path = os.path.join(tgt_dir, f"tgt_{id}.feather")
    tgt_numpy = pd.read_feather(tgt_path).filter(regex='^K_s').to_numpy(dtype=np.float32)
    tgt_tensor = torch.from_numpy(tgt_numpy)
    # Create Data object
    data = Data(
        x= freq_tensor,
        edge_weight = edge_weights,
        edge_index = edges,
        y=tgt_tensor
        # nsp=A_tensor.shape[0]
    )
    return data

# List to store the Data objects
data_list = []
id = ids[1]
for id in ids[:20]:
    data = load_single_data(id, A_dir, odes_dir, tgt_dir)
    data_list.append(data)
        
data = data_list[1]
#=======================  GCN =======================
#SECTION: Graph-Convolutional-Network
class GCNModel(torch.nn.Module):
    def __init__(self, in_features, hidden_features, num_predictions, num_layers, dropout_rate):
        super(GCNModel, self).__init__()
        self.num_layers = num_layers
        self.dropout_rate = dropout_rate
        self.convs = torch.nn.ModuleList()                                        # Create a list to store GCN layers
        self.convs.append(GCNConv(in_features, hidden_features))                  # First-layer
        for _ in range(num_layers-2):
            self.convs.append(GCNConv(hidden_features, hidden_features))          # Append-convolution
        self.convs.append(GCNConv(hidden_features, num_predictions))              # Output-layer 
        self.dropout = torch.nn.Dropout(dropout_rate)                             # Dropout
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_weight
        # Apply GCN layers with ReLU activation and dropout
        for i, conv in enumerate(self.convs):
            x = conv(x, edge_index, edge_weight)
            if i != len(self.convs) - 1:      # not the last layer
                x = F.relu(x)                   # Rectifier
                x = self.dropout(x)             # Dropout
        return x


#======================= Move model and data to the device =======================
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# in_features -> # of input features per node
# hidden_features ->  Hidden layer size (increased for better capacity)
# out_features -> # of prediction features done (prediction. It should be the same number as last layer size.
# num_predictions -> # predictions per node
# num_layers -> 
# droput_rate -> 
model = GCNModel(in_features=-1, hidden_features=40, num_predictions=1, num_layers=5, dropout_rate=0.2).to(device)
criterion = nn.MSELoss()                                                # Loss function for regression
optimizer = optim.Adam(model.parameters(), lr=0.01)                     # Optimizer

#======================================================
# SECTION: Divide-data
n_samples = len(data_list)  
n_val = int(0.2 * n_samples)
shuffled_indices = torch.randperm(n_samples)
train_indices = shuffled_indices[:-n_val]
val_indices = shuffled_indices[-n_val:]
 
train_data = [data_list[i] for i in train_indices.tolist()]
val_data = [data_list[i] for i in val_indices.tolist()]

#======================================================
# SECTION: Training
train_lines = []
def training_loop(n_epochs, train_data):
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
        train_lines.append('>> Epoch %d, Training Loss %f' % (epoch, float(total_train_loss)))

training_loop(10, train_data)

# SECTION: Create-log-TXT
# Create a unique log directory for this training run
log_dir = f"/home/mriveraceron/glv-research/train-logs/GCN-logs"
log_txt =  os.path.join(log_dir, f"{exp}.txt")

current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

with open(log_txt, 'w') as file:
    file.write("Starting training\n")
    file.write(f"{current_time}\n")
    file.write(f"Experiment ID: {exp}\n")
    file.write("-" * 40 + "\n")
    file.writelines(train_lines)


#======================================================
# SECTION: Validation
val_lines = []
def val_loop(n_epochs, val_data):
    for epoch in range(1, n_epochs + 1):
        model.eval()
        total_val_loss = 0
        with torch.no_grad():
            for data in val_data:
                data = data.to(device)
                out = model(data)
                loss = criterion(out, data.y)
                total_val_loss += loss.item()
        val_lines.append('>> Epoch %d, Validation Loss %f' % (epoch, float(total_val_loss)))

val_loop(100, val_data)


with open(log_txt, 'w') as file:
    file.write("-" * 40 + "\n")
    file.write("Starting validation\n")
    file.write(f"{current_time}\n")
    file.write("-" * 40 + "\n")
    file.writelines(val_lines)

#===========================================
for data in data_train: 
    print(data)
    # print("The weights are: ", data.edge_attr.dtype, "|| The inputs are: ", data.x.dtype, "|| The edges are: ", data.edge_index.dtype, " || The preds are: ", data.y.dtype )

#========================================== Set model to evaluation mode ==========================================
data_eval = DataLoader(data_list[int(len(data_list)*.8):],  shuffle=True)     # Use batch_size > 1 if all graphs are compatible
model.eval()  
all_preds = []
all_targets = []

with torch.no_grad():
    for data in data_eval:
        data = data.to(device)
        out = model(data)
        print("== this is the shape of Y-pred: ", out.shape)
        all_preds.append(out.cpu())
        print("== this is the shape of Y-true: ", data.y.shape, "\n")
        all_targets.append(data.y.cpu())


# Concatenate predictions and targets
print("== The shape of preditions are: ",  torch.cat(all_preds, dim=0).numpy().shape, "\n")
print("== The shape of Y is: ",  torch.cat(all_targets, dim=0).numpy().shape, "\n")

# Compute metrics
y_pred = torch.cat(all_preds, dim=0).numpy()
y_true = torch.cat(all_targets, dim=0).numpy()

def rmse_numpy(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred) ** 2))
rmse = rmse_numpy(y_true, y_pred)

print(f"\n Evaluation:")
print(f"RMSE: {rmse:.4f}")

# r2 = r2_score(targets, preds)  # Closer to 1 is better
# print(f"R² Score: {r2:.4f}, Accuracy: {r2*100:.2f}%")

#=======================  Verification =======================

# Check if PyTorch is installed and if CUDA (GPU support) is available
print("PyTorch version:", torch.__version__)            # Prints the PyTorch version
print("CUDA available:", torch.cuda.is_available())     # Checks if CUDA is available (for GPU usage)
print("Number of GPUs:", torch.cuda.device_count())     # Number of available GPUs

# Check device properties (including number of cores and more)
gpu_properties = torch.cuda.get_device_properties(0)
print("GPU Properties:", gpu_properties)


# ========================== Some noise =============
# print("# The number of species are:", Nspecs, " for simulation ", id, "\n")
R_seed = pd.read_csv(tsv_file, sep="\t").query("id == @id")["x0_seed"].iloc[0]            # Extract population seed of simulation
ro.r(f"set.seed({R_seed})")                                                               # Set seed (R)
ro.r(f"x0 <- stats::runif({Nspecs}, min = 0.1, max = 1)")                                 # Generate X0 (R)
x0 = torch.tensor(np.array(ro.r("x0")).round(4), dtype=torch.float32)   