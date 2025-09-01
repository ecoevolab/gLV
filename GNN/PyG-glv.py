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


# Libraries for summary
from torch.utils.tensorboard import SummaryWriter
from torch_geometric.loader import DataLoader
import os
from datetime import datetime

#=======================  Data loading =======================
dir = "/home/mriveraceron/data/2263e52c-8384"                    # Parent-directory
A_paths = os.path.join(dir, "A-mat")                             # Interactions-matrix
preds_paths = os.path.join(dir, "preds-ext")                     # Predictions

# FIX-ME
tsv_file = os.path.join(dir, "2263e52c-8384.tsv")               # Parameters-table
ids =  pd.read_csv(tsv_file, sep="\t")["id"].tolist()           # Extract simulation IDs     

# List to store the Data objects
data_list = []
id = ids[1]
for id in ids[:5]:
    #======================= Adjacency matrix (interactions) =======================
    A_path = os.path.join(A_paths, f"A_{id}.feather")
    A_mat = pd.read_feather(A_path).to_numpy(dtype=np.float32)
    A_tensor = torch.from_numpy(np.round(A_mat, 4))                                            # convert-to-Tensor
    x = A_tensor  
    edges = torch.nonzero(A_tensor, as_tuple=False).t().contiguous().to(torch.int32)                          # Adjacency FROM -> TO
    #======================= Predictions (metrics) =======================
    Nspecs = A_tensor.shape[0]                                                           # Input nodes
    # Direct numpy to tensor conversion
    pred_path = os.path.join(preds_paths, f"Preds_{id}.feather")
    pred_data = pd.read_feather(pred_path).drop(columns="id").to_numpy(dtype=np.float32)
    pred_tensor = torch.from_numpy(np.round(pred_data, 4))
    y = pred_tensor
    # Pre-allocate with correct dtype
    # y = torch.zeros(Nspecs, pred_tensor.shape[1], dtype=torch.float32)
    # y[:pred_tensor.shape[0], :] = pred_tensor                                                     
    #======================= Save Data =======================
    data = Data(x=x, edge_index=edges, y=y, nsp = Nspecs)
    # print("= The input X are: ", data.x, " the dimensions are: ", data.x.shape)
    # print("The weights are: ", data.edge_attr.dtype, "|| The inputs are: ", data.x.dtype, "|| The edges are: ", data.edge_index.dtype, " || The preds are: ", data.y.dtype )
    data_list.append(data)
        

#=======================  GCN =======================
class GCNModel(torch.nn.Module):
    def __init__(self, in_features, hidden_features, num_predictions):
        super(GCNModel, self).__init__()
        # First GCN layer
        self.conv1 = GCNConv(in_features, hidden_features)
        # Second GCN layer  
        self.conv2 = GCNConv(hidden_features, hidden_features)
        # Output layer for predictions
        self.output_layer = torch.nn.Linear(hidden_features, num_predictions)            # Regression 
        # Optional: Add dropout for regularization
        # self.dropout = torch.nn.Dropout(0.2)
    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        # First GCN layer with ReLU activation
        x = F.relu(self.conv1(x, edge_index))
        # Second GCN layer with ReLU activation
        x = F.relu(self.conv2(x, edge_index))
        # Output layer (no activation for regression)
        out = self.output_layer(x)
        return out     
    

#======================= Move model and data to the device =======================
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# in_features -> # of input features per node
# hidden_features ->  Hidden layer size (increased for better capacity)
# out_features -> # of prediction features done (prediction. It should be the same number as last layer size.
# num_predictions -> # predictions per node
model = GCNModel(in_features=20, hidden_features=60, num_predictions=5).to(device)
criterion = nn.MSELoss()                                                # Loss function for regression
optimizer = optim.Adam(model.parameters(), lr=0.01)                     # Optimizer

# FIXME
filtered_data_list = [data for data in data_list if data.nsp == 20]
model = GCNModel(in_features=20, hidden_features=60, num_predictions=5).to(device)


#======================= Train GCN =======================
data_train = DataLoader(data_list[:int(len(data_list)*.8)]) 
data_val = DataLoader(data_list[int(len(data_list)*.8):])     

# Create a unique log directory for this training run
log_dir = f"/home/mriveraceron/glv-research/train-logs/GCN-logs"
writer = SummaryWriter(log_dir)

# Track best validation loss for model saving
best_val_loss = float('inf')

for epoch in range(50):
    # ===== TRAINING PHASE =====
    model.train()                                        # training-mode
    total_train_loss = 0
    num_batches = 0
    for data in data_train:
        data = data.to(device)                          # Move-data-to-gpu
        optimizer.zero_grad()                           # Clear-gradients
        out = model(data)                               # output
        loss = criterion(out, data.y)
        loss.backward()                                 # Backpropagation
        optimizer.step()                                # Fix-weights 
        # measure-performance
        total_train_loss += loss.item()
        num_batches += 1
        avg_train_loss = total_train_loss / num_batches
    # ===== VALIDATION =====
    model.eval()
    total_val_loss = 0
    num_val_batches = 0
    with torch.no_grad():
        for data in data_val:
            data = data.to(device)
            out = model(data)
            val_loss = criterion(out, data.y)
            total_val_loss += val_loss.item()
            num_val_batches += 1
    # measure-loss
    avg_val_loss = total_val_loss / num_val_batches
    # ===== TENSORBOARD LOGGING =====
    # Log epoch-level metrics
    writer.add_scalar('Loss/Train', avg_train_loss, epoch)
    writer.add_scalar('Loss/Validation', avg_val_loss, epoch)
    # Log learning rate (if using scheduler)
    writer.add_scalar('Learning_Rate', optimizer.param_groups[0]['lr'], epoch)
    # Log model parameters (weights and gradients)
    for name, param in model.named_parameters():
        if param.grad is not None:
            writer.add_histogram(f'Parameters/{name}', param.data, epoch)
            writer.add_histogram(f'Gradients/{name}', param.grad.data, epoch)
    # Save best model
    if avg_val_loss < best_val_loss:
        best_val_loss = avg_val_loss
        torch.save(model.state_dict(), f'{log_dir}/best_model.pth')
        writer.add_text('Model_Save', f'Best model saved at epoch {epoch}', epoch)
    # Print progress
    if epoch % 5 == 0:
        print(f'Epoch {epoch:3d} | Train Loss: {avg_train_loss:.6f} | Val Loss: {avg_val_loss:.6f}')

# Close the writer
writer.close()

print(f"Training complete! Results saved at: {log_dir}")


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