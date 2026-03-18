"""
GNN Model Training
=================================================
Purpose:
    Train a Graph Neural Network (GNN) using pre-built tensors. 
    Species with low relative abundance at perturbation time are pre-filtered prior to training.

Input Data:
    - X          : Node feature matrix, where each row represents the statistics of a node in the network.
    - Y          : Target variable representing keystoneness.
    - edge_index : Graph connectivity in COO format, defining species interactions.
    - edge_attr  : Edge weights representing the strength of each interaction.

Dependencies:
    torch==2.8, pandas==2.3.3, numpy==2.0.2

Author: Manuel Rivera
Date:   March 13, 2026
"""
#-------------------------------
"""
Model declaration

5-layer architecture:
    - Layer 1:   13 input features  →  64 channels
    - Layers 2-4: 64 channels       →  64 channels
    - Layer 5:   64 channels        →   1 channel

Output activation: Sigmoid, to match the keystoneness target range.
"""
# Section: Generate model
import torch 
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv

class model(nn.Module):
    def __init__(self, hidden_channels=64, num_layers=5):
        super().__init__()
        self.convs = nn.ModuleList()
        # First layer: 1 -> hidden_channels
        self.convs.append(GraphConv(13, hidden_channels))
        # Middle layers: hidden_channels -> hidden_channels
        for _ in range(num_layers - 2):
            #self.convs.append(GATConv(hidden_channels*heads, hidden_channels, heads=heads))
            self.convs.append(GraphConv(hidden_channels, hidden_channels))
        # Last layer: hidden_channels -> 1
        self.convs.append(GraphConv(hidden_channels, 1))
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_weights
        # Apply all layers except the last
        for i, conv in enumerate(self.convs[:-1]):
            x = conv(x, edge_index, edge_weight)
            x = F.relu(x)
        # Apply last layer with sigmoid
        x = self.convs[-1](x, edge_index, edge_weight)
        x = torch.sigmoid(x)
        return x  # [num_nodes]
    

"""
Training function.

At the final epoch, the model returns:
    - The expected and predicted network node with the highest keystoneness value.
    - The expected and predicted keystoneness values.

Model weights are saved for reproducibility.
"""
import glob
import time
import random
import pandas as pd
from tqdm import tqdm
import os

def training_loop(model_declared, device, batched_paths, weights_path, loss_fn, optimizer, epochs=100):
    #------------------------------------------
    model_declared.train()
    loss_history  =  []             # Loss at epoch
    total_elapsed = 0               # Running time
    metrics_pred, metrics_true = [], []   # Declare list for metrics
    idx_max_true, idx_max_pred = [], []     # Declare list for max values
    #---------------------
    # Sectioin: Load all data once before training
    all_data = []
    for path in batched_paths:
        all_data.extend(torch.load(path, weights_only=False))
    all_data = [data.to(device) for data in all_data]  # move to GPU once
    for epoch in tqdm(range(epochs), desc="Training"):
        start = time.time()
        epoch_loss = 0
        for data in all_data:
            #----------------------
            # Move it to device and run model
            # data = data_list[0]
            optimizer.zero_grad()
            out = model_declared(data)
            loss = loss_fn(out, data.y)
            loss.backward()
            #torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            epoch_loss += loss.item()   # Accumulate loss
            #----------------------
            # Section: At last epochs
            #----------------------
            if epoch == epochs - 1:  
                #----------------------
                # Append values and node with maximum value
                idx_max_true.append(torch.argmax(data.y, dim=0).detach())
                idx_max_pred.append(torch.argmax(out, dim=0).detach())
                metrics_true.append(data.y[:, 0].detach())
                metrics_pred.append(out[:, 0].detach())
        #----------------------
        # Append epoch loss to history
        loss_history.append(epoch_loss)
        elapsed = time.time() - start
        total_elapsed += elapsed
        # Print every n epochs
        if epoch % 10 == 0:
            tqdm.write(f"Epoch {epoch}: Loss = {epoch_loss:.4f}, Elapsed time: {elapsed:.2f}")
    #----------------------
    # Save weights
    #----------------------
    torch.save({
        'epoch': epochs - 1,
        'model_state_dict': model_declared.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'loss': epoch_loss
    }, weights_path)
    # Summary
    print(f'>> the total elapsed time with {epochs} epochs is {total_elapsed:.2f} seconds ( {total_elapsed/60:.2f} minutes)')   
    return  loss_history, metrics_true, metrics_pred, idx_max_true, idx_max_pred


"""
Model preparation

- Seeds are fixed for reproducibility.
- Loss function: MSE.
- Device: CUDA.
- Model is moved to device.
- Optimizer: Adam.
- Number of training epochs is set.
- Learning rate is set to 0.001
"""
import torch.optim as optim
import numpy as np
import random

def seed_fn(seed=42):
    # Set ALL seeds for full reproducibility
    torch.manual_seed(seed)                 # Seed CPU 
    torch.cuda.manual_seed(seed)            # Seed GPU
    np.random.seed(seed)                    # Seed numpy
    random.seed(seed)                       # Seed python random
    torch.backends.cudnn.deterministic = True   # Ensure deterministic behavior
    torch.backends.cudnn.benchmark = False 


# Set seed and run model
loss_fn = nn.MSELoss()
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
n_seed = 42
seed_fn(n_seed)
model_declared = model(hidden_channels=64, num_layers=5).to(device)
optimizer = optim.Adam(model_declared.parameters(), lr=0.01)
epochs = 600
print('The number of epochs will be:', epochs, '\n')
print('The optimizer LR will be:', optimizer.param_groups[0]['lr'], '\n')

#-------------------------------
# Section: Generate paths
#-------------------------------
import os
import glob

# Experiment data
tensors_dir = '/home/mriveraceron/glv-research/data_tensors/'
experiment_name = 'KBoost_v2_filter/'
experiment_data = os.path.join(tensors_dir,experiment_name)
batched_paths = glob.glob(f"{experiment_data}/TrainBatch_*.pt")

# Results path
result_dir = '/home/mriveraceron/glv-research/Results/'
result_path = os.path.join(result_dir, experiment_name, 'LRdefault')
os.makedirs(result_path, exist_ok=True)
print('The results directory will be:', result_path, '\n')

# Declare directory to save model weights
# weights_dir = '/home/mriveraceron/glv-research/model_weights/'
weights_path = os.path.join(result_path,'model_weights.pth')

#-------------------------------
# Section: Run Model
#-------------------------------
loss_history, metrics_true, metrics_pred, idx_max_true, idx_max_pred = training_loop(model_declared, device, batched_paths, weights_path, loss_fn, optimizer, epochs)

# List of node indexes and values
idx_max_true = torch.stack(idx_max_true).cpu().numpy()
idx_max_pred = torch.stack(idx_max_pred).cpu().numpy()
metrics_true = torch.cat(metrics_true).cpu().numpy()
metrics_pred = torch.cat(metrics_pred).cpu().numpy()
loss_history = np.array(loss_history)

# Save
np.savez(f'{result_path}/model_results.npz',
    max_idx_true  = idx_max_true,
    max_idx_pred  = idx_max_pred,
    values_true   = metrics_true,
    values_pred   = metrics_pred,
    loss_history  = loss_history
)
# To load back:
# data = np.load(f'{result_path}/model_results.npz')
# idx_max_true  = data['max_idx_true']
# idx_max_pred  = data['max_idx_pred']
# metrics_true  = data['values_true']
# metrics_pred  = data['values_pred']
# loss_history  = data['loss_history']

#-------------------------------
# Section: Loss plot
#-------------------------------
# Section: Plot loss over time
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 5))
plt.plot(loss_history)
plt.title("Loss over time")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.grid(True)
plt.savefig(f'{result_path}/DataLoss_plot.png')

#-------------------------------
# Section: Expected vs predicted maximum node
#-------------------------------
import matplotlib.pyplot as plt

# Calculate accuracy
accuracy = np.mean(np.array(idx_max_true) == np.array(idx_max_pred))

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.scatter(idx_max_pred, idx_max_true, color='steelblue', s=80)
ax.text(0.95, 0.95, f'Accuracy: {accuracy:.2f}', transform=ax.transAxes,ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax.set_xlabel('Predicted max node', fontsize=13)
ax.set_ylabel('Expected max node', fontsize=13)
ax.set_title('Predicted maximum node', fontsize=15)
ax.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig(f'{result_path}/Indexes_plot.png',dpi=150)

#-------------------------------
# Section: Expected vs predicted values
#-------------------------------
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr

# Generate correlations
correlationP, pvalue = pearsonr(metrics_true, metrics_pred)
correlationS, pvalue = spearmanr(metrics_true, metrics_pred)

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.scatter(x = metrics_true, y = metrics_pred, color='steelblue', s=80)
ax.text(0.95, 0.95, 
    f'Pearson Correlation: {correlationP.item():.4f}\nSpearman Correlation: {correlationS.item():.4f}',
    transform=ax.transAxes, ha='right', va='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
)
ax.set_xlabel('Predicted values', fontsize=13)
ax.set_ylabel('Expected values', fontsize=13)
ax.set_title('Keystoneness values', fontsize=15)
ax.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig(f'{result_path}/Values_plot.png',dpi=150)

#-------------------------------
# Generate Summary
#-------------------------------
summary = f"""
Model Training Summary
=========================
Model: {model_declared}
Optimizer LR:   {optimizer.param_groups[0]['lr']}
Number of epochs: {epochs}
Seed: {n_seed}
Data path: {experiment_data}

Pearson Correlation:  {correlationP.item():.4f}  
Spearman Correlation: {correlationS.item():.4f}  
Maximum node accuracy: {accuracy.item():.4f}  
"""

with open(f'{result_path}/training_summary.txt', 'w') as f:
    f.write(summary)

print(summary)