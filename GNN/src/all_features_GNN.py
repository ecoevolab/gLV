# 24-February-2024
# This code is for running the model using all node features.

#-------------------------------
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
        self.convs.append(GraphConv(hidden_channels, 3))
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
    

#-------------------------------
# Section: Traning loop
import glob
import time
import random
import torch 
import pandas as pd
from tqdm import tqdm

def training_loop(model_declared, device, batched_paths, loss_fn, optimizer, epochs=100):
    model_declared.train()
    loss_history  =  []             # Loss at epoch
    total_elapsed = 0               # Running time
    rows_pred, rows_true = [], []   # Declare list for metrics
    max_true, max_pred = [], []     # Declare list for max values
    for i, iter in enumerate(tqdm(range(epochs), desc="Training")):
        start = time.time()
        epoch_loss = 0
        for path in batched_paths:
            # path = batched_paths[1]
            data_list = torch.load(path, weights_only=False)          
            for data in data_list:
                #----------------------
                # Move it to device and run model
                data = data.to(device)
                optimizer.zero_grad()
                out = model_declared(data)
                loss = loss_fn(out, data.y)
                loss.backward()
                #----------------------
                # Section: Extract node with maximum values
                if iter == epochs - 1:  # Only for the last epoch
                    max_data = torch.argmax(data.y, dim=0).detach().cpu().numpy()
                    max_true.append(max_data)
                    # prop_extinctions, dissimilarity, keystoneness
                    rows_true.append({
                        "prop_extinctions": data.y[:, 0].detach().cpu().numpy(),
                        "dissimilarity": data.y[:, 1].detach().cpu().numpy(),
                        "keystoneness": data.y[:, 2].detach().cpu().numpy()
                        })
                    #----------------------
                    # Append predicted  values
                    max_out = torch.argmax(out, dim=0).detach().cpu().numpy()
                    max_pred.append(max_out)
                    rows_pred.append({
                        "prop_extinctions": out[:, 0].detach().cpu().numpy(),
                        "dissimilarity": out[:, 1].detach().cpu().numpy(),
                        "keystoneness": out[:, 2].detach().cpu().numpy()
                        })
                #torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()
                epoch_loss += loss.item()   # Accumulate loss
        # Append epoch loss to history
        loss_history.append(epoch_loss)
        elapsed = time.time() - start
        total_elapsed += elapsed
        # Print every 25 epochs
        if iter % 10 == 0:
            tqdm.write(f"Epoch {iter}: Loss = {loss:.4f}, Elapsed time: {elapsed:.2f}")
    # Summary
    print(f'>> the total elapsed time with {epochs} epochs is {total_elapsed:.2f} seconds ( {total_elapsed/60:.2f} minutes)')   
    #-----------------------
    # Section: Convert to numpy arrays
    df_pred = pd.concat([ pd.DataFrame(row).assign(graph_id=i) for i, row in enumerate(rows_pred)], ignore_index=True)
    df_true = pd.concat([ pd.DataFrame(row).assign(graph_id=i) for i, row in enumerate(rows_true)], ignore_index=True)  
    # Convert arrays to data frames
    names=['prop_extinctions', 'dissimilarity', 'keystoneness']
    idx_max_true = pd.DataFrame(max_true, columns=names)
    idx_max_pred = pd.DataFrame(max_pred, columns=names)
    return  loss_history, df_true, df_pred, idx_max_true, idx_max_pred


#-------------------------------
# Section: Run model
import torch.optim as optim
import glob
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
seed_fn(42)
model_declared = model(hidden_channels=64, num_layers=5).to(device)
optimizer = optim.Adam(model_declared.parameters(), lr=0.00001)
epochs = 300

# Load file paths
data_dir = '/home/mriveraceron/data/Boosted_keystone'
batched_paths = glob.glob(f"{data_dir}/*.pt")
loss_history, df_true, df_pred, idx_max_true, idx_max_pred= training_loop(model_declared, device, batched_paths, loss_fn, optimizer, epochs)


#-------------------------------
# Section: Save loss and max data/out tensors
import os
results_dir = '/home/mriveraceron/Results'
result_path = f'{results_dir}/{os.path.basename(data_dir)}/Raw-AllFeatures'
os.makedirs(result_path, exist_ok=True)
#-------------------------------
# Save max indexes
idx_max_true.to_feather(f'{result_path}/max_idx_true.feather')
idx_max_pred.to_feather(f'{result_path}/max_idx_pred.feather')
#-------------------------------
# Save metrics
df_true.to_feather(f'{result_path}/Metrics_true.feather')
df_pred.to_feather(f'{result_path}/Metrics_pred.feather')
#-------------------------------
# Save loss history
np.save(f'{result_path}/loss_history.npy', np.array(loss_history))
# np.load(f'{result_path}/Dummy_loss.npy')
# x= np.load(f'{result_path}/Dummy_max_true.npy')


#-------------------------------
# Section: Plot loss over time
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 5))
plt.plot(loss_history)
plt.title("Loss over time")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.grid(True)
plt.savefig(f'{result_path}/Data_loss_plot.png')

#-------------------------------
# Section: Plot values predicted vs true
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
# Plot prop_extinctions
axes[0].scatter(x=df_true['prop_extinctions'], y=df_pred['prop_extinctions'], color='blue')
axes[0].set_title('prop_extinctions True vs predicted')
# Plot dissimilarity
axes[1].scatter(x=df_true['dissimilarity'], y=df_pred['dissimilarity'], color='blue')
axes[1].set_title('dissimilarity True vs predicted')
# Plot keystoneness
axes[2].scatter(x=df_true['keystoneness'], y=df_pred['keystoneness'], color='blue')
axes[2].set_title('keystoneness True vs predicted')
for ax in axes:
    ax.grid(True)
    ax.set_xlabel('True values Raw Data')
    ax.set_ylabel('Predicted values')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)  
    
fig.savefig(f'{result_path}/Values_TP.png')

#-------------------------------
# Section: Plot indexes of max values predicted vs true
v = []
for c in ['prop_extinctions', 'dissimilarity', 'keystoneness']:
    acc = (idx_max_true[c] == idx_max_pred[c]).mean()
    print(f"Accuracy for {c}: {acc:.2f}")
    v.append(acc)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
# Plot prop_extinctions 
axes[0].scatter(x=idx_max_true['prop_extinctions'], y=idx_max_pred['prop_extinctions'], color='blue')
axes[0].text(0.95, 0.95, f'Accuracy: {v[0]:.2f}', transform=axes[0].transAxes,ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
axes[0].set_title('prop_extinctions Max Index True vs predicted')
# Plot dissimilarity
axes[1].scatter(x=idx_max_true['dissimilarity'], y=idx_max_pred['dissimilarity'], color='blue')
axes[1].text(0.95, 0.95, f'Accuracy: {v[1]:.2f}', transform=axes[1].transAxes,ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
axes[1].set_title('dissimilarity Max Index True vs predicted')
# Plot keystoneness
axes[2].scatter(x=idx_max_true['keystoneness'], y=idx_max_pred['keystoneness'], color='blue')
axes[2].text(0.95, 0.95, f'Accuracy: {v[2]:.2f}', transform=axes[2].transAxes,ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
axes[2].set_title('keystoneness Max Index True vs predicted')
for ax in axes:
    ax.grid(True)
    ax.set_xlabel('True Max Index')
    ax.set_ylabel('Predicted Max Index')
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 30)  


fig.savefig(f'{result_path}/Values_MaxIndex_TP.png')
