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
    


#-------------------------------
# Section: Traning loop
import glob
import time
import random
import pandas as pd
from tqdm import tqdm
import os

def training_loop(model_declared, device, batched_paths, checkpoint_path, weights_path, loss_fn, optimizer, epochs=100):
    #------------------------------------------
    # Section: Load checkpoint
    checkpoint = torch.load(checkpoint_path, weights_only=False)
    model_declared.load_state_dict(checkpoint)
    # optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
    # epoch_start = checkpoint['epoch']
    # loss = checkpoint['loss']
    #------------------------------------------
    model_declared.train()
    loss_history  =  []             # Loss at epoch
    total_elapsed = 0               # Running time
    metrics_pred, metrics_true = [], []   # Declare list for metrics
    idx_max_true, idx_max_pred = [], []     # Declare list for max values
    for i, iter in enumerate(tqdm(range(epochs), desc="Training")):
        start = time.time()
        epoch_loss = 0
        for path in batched_paths:
            # path = batched_paths[1]
            data_list = torch.load(path, weights_only=False)          
            for data in data_list:
                #----------------------
                # Section: Move it to device and run model
                # data = data_list[0]
                data = data.to(device)
                optimizer.zero_grad()
                out = model_declared(data)
                loss = loss_fn(out, data.y)
                loss.backward()
                #----------------------
                # Section: Last layeer information
                if iter == epochs - 1:  # Only for the last epoch
                    #----------------------
                    # Append true values
                    max_data = torch.argmax(data.y, dim=0).detach().cpu().numpy()
                    idx_max_true.append(max_data)
                    metrics_true.append(data.y[:, 0].detach().cpu().numpy())
                    #----------------------
                    # Append predicted  values
                    max_out = torch.argmax(out, dim=0).detach().cpu().numpy()
                    idx_max_pred.append(max_out)
                    metrics_pred.append(out[:, 0].detach().cpu().numpy())
                    #----------------------
                    # Save weights
                    torch.save({
                        'epoch': iter,
                        'model_state_dict': model_declared.state_dict(),
                        'optimizer_state_dict': optimizer.state_dict(),
                        'loss': loss
                    }, weights_path)
                #torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()
                epoch_loss += loss.item()   # Accumulate loss
        # Append epoch loss to history
        loss_history.append(epoch_loss)
        elapsed = time.time() - start
        total_elapsed += elapsed
        # Print every 25 epochs
        if iter % 10 == 0:
            tqdm.write(f"Epoch {iter}: Loss = {epoch_loss:.4f}, Elapsed time: {elapsed:.2f}")
    # Summary
    print(f'>> the total elapsed time with {epochs} epochs is {total_elapsed:.2f} seconds ( {total_elapsed/60:.2f} minutes)')   
    return  loss_history, metrics_true, metrics_pred, idx_max_true, idx_max_pred




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
optimizer = optim.Adam(model_declared.parameters(), lr=0.001)
epochs = 400
print('The number of epochs will be:', epochs, '\n')

# Load file paths
data_dir = '/home/mriveraceron/glv-research/data_tensors/Boosted_filtered'

# Weights to keep training
checkpoint_path = '/home/mriveraceron/glv-research/model_weights/Boosted_filtered_V2.pth'
print('Weights checkpoint path:', checkpoint_path, '\n')

# New weights file
weights_path = '/home/mriveraceron/glv-research/model_weights/Boosted_filtered_V3.pth'
batched_paths = glob.glob(f"{data_dir}/*.pt")
loss_history, metrics_true, metrics_pred, idx_max_true, idx_max_pred = training_loop(model_declared, device, batched_paths, checkpoint_path, weights_path, loss_fn, optimizer, epochs)

# Flatten lists of indexes
idx_max_true = np.concatenate(idx_max_true).tolist()
idx_max_pred = np.concatenate(idx_max_pred).tolist()
#  Flatten lists of values
metrics_true = np.concatenate(metrics_true).tolist()
metrics_pred = np.concatenate(metrics_pred).tolist()

#-------------------------------
# Section: Save loss and max data/out tensor
result_path = '/home/mriveraceron/glv-research/Results/Boosted_keystone/Filtered_AllFeats_V3'
os.makedirs(result_path, exist_ok=True)
print('The results directory will be:', result_path, '\n')
#-------------------------------
# Save max indexes
np.save(f'{result_path}/max_idx_true.npy', idx_max_true)
np.save(f'{result_path}/max_idx_pred.npy', idx_max_pred)
#-------------------------------
# Save metrics
np.save(f'{result_path}/values_true.npy', metrics_true)
np.save(f'{result_path}/values_pred.npy', metrics_pred)
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
#plt.xlim(400, 800)
plt.grid(True)
plt.savefig(f'{result_path}/DataLoss_plot.png')

#-------------------------------
# Section: Plot indexes predicted vs true
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
# Section: Plot indexes predicted vs true
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
    f'Pearson Correlation: {correlationP:.4f}\nSpearman Correlation: {correlationS:.4f}',
    transform=ax.transAxes, ha='right', va='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
)
ax.set_xlabel('Predicted values', fontsize=13)
ax.set_ylabel('Expected values', fontsize=13)
ax.set_title('Keystoneness values', fontsize=15)
ax.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig(f'{result_path}/Values_plot.png',dpi=150)