

#-------------------------------
# Section: Generate grid
#-------------------------------
from itertools import product
import pandas as pd

layers = [5,10]
LearnR = [1e-1, 1e-3, 1e-5, 1e-7]
channels = [64,128]
pairs = list(product(layers, LearnR, channels))
names = [f'Variant_{i}' for i in range(1, len(pairs)+1)]

# Create a datafrane
tuning_df = pd.DataFrame({
    'model_id': names,
    'channels': [l[2] for l in pairs],
    'layers': [l[0] for l in pairs],
    'learning_rate': [l[1] for l in pairs],
    'epochs': 700,
    'accuracy_idx': None,
    'pearson_corr': None,
    'spearman_corr': None
})

#-------------------------------
# Section: Declare model
#-------------------------------
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
# Section: Evaluation function at last epoch
#-------------------------------
def collect_metrics(loader, model_declared, device):
    idxt, idxp, mt, mp = [], [], [], []
    try:
        model_declared.eval()
        with torch.no_grad():
            for batch in loader:
                batch = batch.to(device)
                out = model_declared(batch)
                y_list = unbatch(batch.y, batch.batch)
                out_list = unbatch(out, batch.batch)
                for y, o in zip(y_list, out_list):
                    idxt.append(torch.argmax(y, dim=0))
                    idxp.append(torch.argmax(o, dim=0))
                    mt.append(y)
                    mp.append(o)
        # Convert to arrays
        idxt = torch.stack(idxt).cpu().numpy()
        idxp = torch.stack(idxp).cpu().numpy()
        mt = torch.cat(mt).cpu().numpy()
        mp = torch.cat(mp).cpu().numpy()
        return idxt, idxp, mt, mp
    finally:
        model_declared.train()
    

#-------------------------------
# Section: Training function 
#-------------------------------
import glob
import time
import random
import pandas as pd
from tqdm import tqdm
import os

from torch_geometric.utils import unbatch
from torch_geometric.loader import DataLoader

def training_DLloop(model_declared, device, all_data, weights_path, loss_fn, optimizer, epochs=100):
    #------------------------------------------
    model_declared.train()
    loss_history  =  []                   # Loss at epoch
    mt, mp, idxt, idxp = None, None, None, None
    total_elapsed = 0                   # Running time
    patience = np.floor(epochs * 0.1)   # Epochs to wait
    best_loss = float('inf')            # Best loss
    no_improve = 0
    is_last_epoch = False
    #---------------------
    # Section: Create batches of data
    loader = DataLoader(all_data, batch_size=30, shuffle=True)
    for epoch in tqdm(range(epochs), desc="Training"):
        start = time.time()
        epoch_loss = 0
        #--------------------------
        for batch in loader:
            #----------------------
            # Move it to device and run model
            # data = data_list[0]
            batch = batch.to(device)
            optimizer.zero_grad()
            out = model_declared(batch)
            loss = loss_fn(out, batch.y)
            loss.backward()
            #torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            epoch_loss += loss.item()   # Accumulate loss
            if torch.isnan(loss):
                tqdm.write(f"NaN loss detected at epoch {epoch}, stopping.")
                is_last_epoch = True
                break
        #----------------------
        # Append epoch loss to history
        #----------------------
        loss_history.append(epoch_loss)
        elapsed = time.time() - start
        total_elapsed += elapsed
        # Print every n epochs
        if epoch % 10 == 0:
            tqdm.write(f"Epoch {epoch}: Loss = {epoch_loss:.4f}, Elapsed time: {elapsed:.2f}")
        #----------------------
        # Section: At last epochs
        #----------------------
        # If it last epoch...
        if epoch == epochs - 1:
            is_last_epoch = True
        #----------------------
        # Section: Early stopping
        #----------------------
        if epoch_loss < best_loss:
            best_loss = epoch_loss
            no_improve = 0
        else :
            no_improve += 1
        if no_improve >= patience:
            tqdm.write(f"Early stopping at epoch {epoch}, no improvement for {patience} epochs.")
            is_last_epoch = True  
        #----------------------
        # Section: Evaluate model
        #----------------------
        if is_last_epoch:
            loss_history = np.array(loss_history)
            idxt, idxp, mt, mp = collect_metrics(loader, model_declared, device)
            # Save weights
            torch.save({
                'epoch': epoch,
                'model_state_dict': model_declared.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'loss': epoch_loss
            }, weights_path)
            break
    # Summary
    print(f'>> the total elapsed time with {epochs} epochs is {total_elapsed:.2f} seconds ( {total_elapsed/60:.2f} minutes)')   
    return  loss_history, mt, mp, idxt, idxp, total_elapsed

#-------------------------------
# Section: Seeding function
#-------------------------------
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

#------------------------------
# Section: Helper function to handle NaN
#-------------------------------
def fmt(value, decimals=4):
    val = value.item() if hasattr(value, 'item') else value
    if val != val:  # NaN check (NaN != NaN is always True)
        return 'NaN'
    return f'{val:.{decimals}f}'

#------------------------------
# Section: Runner function
#-------------------------------
import torch.optim as optim
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr

def wrapper(experiment_dir, row, all_data):
    #----------------------
    # Get variant information
    name = row.loc['model_id']
    lr = row.loc['learning_rate']
    layers = row.loc['layers']
    channels = row.loc['channels']
    epochs = row.loc['epochs']
    #---------------------
    # Set seed and run model
    loss_fn = nn.MSELoss()
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    n_seed = 42
    seed_fn(n_seed)
    model_declared = model(hidden_channels=int(channels), num_layers=int(layers)).to(device)
    optimizer = optim.Adam(model_declared.parameters(), lr=lr)
    print('Starting training of model:', name, '\n')
    print('The number of epochs will be:', epochs, '\n')
    print('The optimizer LR will be:', optimizer.param_groups[0]['lr'], '\n')
    #----------------------
    # Declare directory to save model weights
    model_results_dir = os.path.join(experiment_dir, name)
    os.makedirs(model_results_dir, exist_ok=True)
    weights_path = os.path.join(model_results_dir, f'{name}-weights.pth')
    loss_history, metrics_true, metrics_pred, idx_max_true, idx_max_pred, total_elapsed = training_DLloop(model_declared, device, all_data, weights_path, loss_fn, optimizer, epochs)
    # Save output
    np.savez(f'{model_results_dir}/{name}-values.npz',
        max_idx_true  = idx_max_true,
        max_idx_pred  = idx_max_pred,
        values_true   = metrics_true,
        values_pred   = metrics_pred,
        loss_history  = loss_history
    )
    #-------------------------------
    # Section: Loss plot
    #-------------------------------
    plt.figure(figsize=(10, 5))
    plt.plot(loss_history)
    plt.title("Loss over time")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.grid(True)
    plt.savefig(f'{model_results_dir}/DataLoss_plot.png')
    plt.close()
    #-------------------------------
    # Section: Expected vs predicted maximum node
    #-------------------------------
    # Calculate accuracy
    accuracy = np.mean(np.array(idx_max_true) == np.array(idx_max_pred))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(idx_max_pred, idx_max_true, color='steelblue', s=80)
    ax.text(0.95, 0.95, f'Accuracy: {accuracy:.2f}', transform=ax.transAxes,ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.set_xlabel('Predicted max node', fontsize=13)
    ax.set_ylabel('Expected max node', fontsize=13)
    ax.set_title('Predicted maximum node', fontsize=15)
    ax.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(f'{model_results_dir}/Indexes_plot.png',dpi=150)
    plt.close()
    #-------------------------------
    # Section: Expected vs predicted values
    #-------------------------------
    # Generate correlations
    correlationP, _ = pearsonr(metrics_true.flatten(), metrics_pred.flatten())
    correlationS, _ = spearmanr(metrics_true.flatten(), metrics_pred.flatten())
    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(x = metrics_pred, y = metrics_true, color='steelblue', s=80)
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
    plt.savefig(f'{model_results_dir}/Values_plot.png',dpi=150)
    plt.close()
    #-------------------------------
    # Generate Summary
    #-------------------------------
    summary = f"""
    Model Training Summary
    =========================
    Model variant: {name}
    Model: {model_declared}
    Optimizer LR:   {optimizer.param_groups[0]['lr']}
    Number of epochs: {epochs}
    Model layers: {layers}
    Model hidden channels: {channels}
    Seed: {n_seed}
    Data path: {experiment_data}
    -----------------------------------------------
    Pearson Correlation:  {fmt(correlationP)}  
    Spearman Correlation: {fmt(correlationS)}  
    Maximum node accuracy: {fmt(accuracy)}  
    Running time seconds: {total_elapsed}
    """
    with open(f'{model_results_dir}/training_summary.txt', 'w') as f:
        f.write(summary)
    print(summary)
    return accuracy, correlationP, correlationS

#-------------------------------
# Section: Generate data
#-------------------------------
import os
import glob

# Experiment data
tensors_dir = '/home/mriveraceron/glv-research/data_tensors/'
experiment_name = 'KBoost_v2_testing'
experiment_data = os.path.join(tensors_dir,experiment_name)
batched_paths = glob.glob(f"{experiment_data}/TrainBatch_*.pt")

all_data = []
for path in batched_paths:
    all_data.extend(torch.load(path, weights_only=False))

#-------------------------------
# Section: Create results directory
#-------------------------------
result_dir = '/home/mriveraceron/glv-research/tuning_results'
experiment_dir = os.path.join(result_dir, experiment_name)
os.makedirs(experiment_dir, exist_ok=True)
print('The experiment result directory will be:', experiment_dir, '\n')


row = tuning_df.iloc[0]
accuracy, correlationP, correlationS = wrapper(experiment_dir, row, all_data)

# In the runner loop
for i, row in tuning_df.iterrows():
    accuracy, corrP, corrS = wrapper(experiment_dir, row, all_data)
    tuning_df.loc[i, 'accuracy_idx'] = accuracy
    tuning_df.loc[i, 'pearson_corr'] = corrP
    tuning_df.loc[i, 'spearman_corr'] = corrS

tuning_df.to_csv(f'{experiment_dir}/tuning_results.csv', index=False)  # save after each run