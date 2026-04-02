# Training script for random networks.
# 


#-------------------------------
# Section: Generate grid
#-------------------------------
"""
    We will be training model with random networks.
    Model to train will be Variant 3 at different sizes.
"""
from itertools import product
import pandas as pd

size = [1000,5000,8000]
names = [f'Variant_{i}' for i in range(1, len(size)+1)]

# Create a datafrane
tuning_df = pd.DataFrame({
    'model_id': names,
    'size': size,
    'channels': 64,
    'layers': 5,
    'learning_rate': 1e-03,
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
from collections import namedtuple
from scipy.stats import pearsonr
from scipy.stats import spearmanr

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
        idxt = torch.cat(idxt).cpu().numpy()
        idxp = torch.cat(idxp).cpu().numpy()
        mt = torch.cat(mt).cpu().numpy()
        mp = torch.cat(mp).cpu().numpy()
        Result = namedtuple('Result', ['idxt', 'idxp', 'mt', 'mp'])
        return Result(idxt, idxp, mt, mp)
    finally:
        model_declared.train()

def compute_metrics(metrics_list):
    idxt, idxp = metrics_list.idxt, metrics_list.idxp
    mt, mp = metrics_list.mt, metrics_list.mp
    accuracy = np.mean(np.array(idxt) == np.array(idxp))
    if np.std(mt) == 0 or np.std(mp) == 0:
        print("Cannot compute correlation: one input is constant.")
        correlationP = correlationS = float('nan')
    else:
        correlationP, _ = pearsonr(mt.flatten(), mp.flatten())
        correlationS, _ = spearmanr(mt.flatten(), mp.flatten())
    Result = namedtuple('Result', ['acc', 'corrP', 'corrS'])
    return Result(accuracy, correlationP, correlationS)

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

def training_DLloop(model_declared, device, data_train, data_eval, weights_path, loss_fn, optimizer, epochs):
    #------------------------------------------
    model_declared.train()
    loss_history  =  []                 # Loss at epoch
    total_elapsed = 0                   # Running time
    is_last_epoch = False
    # Early stopping
    patience = np.floor(epochs * 0.2)   # Epochs to wait
    best_loss = float('inf')            # Best loss
    no_improve = 0
    #---------------------
    # Section: Create batches of data
    loader_train = DataLoader(data_train, batch_size=30, shuffle=False)
    loader_eval = DataLoader(data_eval, batch_size=30, shuffle=False)
    for epoch in tqdm(range(epochs), desc="Training"):
        start = time.time()
        epoch_loss = 0
        #--------------------------
        for batch in loader_train:
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
        loss_history = np.append(loss_history, epoch_loss)
        elapsed = time.time() - start
        total_elapsed += elapsed
        # Print every n epochs
        if epoch % 10 == 0:
            tqdm.write(f"Epoch {epoch}: Loss = {epoch_loss}, Elapsed time: {elapsed:.2f}")
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
            metrics_list = collect_metrics(loader_eval, model_declared, device)
            performance_list = compute_metrics(metrics_list)
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
    return  loss_history, metrics_list, performance_list, total_elapsed

#-------------------------------
# Section: Function to save Summary
#-------------------------------
def summarize(model, optimizer, row, epochs_runned, val_samples, n_seed, data_path, performance_list, result_exp_dir, total_elapsed):
    summary = f"""
    Model Training Summary
    =========================
    Model variant: {row['model_id']}
    Model: {model}
    Samples for training {row['size']}
    Optimizer LR:   {optimizer.param_groups[0]['lr']}
    Number of epochs: {row['epochs']}
    Model layers: {row['layers']}
    Model hidden channels: {row['channels']}
    Seed: {n_seed}
    Training data path: {data_path}
    -----------------------------------------------
    Validation data path: {data_path}
    Validation samples: {val_samples}
    Pearson Correlation:  {performance_list.corrP}    
    Spearman Correlation: {performance_list.corrS}   
    Maximum node accuracy: {performance_list.acc}  
    Running time seconds: {total_elapsed}
    Epochs performed: {epochs_runned}
    """
    with open(f'{result_exp_dir}/training_summary.txt', 'w') as f:
        f.write(summary)
    print(summary)


#-------------------------------
# Section: Run model
#-------------------------------
import os
import glob

# Add values to table
for i, row in tuning_df.iterrows():
    #-------------------------
    # Declare model parameters
    #-------------------------
    # Hyperparameters
    # row = tuning_df.iloc[0]
    size = row['size']
    lr = row['learning_rate']
    model_name = row['model_id']
    epochs = row['epochs']
    # Seeding function
    n_seed = 42
    seed_fn(seed=n_seed)
    # Model parameters
    loss_fn = nn.MSELoss()
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model_declared = model(hidden_channels=64, num_layers=5).to(device)
    optimizer = optim.Adam(model_declared.parameters(), lr=lr)
    #-------------------------
    # Generate data
    #-------------------------
    data_dir = '/home/mriveraceron/glv-research/data_null/d2f93775a813'
    train_paths =  glob.glob(f'{data_dir}/TrainBatch_*.pt')
    eval_paths = glob.glob(f'{data_dir}/ValBatch_*.pt')
    data_train, data_eval = [] , []
    for path in train_paths:
        data_train.extend(torch.load(path, weights_only=False))
    for path in eval_paths:
        data_eval.extend(torch.load(path, weights_only=False))
    # Create result directory
    results_dir = '/home/mriveraceron/glv-research/tuning_results' 
    result_exp_dir = f'{results_dir}/{os.path.basename(data_dir)}/{model_name}'
    os.makedirs(result_exp_dir, exist_ok=True)
    #-------------------------
    # Run model
    #-------------------------
    weights_path = f'{result_exp_dir}/model_weights.pth'
    loss_history, metrics_list, performance_list, total_elapsed = training_DLloop(model_declared, device, data_train[:size], data_eval, weights_path, loss_fn, optimizer, epochs)
    epochs_runned = len(loss_history)
    val_samples = len(data_eval)
    summarize(model_declared, optimizer, row, epochs_runned, n_seed, data_dir, performance_list, result_exp_dir, total_elapsed)
    # Save metrics result_exp_dir
    np.savez(f'{result_exp_dir}/metric-values.npz',
        max_idx_true  = metrics_list.idxt,
        max_idx_pred  = metrics_list.idxp,
        values_true   = metrics_list.mt,
        values_pred   = metrics_list.mp,
        loss_history  = loss_history
    )
    # Add validation to dataframe
    tuning_df.loc[i, 'accuracy_idx'] = performance_list.acc
    tuning_df.loc[i, 'pearson_corr'] = performance_list.corrP
    tuning_df.loc[i, 'spearman_corr'] = performance_list.corrS