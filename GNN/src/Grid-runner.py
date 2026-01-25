# =============================================================================
# GRID AND MODEL GENERATION
# =============================================================================
from sklearn.model_selection import ParameterGrid
import numpy as np 

# Define your parameter grid
param_grid = {
    'hidden_channels': [64],
    'num_layers': [2,5,10],
    'learning_rate': [.0001]
}

# Create all combinations
grid = ParameterGrid(param_grid)

# Import models to namespace
# --------------------------
from pathlib import Path

# List files
src_dir = Path('/home/mriveraceron/glv-research/gLV/GNN/architectures')
mods_files = src_dir.glob('Model-*.py')

for m in mods_files:
    print(f'>> Opening model: {m.name}')
    with open(m) as f:
        exec(f.read())

# Architecture list
# --------------------------
arch_config = {
    1: ('GATConv', ModelGATConv, {'h': 5, 'concat': False}),
    2: ('GCNConv', ModelGCNConv, {}),
    3: ('GraphConv', ModelGraphConv, {})
}

# Generate combinatorial models
# --------------------------
import pandas as pd

rows = []
COLUMNS = ['channels', 'layers', 'arch', 'tr_acc', 'tr_time', 'tr_NaN', 'val_acc', 'val_time']
model_list = []               # List to hold models
for params in grid:
    ch = params['hidden_channels']      # channels
    l = params['num_layers']            # layers
    for n in arch_config:
        arch_name, ModelClass, extra_kwargs = arch_config[n]                # architecture, model class, extra arguments
        m = ModelClass(hidden_channels=ch, num_layers=l, **extra_kwargs)    # Assign model
        model_list.append(m)                                                # Append to model list             
        rows.append({'channels': ch, 'layers': l, 'arch': arch_name, 'tr_acc': None, 'tr_time': None, 'tr_NaN': None, 'val_acc': None, 'val_time': None})

df = pd.DataFrame(rows)  # Convert to DataFrame
df

# =============================================================================
# TRAINING AND VALIDATION
# =============================================================================

import torch
import time

def training_fn(model, device, batched_paths, loss_fn, optimizer, epochs=100):
    model.train()
    loss_history  =  []         # Loss at epoch
    total_elapsed = 0           # Running time
    best_epoch = 0              # Best epoch
    best_loss = float('inf')    # Best loss
    for iter in range(1, epochs+1):
        start = time.time()
        epoch_loss = 0
        for path in batched_paths:
            data_list = torch.load(path, weights_only=False)          
            for data in data_list:
                data = data.to(device)
                optimizer.zero_grad()
                out = model(data)
                loss = loss_fn(out, data.y)
                loss.backward()
                #torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()
                epoch_loss += loss.item()   # Accumulate loss
        # Break loop if NaN encountered
        if torch.isnan(torch.tensor(epoch_loss)):
            print(f'>> NaN encountered at epoch {iter}. Stopping training.')
            flag = False    # Flag indicating NaN encountered
            return None, None, None, None, flag
        # Section: Best loss
        best_loss = float('inf') if iter == 1 else best_loss
        best_loss = min(best_loss, epoch_loss)
        best_epoch = iter if best_loss == epoch_loss else best_epoch
        # Append epoch loss to history
        loss_history.append(epoch_loss)
        elapsed = time.time() - start
        total_elapsed += elapsed
        # Print every 25 epochs
        if iter % 25 == 0:
            print(f"Epoch {iter}: Loss = {loss},  Elapsed time: {elapsed:.2f}")
    # Summary
    print(f'>> the total elapsed time with {epochs} epochs is {total_elapsed:.2f} seconds ( {total_elapsed/60:.2f} minutes)')    
    flag = True     # Flag indicating successful completion  
    return  loss_history, best_loss, best_epoch, total_elapsed, flag

# Declare validation function
# ---------------------------
from torch_geometric.data import Batch
import torch 
import time

def validation_fn(model, files, loss_fn, device):
    model.eval()
    total_elapsed = 0           # Running time
    start = time.time()
    val_loss, total_graphs = 0, 0
    true_idx, pred_idx = [], []
    with torch.no_grad():
        for path in files:
            data_list = torch.load(path, weights_only=False)
            total_graphs += len(data_list)
            # Batch all graphs together
            batch = Batch.from_data_list(data_list).to(device)
            # Forward pass
            out = model(batch)  # shape: [num_graphs * 30, features] or [num_graphs * 30]
            loss = loss_fn(out, batch.y)
            val_loss += loss.item()
            # Process each graph separately
            for i in range(batch.num_graphs):
                mask = batch.batch == i  # nodes belonging to graph i
                graph_y = batch.y[mask]  # 30 targets per graph
                graph_out = out[mask]    # 30 predictions per graph  
                # Get max index in 0-based
                true_idx.append(torch.argmax(graph_y).item() + 1)
                pred_idx.append(torch.argmax(graph_out).item() + 1)
    elapsed = time.time() - start
    total_elapsed += elapsed
    return true_idx, pred_idx, val_loss, total_graphs, total_elapsed


# Seeding function
#--------------------------
import random 

def seed_fn(seed=42):
    # Set ALL seeds for full reproducibility
    torch.manual_seed(seed)                 # Seed CPU 
    torch.cuda.manual_seed(seed)            # Seed GPU
    np.random.seed(seed)                    # Seed numpy
    random.seed(seed)                       # Seed python random
    torch.backends.cudnn.deterministic = True   # Ensure deterministic behavior
    torch.backends.cudnn.benchmark = False 

# Accuracy testing
# ---------------------------
def accuracy(true_idx, pred_idx):
    """Calculate accuracy between true and predicted indices."""
    correct = np.sum(true_idx == pred_idx)
    total = len(true_idx)
    acc_pct = (correct / total) * 100
    print(f'>> Validation Accuracy: {acc_pct:.2f}% ({correct}/{total})')
    return f'{correct}/{total}'

# =============================================================================
# LOAD BATCHES PATHS
# =============================================================================
import glob

dat_dir = Path('/home/mriveraceron/data/exp_20251125')
train_files = glob.glob(f'{dat_dir}/TrainBatch_*.pt')
valid_files = glob.glob(f'{dat_dir}/ValBatch_*.pt')


# Declare loss function and device
# --------------------------
import torch.nn as nn
import torch.optim as optim
import numpy as np

loss_fn = nn.MSELoss()
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
loss_fn = nn.MSELoss()
loss_dir = '/home/mriveraceron/glv-research/GNN-Logs/24Jan26'

for i, model in enumerate(model_list):
    print('-----------------------------\n')
    print(f'>> Training model {i} {model}')
    seed_fn(38)         # Set seed
    model = model.to(device)    # Copy model to device
    optimizer = optim.Adam(model.parameters(), lr=0.0001)
    epochs = 10
    # Calculate TRAINING dataset performance
    # --------------------------
    # returns: loss_history, best_loss, best_epoch, total_elapsed
    loss_history, _, _, total_elapsed, flag = training_fn(model, device, train_files, loss_fn, optimizer, epochs)
    # Skip if NaN encountered
    if flag == False:
        print (f'>> Model {i} encountered NaN during training. Skipping validation.')
        df.loc[i, 'tr_time'] = None
        df.loc[i, 'tr_acc'] = None
        df.loc[i, 'val_time'] = None
        df.loc[i, 'val_acc'] = None
        df.loc[i, 'tr_NaN'] = True
        continue    # Skip to next model
    # Testing 
    print(loss_history)
    # Continue if simulation completed successfully
    df.loc[i, 'tr_time'] = total_elapsed
    # returns: true_idx, pred_idx, val_loss, total_graphs, total_elapsed
    true_idx, pred_idx, val_loss, _, total_elapsed = validation_fn(model, train_files, loss_fn, device)
    acc_num = accuracy(true_idx, pred_idx)
    df.loc[i, 'tr_acc'] = acc_num
    # Calculate VALIDATION dataset performance
    # --------------------------
    true_idx, pred_idx, _, _, total_elapsed = validation_fn(model, valid_files, loss_fn, device)
    acc_num = accuracy(true_idx, pred_idx)
    df.loc[i, 'val_time'] = total_elapsed
    df.loc[i, 'val_acc'] = acc_num
    # Verify if NaNs are present
    # --------------------------
    if bool(pd.isna(loss_history).any()):
        df.loc[i, 'tr_NaN'] =  True
    else:
        df.loc[i, 'tr_NaN'] =  False
    # Save loss history
    file = f'{loss_dir}/loss_model-{i}.npy'  # Just create the path directly
    np.save(file, np.array(loss_history))
    print(f'>> Saved loss history to {file} \n')

# To load later
#loss_history = np.load('loss_history.npy').tolist()

# Save results dataframe
# --------------------------
file = f'{loss_dir}/results_dataframe.feather'
df.to_feather(file)
print(f'>> Saved results dataframe to {file}')
# To load later
# pd.read_feather(file)
