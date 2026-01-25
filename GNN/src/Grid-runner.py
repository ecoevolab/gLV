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
    1: ('GATConv', ModelGATConv, {'h': 5}),
    2: ('GCNConv', ModelGCNConv, {}),
    3: ('GraphConv', ModelGraphConv, {})
}

# Generate combinatorial models
# --------------------------
import pandas as pd

rows = []
COLUMNS = ['channels', 'layers', 'arch', 'tr_acc', 'tr_time', 'train_NaN', 'val_acc', 'val_time']
model_list = []               # List to hold models
for params in grid:
    ch = params['hidden_channels']      # channels
    l = params['num_layers']            # layers
    for n in arch_config:
        arch_name, ModelClass, extra_kwargs = arch_config[n]                # architecture, model class, extra arguments
        m = ModelClass(hidden_channels=ch, num_layers=l, **extra_kwargs)    # Assign model
        model_list.append(m)                                                # Append to model list             
        rows.append({'channels': ch, 'layers': l, 'arch': arch_name, 'tr_acc': None, 'tr_time': None, 'train_NaN': None, 'val_acc': None, 'val_time': None})

df = pd.DataFrame(rows)  # Convert to DataFrame
df

# =============================================================================
# TRAINING AND VALIDATION
# =============================================================================

src_dir = Path('/home/mriveraceron/glv-research/gLV/GNN/architectures')
mods_files = src_dir.glob('training_fun.py')

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

# =============================================================================
# LOAD BATCHES PATHS
# =============================================================================
import glob

dat_dir = Path('/home/mriveraceron/data/exp_20251125')
train_files = glob.glob(f'{dat_dir}/TrainBatch_*.pt')
valid_files = glob.glob(f'{dat_dir}/ValBatch_*.pt')

# =============================================================================
# PLOTTERS AND METRICS
# =============================================================================
# Loss plotter
# ---------------------------
import matplotlib.pt as plt
import numpy as np 

def loss_plotter(loss_epochs = None, epochs = None):
    # After collecting your data
    y = np.round(loss_epochs, 10)
    x = list(range(0, epochs))
    # Create scatter plot
    fig = plt.figure(figsize=(8, 8))
    # y = np.log1p(y)  # Log scale for better visualization
    plt.plot(x, y, alpha=0.5)
    # Add labels and title
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title('Loss over epochs')
    plt.grid(True, alpha=0.3)
    plt.tight_layout() 
    plt.ylim(0, max(y))
    plt.xlim(0, max(x))
    return fig

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
import numpy as np

def accuracy(true_idx, pred_idx):
    """Calculate accuracy between true and predicted indices."""
    correct = np.sum(true_idx == pred_idx)
    total = len(true_idx)
    acc_pct = (correct / total) * 100
    print(f'>> Validation Accuracy: {acc_pct:.2f}% ({correct}/{total})')
    return f'{correct}/{total}'

# Load batches paths
# --------------------------
directory = '/home/mriveraceron/data/exp_20251125'
train_files = glob.glob(f'{directory}/TrainBatch_*.pt')
valid_files = glob.glob(f'{directory}/ValBatch_*.pt')

# Declare loss function and device
# --------------------------
import torch.nn as nn

loss_fn = nn.MSELoss()
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Run model
import torch.optim as optim

for i,model in enumerate(model_list):
    seed_fn(38)         # Set seed
    test_model = model.to(device)
    optimizer = optim.Adam(test_model.parameters(), lr=0.00001)
    epochs = 200
    loss_history, _, _, total_elapsed = training_fn(test_model, device, train_files, loss_fn, optimizer, epochs)
    df.loc[i, 'train_time'] = total_elapsed
    # Calculate TRAINING dataset performance
    # --------------------------
    true_idx, pred_idx, val_loss, total_graphs = validation_fn(test_model, train_files, loss_fn, device)
    acc_num = accuracy(true_idx, pred_idx)
    df.loc[i, 'train_acc'] = acc_num
    # Calculate VALIDATION dataset performance
    # --------------------------
    true_idx, pred_idx, _, _, total_elapsed = validation_fn(test_model, valid_files, loss_fn, device)
    acc_num = accuracy(true_idx, pred_idx)
    df.loc[i, 'val_time'] = total_elapsed
    df.loc[i, 'val_acc'] = acc_num
    # Verify if NaNs are present
    # --------------------------
    if bool(pd.isna(loss_history).any()):
        df.loc[i, 'train_NaN'] =  True
    else:
        df.loc[i, 'train_NaN'] =  False
    # Plot the loss over epochs
    # --------------------------
    fig = loss_plotter(loss_history, epochs)
    name = f'/home/mriveraceron/glv-research/plots/Mod{i}-training_loss.png'
    fig.savefig(name, dpi=150, bbox_inches='tight')
