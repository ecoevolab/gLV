
#---------------------------------------
# Section: Set seeds
import torch 
import numpy as np 
import random 

# Define function for seeding
def set_seed(seed=42):
    # Set ALL seeds for full reproducibility
    torch.manual_seed(seed)                     # Seed CPU 
    torch.cuda.manual_seed(seed)                # Seed GPU
    np.random.seed(seed)                        # Seed numpy
    random.seed(seed)                           # Seed python random
    torch.backends.cudnn.deterministic = True   # Ensure deterministic behavior
    torch.backends.cudnn.benchmark = False   

set_seed(seed=54)  # Ensure reproducibility

#---------------------------------------
# Section: Load data
# Load data
import glob

exp = "4379fd40-9f0a"
dpath = "/home/mriveraceron/data/4379fd40-9f0a"
train_path = glob.glob(f'{dpath}/TrainBatch_*.pt')
val_path = glob.glob(f'{dpath}/ValBatch_*.pt')

#---------------------------------------
# SECTION: Define-GNN
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv
import torch.optim as optim

class simple_gnn_gcn(nn.Module):
    def __init__(self, hidden_channels=64):
        super().__init__()
        self.conv1 = GraphConv(1, hidden_channels)
        self.conv2 = GraphConv(hidden_channels, 1)
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_weights
        x = self.conv1(x, edge_index, edge_weight)
        x = F.relu(x)
        x = self.conv2(x, edge_index, edge_weight)
        x = torch.sigmoid(x)  # Outputs between 0-1
        return x  # [num_nodes]
    
# Declare other model parameters
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = simple_gnn_gcn(hidden_channels=72).to(device)
loss_fn = nn.MSELoss()                                                # Loss function for regression
optimizer = optim.Adam(model.parameters(), lr=0.01) 

#------------------------------------------
# Section: Training loop
from tqdm import tqdm
import time
import torch
from torch_geometric.loader import DataLoader

def train_fn(epochs, model, train_path):
    start = time.time()
    model.train()
    # Empty lists for predictions, targets, loss at each epoch
    x_train, y_train, loss_epochs  = [], [], []
    for i in tqdm(range(1, epochs+1), total=epochs, desc="Training model:"):
        epoch_loss = 0  # Reset for each epoch
        num_batches = 0  # Track number of batches
        for path in train_path:
            data = torch.load(path, weights_only = False)
            Dloaded = DataLoader(data, batch_size=200, shuffle=True)
            # Iterate over all data loaded
            for data in Dloaded:
                optimizer.zero_grad()
                data = data.to(device)
                out = model(data)
                loss = loss_fn(out, data.y)
                loss.backward()
                optimizer.step()
                epoch_loss += loss.item()   # Accumulate loss
                num_batches += 1
        # Store average loss per batch for this epoch
        avg_loss = epoch_loss / num_batches
        loss_epochs.append(avg_loss)
        if i==(epochs):
            x_train.append(out.cpu().detach().numpy()) 
            y_train.append(data.y.cpu().detach().numpy())
            path = f'/home/mriveraceron/glv-research/GNN-params/{exp}/'
        elapsed_time = time.time() - start
        print(f">> Elapsed time for training epoch {i}: {elapsed_time:.2f}")
    # Elapsed time for training 
    elapsed_time = time.time() - start
    print(f">> Elapsed time for training: {elapsed_time:.2f}")
    return x_train, y_train, loss_epochs

epochs = 100
x_train, y_train, loss_epochs = train_fn(epochs = epochs, model = model, train_path = train_path)
loss = [ '%.5f' % elem for elem in loss_epochs ]
#----------------------------------------------------------
# SECTION: Validation-loop
def val_fn(model, paths):
    total_loss = 0
    model.eval()  # Set to evaluation mode
    start_time = time.time()
    x_val, y_val = [], []
    # Disable gradient computation
    with torch.no_grad():  
        for path in paths:
                    data = torch.load(path, weights_only = False)
                    Dloaded = DataLoader(data, batch_size=100, shuffle=True)
                    for data in Dloaded:
                        optimizer.zero_grad()
                        data = data.to(device)
                        out = model(data)
                        loss = loss_fn(out, data.y)
                        total_loss += loss.item()   # Accumulate loss
                        x_val.append(out.cpu().numpy())
                        y_val.append(data.y.cpu().numpy())
    print(f"Validation Loss = {total_loss:.4f}")
    elapsed_time = time.time() - start_time
    print(f"Validation completed in {elapsed_time:.2f} seconds")
    return x_val, y_val

x_val, y_val = val_fn(model= model, paths = val_path)


#-----------------------------------------------------
# Section: Loss vs epochs
import matplotlib.pyplot as plt
import numpy as np

def loss_plt(loss = None, epochs = None, path = None ):
    # After collecting your data
    y = np.array(loss, dtype=float)
    x = list(range(1, epochs + 1))
    # Create scatter plot
    fig = plt.figure(figsize=(8, 8))
    plt.plot(x, y, alpha=0.8)
    # Add labels
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title('Loss over epochs')
    plt.grid(True, alpha=0.3)
    # Set plot axis limits
    plt.ylim(min(y) * 0.9, max(y) * 1.1)  # Add 10% padding
    plt.xlim(0, max(x) + 1)
    plt.tight_layout()
    return fig

fig = loss_plt(loss = loss, epochs = epochs)
fig.savefig('/home/mriveraceron/glv-research/plots/10-Oct-Loss.jpg', dpi=150, bbox_inches='tight')

#-----------------------------------------------------
# Section: Preds vs Tgts
import matplotlib.pyplot as plt

def preds_plt(preds = None, tgts = None):
    # Create scatter plot
    fig = plt.figure(figsize=(8, 8))
    x = np.concatenate(preds)   # predictions
    y = np.concatenate(tgts)    # True values
    plt.scatter(x, y, alpha=0.5)
    # Add perfect prediction line (y=x)
    max_val = max(x.max(), y.max())
    plt.plot([0, max_val], [0, max_val], 'r--', label='Perfect prediction')
    # Add labels
    plt.xlabel('Predictions')
    plt.ylabel('True Values')
    plt.title('Predictions vs True Values')
    plt.legend()
    plt.grid(True, alpha=0.3)
    # Set plot axis limits
    plt.xlim(x.min(), y.max())
    plt.ylim(x.min(), y.max())
    plt.tight_layout()
    return fig

fig = preds_plt(preds = x_train, tgts = y_train)
fig.savefig('/home/mriveraceron/glv-research/plots/10-Oct-preds.jpg', dpi=150, bbox_inches='tight')