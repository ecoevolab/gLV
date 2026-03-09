# SECTION: Define-GNN

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv, GCNConv, GATConv
import torch.optim as optim
import numpy as np
import random
import matplotlib.pyplot as plt
import numpy as np

class model(nn.Module):
    def __init__(self, hidden_channels=64, num_layers=5, dropout=0.5, heads = 4 ):
        super().__init__()
        self.convs = nn.ModuleList()
        self.dropout = dropout
        # First layer: 1 -> hidden_channels
        self.convs.append(GraphConv(1, hidden_channels))
        # self.convs.append(GATConv(1, hidden_channels, heads=heads))
        # Middle layers: hidden_channels -> hidden_channels
        for _ in range(num_layers - 2):
            #self.convs.append(GATConv(hidden_channels*heads, hidden_channels, heads=heads))
            self.convs.append(GraphConv(hidden_channels, hidden_channels))
        # Last layer: hidden_channels -> 1
        # self.convs.append(GATConv(hidden_channels*heads, 1, heads=heads))
        self.convs.append(GraphConv(hidden_channels, 1))
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_weights
        # Apply all layers except the last
        for i, conv in enumerate(self.convs[:-1]):
            x = conv(x, edge_index, edge_weight)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout)
        # Apply last layer with sigmoid
        x = self.convs[-1](x, edge_index, edge_weight)
        # x = torch.sigmoid(x)
        return x  # [num_nodes]


def seed_fn(seed=42):
    # Set ALL seeds for full reproducibility
    torch.manual_seed(seed)                 # Seed CPU 
    torch.cuda.manual_seed(seed)            # Seed GPU
    np.random.seed(seed)                    # Seed numpy
    random.seed(seed)                       # Seed python random
    torch.backends.cudnn.deterministic = True   # Ensure deterministic behavior
    torch.backends.cudnn.benchmark = False 


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

#============================================================
import glob
from torch_geometric.loader import DataLoader
import time
import random
import torch 


def training_loop(model, device, batched_paths, loss_fn, optimizer, epochs=100):
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
    return  loss_history, best_loss, best_epoch


#============================================================
from torch_geometric.data import Batch

def validation_fn(model, files, loss_fn, device):
    model.eval()
    # Assign variables
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
    return true_idx, pred_idx, val_loss, total_graphs


# =========================
# Model testing
def accuracy(true, preds):
    true = true_idx
    preds = pred_idx
    result = [a == b for a, b in zip(true, preds)]
    correct = sum(result)
    accuracy = (correct / len(true)) * 100
    return(f'>> Validation Accuracy: {accuracy:.2f}% ({correct}/{len(true)})')

# Declare data paths
directory = '/home/mriveraceron/data/exp_20251125'
train_files = glob.glob(f'{directory}/TrainBatch_*.pt')
valid_files = glob.glob(f'{directory}/ValBatch_*.pt')

# Loss function and device
loss_fn = nn.MSELoss()
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Set seed and run model
seed_fn(38)
mymodel = model(hidden_channels=64, num_layers=5, dropout=0.5, heads=4).to(device)
optimizer = optim.Adam(mymodel.parameters(), lr=0.00001)
epochs = 200

mod1_loss, mod1_BestLoss, mod1_BestEpoch = training_loop(mymodel,device, train_files, loss_fn, optimizer, epochs)
true_idx, pred_idx, mod1_ValLoss, total_graphs = validation_fn(mymodel, train_files, loss_fn, device)
pf1 = accuracy(true_idx, pred_idx)
print(pf1)

# Plot the loss over epochs
print(mod1_loss)
fig = loss_plotter(mod1_loss, epochs)
fig.savefig('/home/mriveraceron/glv-research/plots/9Dec-Test-Loss-V4.png', dpi=150, bbox_inches='tight')