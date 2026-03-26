


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
# Section: Evaluation function 
#-------------------------------
from torch_geometric.utils import unbatch
from datetime import datetime
import time

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


def compute_metrics(idxt, idxp, mt, mp):
    accuracy = np.mean(np.array(idxt) == np.array(idxp))
    if np.std(mt) == 0 or np.std(mp) == 0:
        print("Cannot compute correlation: one input is constant.")
        correlationP = correlationS = float('nan')
    else:
        correlationP, _ = pearsonr(mt.flatten(), mp.flatten())
        correlationS, _ = spearmanr(mt.flatten(), mp.flatten())
    return accuracy, correlationP, correlationS

def save_eval_summary(results_dir, model_dir, data_dir, accuracy, correlationP, correlationS, elapsed, loader):
    summary = (
        f"Evaluation Summary\n"
        f"{'=' * 40}\n"
        f"Date:                  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"Model evaluated:       {model_dir}\n"
        f"Data used for evaluation:  {data_dir}\n"
        f"Results dir:           {results_dir}\n"
        f"{'=' * 40}\n"
        f"Dataset\n"
        f"  Samples evaluated:   {len(loader.dataset)}\n"
        f"  Batch size:          {loader.batch_size}\n"
        f"{'=' * 40}\n"
        f"Metrics\n"
        f"  Accuracy:            {accuracy:.4f}\n"
        f"  Pearson correlation: {correlationP:.4f}\n"
        f"  Spearman correlation:{correlationS:.4f}\n"
        f"{'=' * 40}\n"
        f"Elapsed time:          {elapsed:.2f} seconds\n"
    )
    path = os.path.join(results_dir, 'eval_summary.txt')
    with open(path, 'w') as f:
        f.write(summary)
    print(summary)

def generate_plots(loader, model_declared, device, results_dir):
    # Evaluate model
    start = time.time()
    idxt, idxp, mt, mp = collect_metrics(loader, model_declared, device)
    accuracy, correlationP, correlationS = compute_metrics(idxt, idxp, mt, mp)
    #-----------------------------------
    # Save values
    #-----------------------------------
    np.savez(f'{results_dir}/eval-results.npz',
        max_idx_true  = idxt,
        max_idx_pred  = idxp,
        values_true   = mt,
        values_pred   = mp,
    )
    #-----------------------------------
    # Maximum values plot
    #-----------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(idxp, idxt, color='steelblue', s=80)
    ax.text(0.95, 0.95, f'Accuracy: {accuracy:.2f}', transform=ax.transAxes,ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.set_xlabel('Predicted max node', fontsize=13)
    ax.set_ylabel('Expected max node', fontsize=13)
    ax.set_title('Predicted maximum node', fontsize=15)
    ax.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(f'{results_dir}/Indexes_plot.png',dpi=150)
    plt.close()
    #-----------------------------------
    # Predicted values plot
    #-----------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(x = mp, y = mt, color='steelblue', s=80)
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
    plt.savefig(f'{results_dir}/Values_plot.png',dpi=150)
    plt.close()
    #-----------------------------------
    # Predicted values plot
    #-----------------------------------
    elapsed = time.time() - start
    save_eval_summary(results_dir, model_dir, data_dir, accuracy, correlationP, correlationS, elapsed, loader)
    
#-------------------------------
# Section: Load weights and validate
#------------------------------
import numpy as np 
from scipy.stats import pearsonr
from scipy.stats import spearmanr

# Enable GPU
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Load model weights and model
model_dir = '/home/mriveraceron/glv-research/tuning_results/91074c4e25b4/Variant_3'
weights_path = f'{model_dir}/model-weights.pth'
model_declared = model()
checkpoint = torch.load(weights_path, weights_only=True)
model_declared.load_state_dict(checkpoint['model_state_dict'])
model_declared = model_declared.to(device)

#-------------------------------
# Section: Load data
#------------------------------
import os
import glob
from torch_geometric.loader import DataLoader

# Experiment data
# data_dir = '/home/mriveraceron/glv-research/data_null/d2f93775a813'
data_dir = '/home/mriveraceron/glv-research/data_tensors/KBoost_v2_filter'
data_paths = glob.glob(f'{data_dir}/ValBatch_*.pt')

validation_data = []
for path in data_paths:
    validation_data.extend(torch.load(path, weights_only=False))

    
loader = DataLoader(validation_data, batch_size=30, shuffle=False)

#-------------------------------
# Section: Generate directory for validation
#-------------------------------
import os
import matplotlib.pyplot as plt

# base_results_dir = f'/home/mriveraceron/glv-research/Validation/{os.path.basename(model_dir)}'
base_results_dir = f'/home/mriveraceron/glv-research/Validation/91074c4e25b4'
os.makedirs(base_results_dir, exist_ok=True)
val_results_dir = os.path.join(base_results_dir, os.path.basename(data_dir))
os.makedirs(val_results_dir, exist_ok=True)

generate_plots(loader, model_declared, device, val_results_dir)

# For running modify:
#   model directory where the weights are located
#   data_dir: where the batches used for evaluation are stored
#   Results directory

