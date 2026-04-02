# Script to generate evaluation of model

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
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import spearmanr
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

#-------------------------------
# Section: Load variants table
#-------------------------------
import pandas as pd
import numpy as np

# Load variants table
exp_dir = '/home/mriveraceron/glv-research/tuning_results/91074c4e25b4'
tuning_df = pd.read_csv(f'{exp_dir}/tuning_results.csv')
tuning_df.iloc[:, 5:8] = np.nan  

# Load model, device 
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model_declared = model()

#-------------------------------
# Section: Load data
#-------------------------------
import os
import glob
from torch_geometric.loader import DataLoader

data_paths = glob.glob(f'/home/mriveraceron/glv-research/data_tensors/91074c4e25b4/ValBatch_*.pt')
validation_data = []
for path in data_paths:
    validation_data.extend(torch.load(path, weights_only=False))

    
loader = DataLoader(validation_data, batch_size=30, shuffle=False)

#-------------------------------
# Section: Generate results
#-------------------------------
import os
val_results = '/home/mriveraceron/glv-research/Validation/91074c4e25b4'
os.makedirs(val_results, exist_ok=True)

for i, row in tuning_df.iterrows():
    #--------------------
    # Generate model
    model_declared = model(hidden_channels=row['channels'], num_layers=row['layers'])
    #--------------------
    # Load model variant
    variant = row['model_id']
    weights_path = f'{exp_dir}/{variant}/model-weights.pth'
    print(f'>> The weights path are: {weights_path}')
    # Load model weights
    checkpoint = torch.load(weights_path, weights_only=True)
    model_declared.load_state_dict(checkpoint['model_state_dict'])
    model_declared = model_declared.to(device)
    #--------------------
    # Evaluate model
    idxt, idxp, mt, mp = collect_metrics(loader, model_declared, device)
    # Get metrics
    accuracy, corrP, corrS = compute_metrics(idxt, idxp, mt, mp)
    # Assgin values
    tuning_df.loc[i, 'accuracy_idx'] = accuracy
    tuning_df.loc[i, 'pearson_corr'] = corrP
    tuning_df.loc[i, 'spearman_corr'] = corrS

# Save table
tmp = tuning_df.sort_values("pearson_corr", ascending=False)
tmp.to_csv(f'{val_results}/validation_results.csv', index=False)  # save after each run
print(tmp.to_markdown(index=False))