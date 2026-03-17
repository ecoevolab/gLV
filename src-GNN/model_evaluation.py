"""
GNN model evaluation
=================================================
Purpose:
    Evaluate GNN performance using pre-trained model weights.

Dependencies:
    torch==2.8, pandas==2.3.3, numpy==2.0.2

Author: Manuel Rivera
Date:   March 16, 2026
"""

#-------------------------------
# Model definition
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
# Load model weights (.pth)
#-------------------------------
import os

weights_dir = '/home/mriveraceron/glv-research/model_weights'
name = 'Boosted_filtered_V3.pth'
model_weights = os.path.join(weights_dir, name)
# Load weights
model_declared = model()
checkpoint = torch.load(model_weights, weights_only=True)
model_declared.load_state_dict(checkpoint['model_state_dict'])

#-------------------------------
# Evaluate model
#-------------------------------
import glob

# Load data
tensors_dir = '/home/mriveraceron/glv-research/data_tensors/'
experiment_name = 'KBoost_v2_filter'
experiment_data = os.path.join(tensors_dir,experiment_name)
val_paths = glob.glob(f"{experiment_data}/ValBatch_*.pt")
result_path = os.path.join('/home/mriveraceron/glv-research/Validation',experiment_name)

model_declared.eval()
# Run inference without computing gradients
predicted_vector, expected_vector = [], []
idx_max_expected, idx_max_predicted= [], []     

with torch.no_grad():
    for path in val_paths:
        data_list = torch.load(path, weights_only=False)
        for data in data_list:
            # Generate prediction
            y_pred = model_declared(data)
            # Append values
            expected_vector.append(data.y)
            predicted_vector.append(y_pred)
            # Append maximum node
            idx_max_expected.append(torch.argmax(data.y, dim=0))
            idx_max_predicted.append(torch.argmax(y_pred, dim=0))

# Convert once at the end
expected_vector  = torch.cat(expected_vector).cpu().numpy()
predicted_vector = torch.cat(predicted_vector).cpu().numpy()
idx_expected_vector = torch.cat(idx_max_expected).cpu().numpy()
idx_predicted_vector = torch.cat(idx_max_predicted).cpu().numpy()

#-------------------------------
# Save performance
#-------------------------------
import numpy as np

os.makedirs(result_path, exist_ok=True)
print('The validation directory will be:', result_path, '\n')
#-------------------------------
# Save max indexes
np.save(f'{result_path}/max_idx_true.npy', idx_expected_vector)
np.save(f'{result_path}/max_idx_pred.npy', idx_predicted_vector)
#-------------------------------
# Save metrics
np.save(f'{result_path}/values_true.npy', expected_vector)
np.save(f'{result_path}/values_pred.npy', predicted_vector)

#-------------------------------
# Generate correlation plot
#-------------------------------
import matplotlib.pyplot as plt

# Calculate accuracy
accuracy = np.mean(np.array(idx_expected_vector) == np.array(idx_predicted_vector))

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.scatter(idx_expected_vector, idx_predicted_vector, color='steelblue', s=80)
ax.text(0.95, 0.95, f'Accuracy: {accuracy:.2f}', transform=ax.transAxes,ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax.set_xlabel('Predicted max node', fontsize=13)
ax.set_ylabel('Expected max node', fontsize=13)
ax.set_title('Predicted maximum node', fontsize=15)
ax.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig(f'{result_path}/Indexes_plot.png',dpi=150)

#-------------------------------
# Generate values plot
#-------------------------------
from scipy.stats import pearsonr
from scipy.stats import spearmanr

# Generate correlations
correlationP, _ = pearsonr(predicted_vector, expected_vector)
correlationS, _ = spearmanr(predicted_vector, expected_vector)

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.scatter(x = expected_vector, y = predicted_vector, color='steelblue', s=80)
ax.text(0.95, 0.95, 
    f'Pearson Correlation: {correlationP.item():.4f}\nSpearman Correlation: {correlationS.item():.4f}',
    transform=ax.transAxes, ha='right', va='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
)
ax.set_xlabel('Predicted values', fontsize=13)
ax.set_ylabel('Expected values', fontsize=13)
ax.set_title('Keystoneness values', fontsize=15)
ax.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig(f'{result_path}/Values_plot.png',dpi=150)

#-------------------------------
# Generate Summary
#-------------------------------
summary = f"""
Model Evaluation Summary
=========================
Model: {model_weights}
Validation samples:   {len(predicted_vector)}
Validation data: {experiment_data}
Validation results: {result_path}

Pearson Correlation:  {correlationP.item():.4f}  
Spearman Correlation: {correlationS.item():.4f}  
Maximum node accuracy: {accuracy.item():.4f}  
"""

with open(f'{result_path}/evaluation_summary.txt', 'w') as f:
    f.write(summary)

print(summary)