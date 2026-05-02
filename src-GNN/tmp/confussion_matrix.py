"""
    Generating old correlation plots and confussion matrices.
"""



import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix

# Load data
DIR = '/home/mriveraceron/glv-research/Results/Boosted_keystone/Filtered_AllFeats_V3'
mt = np.load(f'{DIR}/values_true.npy')
mp = np.load(f'{DIR}/values_pred.npy')
idx_true = np.load(f'{DIR}/max_idx_true.npy')
idx_pred = np.load(f'{DIR}/max_idx_pred.npy')

# Log-transform values 
x = np.log1p(mp)
y = np.log1p(mt)

# Compute correlations
r_p, _ = pearsonr(mt, mp)
r_s, _ = spearmanr(mt, mp)

# ----------- Values plot -----------
fig, ax = plt.subplots(figsize=(4.5, 4.5))
bbox=dict(boxstyle='round,pad=0.4,rounding_size=0.5', facecolor='white', edgecolor='gray', linewidth=0.8, alpha=0.85)
hb = ax.hexbin(x, y, gridsize=40, cmap='YlOrRd', mincnt=1, bins='log', linewidths=0.2)
cb = plt.colorbar(hb, ax=ax, pad=0.02, shrink=0.85)
cb.set_label('Count (log scale)', fontsize=10)
cb.ax.tick_params(labelsize=9)
ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.8, zorder=5)
ax.set_xlabel('Predicted keystoneness', labelpad=8)
ax.set_ylabel('Expected keystoneness', labelpad=8)
ax.set_title('Boosted data keystoneness values', pad=10)
ax.text(
    0.5, 1.01,                  # x=center, y=just below title
    'GraphConv model',
    transform=ax.transAxes,
    ha='center', va='bottom',
    fontsize=10, color='gray'
)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.grid(False)
ax.text(0.05, 0.95,
    f'$r$ = {r_p:.3f}  (Pearson)\n$\\rho$ = {r_s:.3f}  (Spearman)',
    transform=ax.transAxes, ha='left', va='top', fontsize=9,
    bbox=dict(boxstyle='round,pad=0.4,rounding_size=0.5',facecolor='white', edgecolor='gray', linewidth=0.8, alpha=0.85),
    zorder=10
)
plt.tight_layout()
# plt.savefig(f'{DIR}/Hexbin_plot.pdf', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(f'{DIR}/Hexbin_plot.png', dpi=300, bbox_inches='tight')
plt.show()

# ----------- Generate confusion matrix -----------
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import os
import torch

# Load number of nodes per graphs
data_dir = '/home/mriveraceron/glv-research/data_tensors/Boosted_filtered'
list_files = sorted([f for f in os.listdir(data_dir) if f.endswith('.pt')])
file = list_files[0]
num_nodes = 0
for file in list_files:
    data_list = torch.load(f'{data_dir}/{file}')
    for data in data_list:
        num_nodes += data.x.shape[0]

# Create confussion matrix
cm = confusion_matrix(idx_true, idx_pred)
cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)  # row-normalized

fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(
    cm_norm,
    annot=True,
    fmt='d',
    cmap='YlOrRd',
    ax=ax,
    linewidths=0.5,
    linecolor='gray'
)

ax.set_title('Confusion Matrix', fontsize=15, pad=12, fontweight = 'bold')
ax.set_xlabel('Predicted Label', fontsize=12, fontweight = 'bold')
ax.set_ylabel('True Label', fontsize=12, fontweight = 'bold')
# --- Accuracy box ---
accuracy = cm.diagonal().sum() / cm.sum()
ax.text(
    0.98, 0.98,
    f'Accuracy: {accuracy:.2%}',
    transform=ax.transAxes,
    ha='right', va='top',
    fontsize=11, fontweight='bold', color='#222222',
    bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#aaaaaa', alpha=0.85)
)
plt.tight_layout()
plt.savefig(f'{DIR}/confusion_matrix.jpg', dpi=150, bbox_inches='tight')
plt.show()