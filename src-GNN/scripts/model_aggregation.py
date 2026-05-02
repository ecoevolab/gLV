"""
    Model aggregation script
    =======================================
    Description:
    This script is geenerated to perform the model aggregation of best SAGE and GraphConv model.
"""

from scipy.stats import pearsonr, spearmanr
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Ridge, LinearRegression
from itertools import combinations
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

#--------------- Load data --------------- 
# Data for training
graph_dir = '/home/mriveraceron/glv-research/best_models/GraphConv'
sage_dir = '/home/mriveraceron/glv-research/best_models/SAGE'

sage_data = np.load(f'{sage_dir}/metric-values.npz')
graph_data = np.load(f'{graph_dir}/metric-values.npz')

#--------------- Select hyperparameters --------------- 
X_meta = np.column_stack([sage_data['values_pred'], graph_data['values_pred']])
# X_train, X_test, y_train, y_test = train_test_split(X_meta, graph_data['values_true'], test_size=0.3)

model = LinearRegression()
expected =  graph_data['values_true']
model.fit(X_meta,expected)
preds_stacked = model.predict(X_meta)  # Regression
preds_mean = X_meta.mean(axis=1)       # Mean

#------------------- Max node prediction --------------
import glob
import torch 
import os 
# Data for training
data_dir = '/home/mriveraceron/glv-research/data_null'

def data_generator(data_dir, split='train'):
    data_list = []
    paths = glob.glob(f'{data_dir}/*_{split}/*.pt')
    if not paths:
        raise FileNotFoundError(f"No .pt files found under {data_dir}/*_{split}/")
    for path in paths:
        data = torch.load(path, weights_only=False)
        data_list.extend(data)
    print(f"Total samples for {split}: {len(data_list)}")
    # Get directories of data for training
    used_dirs = '\n'.join(f'  {p}' for p in set(os.path.dirname(p) for p in paths))
    return data_list, used_dirs


eval_data, eval_dirs = data_generator(data_dir, split='eval')
data = eval_data[0]
true_values, graph_id, test_maxidx = [], [], []
for i, data in enumerate(eval_data):
    y = data.y.cpu().numpy().flatten()
    true_values.append(y)
    test_maxidx.append(np.argmax(y))
    graph_id.append(np.full(data.x.shape[0], i))

true_values = np.concatenate(true_values)
test_maxidx = np.array(test_maxidx)
graph_id    = np.concatenate(graph_id)

# Verify if predicted node with the max value is the same
# This is just a control
np.array_equal(test_maxidx, graph_data['max_idx_true']) 
np.array_equal(true_values, graph_data['values_true'].squeeze())

#--------------- Correlations ---------------
# Model aggregation correlation
def corr(pred, true, graph_ids, true_ids):
    pred = np.asarray(pred).ravel()
    true = np.asarray(true).ravel()
    graph_ids = np.asarray(graph_ids)
    # Grouped argmax using numpy
    unique_ids, starts = np.unique(graph_ids, return_index=True)
    splits = np.split(pred, starts[1:])
    idx_max_local = np.array([np.argmax(s) for s in splits])
    ppv = np.mean(idx_max_local == true_ids)
    return pearsonr(pred, true)[0], spearmanr(pred, true)[0], ppv

models_data = {
    'aggr_regression':  preds_stacked.squeeze(),
    'aggr_mean': preds_mean,
    'graph': graph_data['values_pred'].squeeze(),
    'sage':  sage_data['values_pred'].squeeze(),
}

y = graph_data['values_true']
true_ids = graph_data['max_idx_true']
graph_ids = graph_id
results = {name: corr(pred, y, graph_ids, true_ids) for name, pred in models_data.items()}
print(results)

# Save values into a npz
np.savez('/home/mriveraceron/glv-research/best_models/aggr_preds.npz', **models_data)
# np.load('/home/mriveraceron/glv-research/best_models/aggr_preds.npz')
#--------------- Summary ---------------
summary = f"""
    Model aggregation summary
    =========================
    Aggregation regression model:
    Pearson correlations {results['aggr_regression'][0]:.5f}
    Spearman correlations {results['aggr_regression'][1]:.5f}
    Positive predicted value: {results['aggr_regression'][2]:.5f}
    Model coefficientes SAGE: {model.coef_[0][0].item():.4f}
    Model coefficientes GraphConv: {model.coef_[0][1].item():.4f}
    Model bias: {model.intercept_.item():.4f}
    -------------------------
    Aggregation mean model:
    Pearson correlations {results['aggr_mean'][0]:.5f}
    Spearman correlations {results['aggr_mean'][1]:.5f}
    Positive predicted value: {results['aggr_mean'][2]:.5f}
    -------------------------
    GraphConv Model:
    Pearson correlations {results['graph'][0]:.5f}
    Spearman correlations {results['graph'][1]:.5f}
    Positive predicted value: {results['graph'][2]:.5f}
    -------------------------
    SAGE Model:
    Pearson correlations {results['sage'][0] :.5f}
    Spearman correlations {results['sage'][1] :.5f}
    Positive predicted value: {results['sage'][2]:.5f}
    """
print(summary)

with open(f'/home/mriveraceron/glv-research/best_models/aggregation_summary.txt', 'w') as f:
    f.write(summary)


#-----------------Generate plots---------------
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

keys = list(models_data.keys())
n = len(keys)
box_style = dict(boxstyle='square,pad=0.4', facecolor='#DDEEFF',edgecolor='#AABBCC', linewidth=0.8)

plt.clf()
fig, axes = plt.subplots(n, n, figsize=(12, 10))
for i, key_y in enumerate(keys):
    ax = axes[i, n - 1]
    pos = ax.get_position()  # gets the axes bounding box in figure coords
    # Add right titles
    fig.text(
        pos.x1 + 0.05,          # just to the right of the last column
        pos.y0 + pos.height / 2, # vertically centered on the row
        key_y,
        fontsize=10, fontweight='bold',
        ha='left', va='center',
        rotation=270,
        bbox=box_style,
        transform=fig.transFigure
    )
    for j, key_x in enumerate(keys):
        ax = axes[i, j]
        if i >= j:
            # Lower triangle: scatter + correlation annotation
            x = models_data[key_x].numpy() if hasattr(models_data[key_x], 'numpy') else np.array(models_data[key_x])
            y = models_data[key_y].numpy() if hasattr(models_data[key_y], 'numpy') else np.array(models_data[key_y])
            # plotting method
            # x_log = np.log2(x + 1)
            # y_log = np.log2(y + 1)
            hb = ax.hexbin(x, y, gridsize=40, cmap='YlOrRd', mincnt=1, bins='log', linewidths=0.2)
            cb = plt.colorbar(hb, ax=ax, pad=0.15, shrink=0.85)
            # cb.set_label('Count (log scale)', fontsize=10, labelpad=8)
            cb.ax.tick_params(labelsize=9)
            # perfect correlation line
            ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.8, zorder=5)
            # Correlation annotation top-right corner of the scatter
            correlationP, _ = stats.pearsonr(x, y)
            correlationS, _ = stats.spearmanr(x, y)
            ax.text(0.97, 0.97,
                f'$r$ = {correlationP:.3f}\n$\\rho$ = {correlationS:.3f}',
                transform=ax.transAxes, ha='right', va='top', fontsize=8,
                bbox=dict(boxstyle='round,pad=0.4,rounding_size=0.5', facecolor='white',
                            edgecolor='gray', linewidth=0.8, alpha=0.85),
                zorder=10
            )
            ax.tick_params(labelsize=7)
        else:
            # Hide visual elements but keep axes visible so title shows
            ax.set_facecolor('none')
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
        # title
        if i == 0:
            ax.text(0.5, 1.05, key_x, fontsize=10, fontweight='bold',transform=ax.transAxes, ha='center', va='bottom', bbox=box_style)
        # Add Y-axis label to only the first column
        if j == 0:
            ax.set_ylabel(r'Predicted value $\log(1+p)$', fontsize=9, labelpad=4)
        else:
            ax.set_ylabel('')
        # Add X-axis label to only thee last row
        if i == n - 1: 
            ax.set_xlabel(r'Predicted value $\log(1+p)$', fontsize=9, labelpad=4)
        else:
            ax.set_xlabel('')
        # Add gradient label
        if i == j:
            cb.set_label('Count (log scale)', fontsize=10, labelpad=8)
        else:
            cb.ax.tick_params(labelsize=9)
        

            
plt.suptitle('Pairwise Model Keystoneness Predictions', fontsize=14, fontweight='bold', y=0.95)
plt.savefig(f'/home/mriveraceron/glv-research/best_models/correlations_plots_V2.png', dpi=300, bbox_inches='tight')
# plt.show()
# plt.close('all')



#-----------------Generate plots---------------

keys = list(models_data.keys())
n = len(keys)
box_style = dict(boxstyle='round,pad=0.4,rounding_size=0.5', facecolor='#DDEEFF',edgecolor='#AABBCC', linewidth=0.8)

plt.clf()
fig, axes = plt.subplots(n, 1, figsize=(12, 10), sharex=True)
for i, key_x in enumerate(keys):
    ax = axes[i]
    # Predicted and expected values
    y = models_data[key_x].numpy() if hasattr(models_data[key_x], 'numpy') else np.array(models_data[key_x])
    x = expected.squeeze()
    # plotting method
    x_log = np.log2(x + 1)
    y_log = np.log2(y + 1)
    # plotting method
    hb = ax.hexbin(x, y, gridsize=40, cmap='YlOrRd', mincnt=1, bins='log', linewidths=0.2)
    # Color gradient bar
    cb = plt.colorbar(hb, ax=ax, pad=0.15, shrink=0.85)
    cb.set_label('Count (log scale)', fontsize=10, labelpad=2)
    cb.ax.tick_params(labelsize=9)
    # perfect correlation line
    ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.8, zorder=5)
    # Correlation annotation top-right corner of the scatter
    correlationP, _ = stats.pearsonr(x, y)
    correlationS, _ = stats.spearmanr(x, y)
    ax.text(0.97, 0.97,
        f'$r$ = {correlationP:.3f}\n$\\rho$ = {correlationS:.3f}',
        transform=ax.transAxes, ha='right', va='top', fontsize=8,
        bbox=dict(boxstyle='round,pad=0.4,rounding_size=0.5', facecolor='white',
                    edgecolor='gray', linewidth=0.8, alpha=0.85),zorder=10
    )
    ax.tick_params(labelsize=7)
    # title
    ax.text(0.5, 1.05, key_x, fontsize=10, fontweight='bold',transform=ax.transAxes, ha='center', va='bottom', bbox=box_style)
    # Add Y-axis label 
    ax.set_ylabel(r'Predicted value $\log(1+p)$', fontsize=9, labelpad=4)
    #ax.set_xlabel(r'Expected value $\log(1+p)$', fontsize=9, labelpad=15)
        
axes[-1].set_xlabel("Expected value log(1 + p)", labelpad=10)
fig.subplots_adjust(bottom=0.08)  # increase to push everything up          
plt.suptitle('Expected vs. predicted Keystoneness', fontsize=14, fontweight='bold', y=0.95)
plt.savefig(f'/home/mriveraceron/glv-research/best_models/tmp_plots.png', dpi=300, bbox_inches='tight')
# plt.show()
# plt.close('all')

# Data analysis
for mod in models_data:
    arr = models_data[mod]
    print(f"Min: {arr.min()}, Max: {arr.max()}, Mean: {arr.mean():.2f} for {mod}")
    print(f'number of zeros: {np.sum(arr == 0)} for {mod}')