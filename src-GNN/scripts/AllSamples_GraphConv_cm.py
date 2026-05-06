"""
Generate confussion matrix for AllSamples GraphConv model
--------------------------
Previously best GraphConv model was trained using all available samples.
This script purpose is to generate the confussion matrix of the model
"""

# Import libraries
import glob
import torch
import os
import numpy as np
import pandas as pd

# We will regenerate the evaluation data for counting the nodes number of each graph
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

# Generate validation data
eval_data, eval_dirs = data_generator(data_dir, split='eval')


# Generate number of nodes
graph_nodes = []
for d in eval_data:
    graph_nodes.append(d.x.shape[0])

graph_nodes = np.array(graph_nodes)
print(f'Prediction expected by chance: {1 / graph_nodes.mean()}')


# Load predictions
data = np.load('/home/mriveraceron/glv-research/best_models/GraphConv/metric-values.npz')
idxt        = data['max_idx_true']
idxp        = data['max_idx_pred']
print(f'Lengths should match idxt {len(idxt)}| idxp {len(idxp)}| nodes {len(graph_nodes)}')

# Confusion matrix
tp = sum((idxt == idxp).astype(int)) # Predicted correctly
fp = sum((idxt != idxp).astype(int)) # Failed to predict
fn = fp
tn = sum((graph_nodes - 1)) - fp    # Predicted correctly (without top node)
acc = (tp + tn) / (tp + tn + fp + fn)
ppv = tp / (tp + fp)
print(f'Model had an accuracy of {acc} and precision of {ppv}')

# Convert to df
df_cm = pd.DataFrame(
    [[tp, fp], [fn, tn]],
    index   = ['Expected_Positive',  'Expected_Negative'],
    columns = ['Predicted_Positive','Predicted_Negative']
)
# Save it 
df_cm.to_csv(f'/home/mriveraceron/glv-research/best_models/GraphConv/confussion_matrix.csv')
   