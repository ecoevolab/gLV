"""
Generate new single model variant and append it to the list of variants.

Previously hyperparameters were already tested and the best ones were selected.
After determining the hyperparameters, sample size was tested and a tendency was observed
where more samples increased the performance of the model.

Therefore, new variant will include all available traning samples.

14-March-2026
"""


# Imports
import pandas as pd
import logging
import sys
import os 
import glob

# For package versions
from importlib.metadata import packages_distributions, version
import importlib.metadata as metadata

# For model
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv

# For metrics
from collections import namedtuple
from scipy.stats import pearsonr
from scipy.stats import spearmanr

import torch.optim as optim
import numpy as np
import random

# For training loop
import time
from tqdm import tqdm


# For memory usage
import glob
import pickle
import pkg_resources

# Import functions from FUN.py
import sys      
sys.path.append("/home/mriveraceron/glv-research/gLV/src-GNN/")     # Specify functions directory path
from FUN import collect_metrics, compute_metrics, training_fn, summarize, seed_fn
#-----------------------
# Set up logging
#-----------------------
results_dir = '/home/mriveraceron/glv-research/tuning_results/GraphConv_sample_size/Variant_9'
os.makedirs(results_dir, exist_ok=True)

def make_logger(name, filepath):
    logger = logging.getLogger(name)
    # info, warning, and error messages, but ignores debug.
    logger.setLevel(logging.INFO)   
    # Handler: Where does a log goes to? (file, console, etc.)
    logger.handlers.clear()  # removes existing handlers 
    # Propagate: when a logger handles a message, it also passes it up to its parent logger.
    # With this messages do not appear twice
    logger.propagate = False  # don't bubble up to root logger
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    # File handler
    file_handler = logging.FileHandler(filepath, mode='w')
    file_handler.setFormatter(formatter)
    # Console handler
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)  # Add timestamps
    logger.addHandler(file_handler)         # Log to file
    logger.addHandler(stream_handler)       # Log to console
    return logger


# Two independent loggers
log     = make_logger('run_log',  f'{results_dir}/run_log.txt')
pkglog  = make_logger('pkg_log',  f'{results_dir}/pkgs_log.txt')

#  Route FUN.py logs into run_log 
fun_logger = logging.getLogger("FUN")
fun_logger.setLevel(logging.INFO)
fun_logger.handlers.clear()
fun_logger.propagate = False
for handler in log.handlers:       # reuse run_log handlers
    fun_logger.addHandler(handler)

#-----------------------
# Section: Print imported packages and versions
#-----------------------
# Print imported packages and versions 
installed = {dist.metadata['Name']: dist.metadata['Version'] for dist in metadata.distributions()}

for package, version in sorted(installed.items()):
    pkglog.info(f"{package}=={version}")
#-----------------------
# Load the directory where model was stored
#-----------------------
model_dir = "/home/mriveraceron/glv-research/tuning_results/GraphConv_sample_size"
variant_table = pd.read_csv(f'{model_dir}/tuning_results.csv')

# Generate new row for new variant
new_row = pd.DataFrame([{'model_id': 'Variant_9', 'channels':64, 'layers':5, 'learning_rate':0.001, 'epochs':700}])

# Append new row to the existing table
new_df = pd.concat([variant_table, new_row], ignore_index=True)

#-----------------------
# Load data
#-----------------------
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
    log.info(f"Total samples for {split}: {len(data_list)}")
    # Get directories of data for training
    used_dirs = '\n'.join(f'  {p}' for p in set(os.path.dirname(p) for p in paths))
    return data_list, used_dirs


train_data, train_dirs = data_generator(data_dir, split='train')
eval_data, eval_dirs = data_generator(data_dir, split='eval')

#-------------------------------
# Section: Declare model
#-------------------------------
class GraphConvModel(nn.Module):
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
# Section: Run model
#-------------------------------

# Constant model parameters
device  = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
loss_fn = nn.MSELoss()

# Declare model hyperparameters
size = len(train_data)  # Use all available training samples for this variant
lr = new_row['learning_rate'].iloc[0]
model_name = new_row['model_id'].iloc[0]
epochs = new_row['epochs'].iloc[0]
channels = int(new_row['channels'].iloc[0])
layers = int(new_row['layers'].iloc[0])

# Early stopping parameters
eval_interval = 50
patience = 2
batch_size = 30
log.info(f'>> Starting model {model_name} with train size {size} and learning rate {lr}')
#------------------------
# Section: Run model
#------------------------
n_seed = 42
seed_fn(seed=n_seed)
# Model parameters
model_declared = GraphConvModel(hidden_channels=channels,num_layers=layers).to(device)
optimizer = optim.Adam(model_declared.parameters(), lr=lr)
log.info(f'>> Variant results will be saved at: {results_dir}')
# Slice data for training
random.shuffle(train_data)
# Run training function
weights_path = f'{results_dir}/model_weights.pth'
try:
    loss_history, metrics_list, performance_list, total_elapsed = training_fn(model_declared, device, train_data, eval_data, weights_path, loss_fn, optimizer, epochs, eval_interval, patience, batch_size)
except ValueError as e:
    log.info(f"Error occurred: {e}")
#------------------------
# Add results to dataframe
#------------------------
total_bytes = len(pickle.dumps(train_data))
row = new_df.loc[new_df['model_id'] == 'Variant_9']
row.loc[:, 'train_size'] = len(train_data)
row.loc[:, 'accuracy_idx'] = performance_list.acc
row.loc[:, 'pearson_corr'] = performance_list.corrP
row.loc[:, 'spearman_corr'] = performance_list.corrS
row.loc[:, 'mem_usage(mb)'] = total_bytes / 1e6
row.loc[:, 'eval_size'] = len(eval_data)
# Add row to dataframe
new_df.loc[new_df['model_id'] == 'Variant_9'] = row
# Save table every row
# table_dir = '/home/mriveraceron/glv-research/tuning_results/GraphConv_sample_size'
new_df.to_csv(f'{results_dir}/tuning_results_V2.csv', index=False)

#------------------------
# Section: Generate summary
#------------------------
# Define namedtuple for results
ExtraInfo = namedtuple('ExtraInfo', ['epochs_runned', 'n_seed', 'total_elapsed', 'validation_samples', 'eval_interval', 'patience', 'batch_size'])

# ['epochs_runned', 'n_seed', 'total_elapsed', 'validation_samples', 'eval_interval', 'patience', 'batch_size']
extra_info = ExtraInfo(len(loss_history), n_seed, total_elapsed, len(eval_data), eval_interval, patience, batch_size)
summarize(model_declared, optimizer, row, train_dirs, eval_dirs, performance_list, results_dir, extra_info)

#------------------------
# Section: Save metrics result_exp_dir
#------------------------
np.savez(f'{results_dir}/metric-values.npz',
    max_idx_true  = metrics_list.idxt,
    max_idx_pred  = metrics_list.idxp,
    values_true   = metrics_list.mt,
    values_pred   = metrics_list.mp,
    loss_history  = loss_history
)
