"""
Continue gnn model training from pre-trained weights
=================================================
Purpose:
    Continue trainining of Graph Neural Network (GNN) using pre-trained weights.

Input Data:
    - X          : Node feature matrix, where each row represents the statistics of a node in the network.
    - Y          : Target variable representing keystoneness.
    - edge_index : Graph connectivity in who->whom format, defining species interactions.
    - edge_attr  : Edge weights representing the strength of each interaction.

Dependencies:
    torch==2.8, pandas==2.3.3, numpy==2.0.2

Author: Manuel Rivera
Date:   June 18, 2026
"""

# Declare libraries 
import os
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"  # Must be before torch import
import torch 
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv, SAGEConv
import logging
import importlib.metadata as metadata
import sys
from collections import namedtuple
import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import torch.nn.functional as F
import torch.optim as optim
from torch_geometric.loader import DataLoader
from torch_geometric.utils import unbatch
from tqdm import tqdm
import time
from scipy.stats import pearsonr, spearmanr
import random
import glob
from functools import partial
import textwrap
import shutil
import pickle


# Section: Declare logging configuration and handler to redirect the logs to a file.

# Create result directory
results_dir = '/home/mriveraceron/glv-research/tuning_results/testing_funs'
shutil.rmtree(results_dir) if os.path.exists(results_dir) else None
os.makedirs(results_dir, exist_ok=True)


def make_logger(name, filepath):
    # Generate a log file to keep track of training progress and results. 
    # Handlers determine where the logs are sent (e.g., file, console), and formatters specify the layout of log messages.
    # Propagate: when a logger handles a message, it also passes it up to its parent logger, resulting
    # in the message appearing twice.
    logger = logging.getLogger(name)
    # info, warning, and error messages, but ignores debug.
    logger.setLevel(logging.INFO)   
    logger.handlers.clear()  # removes existing handlers 
    logger.propagate = False 
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

def reroute_logger(source_name: str, target_logger: logging.Logger) -> None:
    # Redirect functions logger's output into another logger's handlers.
    source_logger = logging.getLogger(source_name)
    source_logger.setLevel(logging.INFO)
    # remove existing handlers and prevent propagation to avoid duplicate logs
    source_logger.handlers.clear() 
    source_logger.propagate = False
    for handler in target_logger.handlers:
        source_logger.addHandler(handler)
 
# Two independent loggers
log     = make_logger('run_log',  f'{results_dir}/run_log.txt')

# Reroute evaluation logger to the main logger to keep all logs in a single file.
reroute_logger('FUN', log)

#-------------------------------
# Section: Declare model
# Note: For continuing training, the model architecture must match the one used for pre-training.

# This researche employed GraphConv and SAGEConv, but other models can be used as well.
# Both models will be declares here, but only one will be used for training.
class GraphConv_model(nn.Module):
    def __init__(self, in_channels, hidden_channels=64, num_layers=5):
        super().__init__()
        dims = [in_channels] + [hidden_channels] * (num_layers - 1) + [1]
        self.convs = nn.ModuleList(GraphConv(dims[i], dims[i+1]) for i in range(num_layers))
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_weights
        for conv in self.convs[:-1]:
            x = F.relu(conv(x, edge_index, edge_weight))
        return torch.sigmoid(self.convs[-1](x, edge_index, edge_weight))
    
class SAGE_model(nn.Module):
    def __init__(self, in_channels, hidden_channels=64, num_layers=5):
        super().__init__()
        dims = [in_channels] + [hidden_channels] * (num_layers - 1) + [1]
        self.convs = nn.ModuleList(SAGEConv(dims[i], dims[i+1]) for i in range(num_layers))
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_weights
        for conv in self.convs[:-1]:
            x = F.relu(conv(x, edge_index, edge_weight))
        return torch.sigmoid(self.convs[-1](x, edge_index, edge_weight))

#-------------------------
# Section: Generate data for continuin training
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
    return data_list

train_data = data_generator(data_dir, split='train')
eval_data = data_generator(data_dir, split='eval')

print(f"Loaded data from: {data_dir} | Training samples: {len(train_data)} | Evaluation samples: {len(eval_data)}")


# Section: Seed function for reproducibility
def seed_fn(seed=42):
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)        
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True, warn_only=True) 

# Section: Summary function to log training results
def summarize(model_name, model, tr_size, ev_size, perf_tr, saving_path, elapsed):
    summary = textwrap.dedent(f"""
        \n-----------------------------------------------
        Model name:          {model_name}
        Model declared:\n      {model}
        Training samples:    {tr_size}
        Validation samples:  {ev_size} \n
        Pearson Correlation:   {perf_tr.corrP}
        Spearman Correlation:  {perf_tr.corrS}
        Maximum node PPV:      {perf_tr.ppv} \n
        Running time:   {elapsed:.2f}s ({elapsed/60:.2f} min)
        Results saved:  {saving_path}\n
    """).strip()
    log.info(summary)
    return summary


# Section: Source training and evaluation functions.


# Section: Train model and evaluate it
