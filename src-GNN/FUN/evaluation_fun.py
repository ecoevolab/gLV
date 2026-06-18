"""
Evaluation function on trined GNN model 
=================================================
Purpose:
    Declare evaluation function for gnn evaluation. Evaluation data must be in the form of a DataLoader for faster evaluation time.

Input Data:
    - model (torch.nn.Module): Previously declared model architecture (e.g., GCN, GraphConv, SAGEConv).
    - device (str | torch.device): Device to evaluate the model on (e.g., 'cuda' or 'cpu').
    - eval_data (List[torch_geometric.data.Data]): Evaluation dataset, expected to be a list of Data() objects 
        with x being the node features and y being the target variable.


Dependencies:
    torch==2.4.0+cu121, pandas==2.0.1, tqdm==4.66.5

Returns:
    - metrics (namedtuple): Named tuple containing:
        - idxt (np.ndarray): True node indices of the maximum target value for each graph.
        - idxp (np.ndarray): Predicted node indices of the maximum target value for each graph.
        - mt (np.ndarray): True target values for all nodes across all graphs.
        - mp (np.ndarray): Predicted target values for all nodes across all graphs.
        - nodes (np.ndarray): Number of nodes in each graph.
        
    - performance (namedtuple): Named tuple containing:
        - ppv (float): Positive Predictive Value.
        - corrP (float): Pearson's correlation coefficient between true and predicted target values.
        - corrS (float): Spearman's rank correlation coefficient between true and predicted target values

Author: Manuel Rivera
Date:   June 18, 2026
"""

# Import libraries
import os
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"  # Must be before torch import
import torch
import numpy as np
from collections import namedtuple
from scipy.stats import pearsonr, spearmanr
import logging
from torch_geometric.utils import unbatch
from torch_geometric.loader import DataLoader

# Declare named tuples for metrics and performance results
MetricsResult     = namedtuple('MetricsResult', ['idxt', 'idxp', 'mt', 'mp', 'nodes'])
PerformanceResult = namedtuple('PerformanceResult', ['ppv', 'corrP', 'corrS'])

# Function to generate predictions on the evaluation set.
def collect_metrics(eval_data, model, device):
    # Declare lists to collect metrics
    idxt, idxp, mt, mp, nodes = [], [], [], [], []
    try:
        model.eval()  # Set model to evaluation mode
        with torch.no_grad():  # diable gradient 
            for batch in eval_data:
                batch = batch.to(device)      # move data to device
                out = model(batch)   # forward pass
                for y, o in zip(unbatch(batch.y, batch.batch), unbatch(out, batch.batch)):
                    idxt.append(torch.argmax(y, dim=0))   # TRUE node-index for maximum target value (0-dim tensor)
                    idxp.append(torch.argmax(o, dim=0))   # PREDICTED node-index for maximum target value (0-dim tensor)
                    # Append true and predicted target values for correlation computation
                    mt.append(y)   
                    mp.append(o)       
                    # Number of nodes in each graph (for later analysis of performance vs. graph size)
                    nodes.append(y.shape[0])              
    finally:
        model.train()
    return MetricsResult(
            torch.cat(idxt).cpu().numpy(),   
            torch.cat(idxp).cpu().numpy(),   
            torch.cat(mt).cpu().numpy(),
            torch.cat(mp).cpu().numpy(),
            np.array(nodes),                   
        )

# Function to compute performance metrics (PPV, Pearson's r, Spearman's rho) based on collected metrics.
def compute_metrics(metrics_list):
    idxt, idxp = metrics_list.idxt, metrics_list.idxp
    mt, mp = metrics_list.mt, metrics_list.mp
    ppv = np.mean(np.array(idxt) == np.array(idxp))
    # Check for constant arrays before computing correlations to avoid errors
    if np.std(mt) == 0 or np.std(mp) == 0:
        log.warning("Cannot compute correlation: one input is constant.")
        correlationP = correlationS = float('nan')
    else:
        correlationP, _ = pearsonr(mt.flatten(), mp.flatten())
        correlationS, _ = spearmanr(mt.flatten(), mp.flatten())
    return PerformanceResult(ppv, correlationP, correlationS)

# Wrapper function to evaluate the model with evaluation data and compute performance metrics.
def evaluate_split(eval_data, model, device, batch_size=30):
    eval_loader = DataLoader(eval_data, batch_size=batch_size, shuffle=False)
    metrics     = collect_metrics(eval_loader, model, device)
    performance = compute_metrics(metrics)
    return metrics, performance