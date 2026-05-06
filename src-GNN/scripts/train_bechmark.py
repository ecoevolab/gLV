"""
DataLoader vs. Standard Training Method
========================================
This script benchmarks and compares the training time of a model
using identical architecture but different data loading strategies.
"""


# Imports
import os
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"  # Must be before torch import
import logging
import importlib.metadata as metadata
import sys
from collections import namedtuple
import torch
from torch_geometric.nn import GraphConv
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

#-------------- Set up logging --------------
results_dir = '/home/mriveraceron/glv-research/Results/train_benchmark'
shutil.rmtree(results_dir) if os.path.exists(results_dir) else None
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


# Print imported packages and versions 
installed = {dist.metadata['Name']: dist.metadata['Version'] for dist in metadata.distributions()}

for package, version in sorted(installed.items()):
    pkglog.info(f"{package}=={version}")

#--------------- Running ----------------------

def training_fn(model, device, data_train, weights_path, loss_fn, optimizer, epochs, shuffle = True, method = 'DL', batch_size=30):
    model.train()
    if (method == 'DL'):
        log.info(f'Training with DataLoader method \n')
        data = DataLoader(data_train, batch_size, shuffle = shuffle)
    else:
        log.info(f'Training with standard (1x1) method \n')
        data = data_train
    loss_history, elapsed_history = [], []
    total_elapsed = 0
    for epoch in tqdm(range(epochs), desc="Training"):
        start = time.time()
        epoch_loss = 0
        for d in data:
            d = d.to(device)
            optimizer.zero_grad()
            out = model(d)
            loss = loss_fn(out, d.y)
            if torch.isnan(loss):
                raise ValueError(f"NaN loss detected at epoch {epoch}")
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            epoch_loss += loss.item()
        epoch_elapsed = time.time() - start		# Epoch elapsed time
        total_elapsed += epoch_elapsed			# Add it to model elapsed time
        loss_history.append(epoch_loss)			# Append loss
        elapsed_history.append(epoch_elapsed)	# Add epoch elpsed time
        if epoch % 50 == 0:
            log.info(f"Epoch {epoch}: Loss={epoch_loss:.6f} | Time={epoch_elapsed:.2f}s")
    # Save after all epochs
    torch.save({
        'epoch': epochs - 1,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'loss': epoch_loss
    }, weights_path)
    log.info(f"Total elapsed: {total_elapsed:.2f}s ({total_elapsed/60:.2f} min)")
    return np.array(loss_history), total_elapsed, elapsed_history


# Namedtuple
MetricsResult     = namedtuple('MetricsResult', ['idxt', 'idxp', 'mt', 'mp', 'nodes'])
PerformanceResult = namedtuple('PerformanceResult', ['ppv', 'corrP', 'corrS'])

def collect_metrics(loader, model_declared, device):
    idxt, idxp, mt, mp, nodes = [], [], [], [], []
    try:
        model_declared.eval()
        with torch.no_grad():
            for batch in loader:
                batch = batch.to(device)
                out = model_declared(batch)
                y_list  = unbatch(batch.y, batch.batch)
                out_list = unbatch(out, batch.batch)
                for y, o in zip(y_list, out_list):
                    idxt.append(torch.argmax(y, dim=0))   # 0-dim tensor
                    idxp.append(torch.argmax(o, dim=0))   # 0-dim tensor
                    mt.append(y)
                    mp.append(o)
                    nodes.append(y.shape[0])              # plain int
        return MetricsResult(
            torch.cat(idxt).cpu().numpy(),   # stack, not cat
            torch.cat(idxp).cpu().numpy(),   # stack, not cat
            torch.cat(mt).cpu().numpy(),
            torch.cat(mp).cpu().numpy(),
            np.array(nodes),                   # np.array, not torch.cat
        )
    finally:
        model_declared.train()

def compute_metrics(metrics_list):
    idxt, idxp = metrics_list.idxt, metrics_list.idxp
    mt, mp = metrics_list.mt, metrics_list.mp
    ppv = np.mean(np.array(idxt) == np.array(idxp))
    if np.std(mt) == 0 or np.std(mp) == 0:
        log.warning("Cannot compute correlation: one input is constant.")
        correlationP = correlationS = float('nan')
    else:
        correlationP, _ = pearsonr(mt.flatten(), mp.flatten())
        correlationS, _ = spearmanr(mt.flatten(), mp.flatten())
    return PerformanceResult(ppv, correlationP, correlationS)

def evaluate_split(data, model, device, batch_size):
    loader = DataLoader(data, batch_size, shuffle=False)
    metrics     = collect_metrics(loader, model, device)
    performance = compute_metrics(metrics)
    return metrics, performance

def seed_fn(seed=42):
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)        
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True, warn_only=True) 


def summarize(model_name, model, tr_size, ev_size, perf_tr, saving_path, elapsed):
    summary = textwrap.dedent(f"""
        \n-----------------------------------------------
        Model name:          {model_name}
        Model declared:\n      {model}
        Training samples:    {tr_size}
        Validation samples (same training data):  {ev_size}

        Pearson Correlation:   {perf_tr.corrP}
        Spearman Correlation:  {perf_tr.corrS}
        Maximum node PPV:      {perf_tr.ppv}

        Running time:   {elapsed:.2f}s ({elapsed/60:.2f} min)
        Results saved:  {saving_path}\n
    """).strip()
    log.info(summary)
    return summary
    
#-------------- Setup model-------------------
class GraphConvModel(nn.Module):
    def __init__(self, in_channels, hidden_channels=64, num_layers=5):
        super().__init__()
        dims = [in_channels] + [hidden_channels] * (num_layers - 1) + [1]
        self.convs = nn.ModuleList(GraphConv(dims[i], dims[i+1]) for i in range(num_layers))
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_weights
        for conv in self.convs[:-1]:
            x = F.relu(conv(x, edge_index, edge_weight))
        return torch.sigmoid(self.convs[-1](x, edge_index, edge_weight))

# Model hyperparameters
variants_table = pd.read_csv(f'/home/mriveraceron/glv-research/tuning_results/91074c4e25b4/tuning_results.csv').sort_values(by="pearson_corr", ascending=False)
hyperparams = variants_table.iloc[0][['channels', 'layers', 'learning_rate', 'epochs']].to_frame().T
log.info(f"Selected hyperparameters for new variant: {hyperparams.to_dict('records')[0]}")

# Device and loss
device  = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
loss_fn = nn.MSELoss()

# Hyperparameters
channels   = int(hyperparams['channels'].iloc[0])
layers     = int(hyperparams['layers'].iloc[0])
lr         = hyperparams['learning_rate'].iloc[0]
epochs     = int(hyperparams['epochs'].iloc[0])   # removed the hardcoded override
batch_size = 50

#----------------------- Generate data--------------------
# Data for training
data_dir = '/home/mriveraceron/glv-research/data_tensors/gnn_proof_full'

paths = glob.glob(f'{data_dir}/*.pt')
data_list = []
for path in paths:
    data_list.extend(torch.load(path, weights_only=False))
    

log.info(f"Loaded data from: {data_dir} | Total samples: {len(data_list)}")
tr_loop = partial(training_fn,
		device         = device,
		loss_fn        = loss_fn,
		epochs         = epochs,
		batch_size     = batch_size,
	)
METHOD_MAP = {'dataloader': 'DL', 'normal': 'normal'}

def wrapper(method='dataloader', suffix = None, shuffle=True):
    seed_fn(42)
    model = GraphConvModel(in_channels=13, hidden_channels=channels, num_layers=layers).to(device)
    loss_history, total_elapsed, elapsed_history = tr_loop(
        model=model,
        weights_path=f'{results_dir}/weights_{suffix}.pth',
        data_train=data_list,
        optimizer=optim.Adam(model.parameters(), lr=lr),
        method=METHOD_MAP[method],
        shuffle=shuffle,
    )
    log.info(f'Finished training in {total_elapsed:.2f}s')
    model_metrics, model_performance = evaluate_split(data_list, model, device, batch_size)
    np.savez(f'{results_dir}/preds_{suffix}.npz',
        **model_metrics._asdict(),
        loss_history=loss_history,
        epochs_history=elapsed_history,
    )
    txt_summary = summarize(suffix, model, len(data_list), len(data_list), model_performance, results_dir, total_elapsed)
    path = f'{results_dir}/summary_notes.txt'
    with open(path, 'a') as f:
        f.write(txt_summary)

wrapper(method='dataloader', suffix = 'dl_shuffled', shuffle=True)
wrapper(method='dataloader', suffix = 'dl_unshuffled', shuffle=False)
wrapper(method='normal', suffix = 'std_ordered') 