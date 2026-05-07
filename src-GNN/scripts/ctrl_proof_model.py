"""
Generate proof of concept models

2-May-2026
"""

# Imports
import os
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"  # Must be before torch import
import sys
import shutil
import glob
import time
import random
import logging
from collections import namedtuple
from functools import partial
import importlib.metadata as metadata

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from tqdm import tqdm

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GraphConv
from torch_geometric.utils import unbatch

#-----------------------
# Set up logging
#-----------------------
results_dir = '/home/mriveraceron/glv-research/Results/gnn_proof'
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


#-------------------------------
# Section: Declare model
#-------------------------------
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


def training_fn(model, device, data_train, weights_path, loss_fn, optimizer, epochs, batch_size=30):
    model.train()
    loader_train = DataLoader(data_train, batch_size, shuffle=True)
    loss_history = []
    total_elapsed = 0
    for epoch in tqdm(range(epochs), desc="Training"):
        start = time.time()
        epoch_loss = 0
        for batch in loader_train:
            batch = batch.to(device)
            optimizer.zero_grad()
            out = model(batch)
            loss = loss_fn(out, batch.y)
            if torch.isnan(loss):
                raise ValueError(f"NaN loss detected at epoch {epoch}")
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            epoch_loss += loss.item()
        elapsed = time.time() - start
        total_elapsed += elapsed
        loss_history.append(epoch_loss)
        if epoch % 100 == 0:
            log.info(f"Epoch {epoch}: Loss={epoch_loss:.6f} | Time={elapsed:.2f}s")
    # Save after all epochs
    torch.save({
        'epoch': epochs - 1,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'loss': epoch_loss
    }, weights_path)
    log.info(f"Total elapsed: {total_elapsed:.2f}s ({total_elapsed/60:.2f} min)")
    return np.array(loss_history), total_elapsed

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
    # Set ALL seeds for full reproducibility
    torch.manual_seed(seed)                 # Seed CPU 
    torch.cuda.manual_seed(seed)            # Seed GPU
    np.random.seed(seed)                    # Seed numpy
    random.seed(seed)                       # Seed python random
    torch.backends.cudnn.deterministic = True   # Ensure deterministic behavior
    torch.backends.cudnn.benchmark = False 

#-----------------------
# Load data
#-----------------------
# Data for training
full_dir = '/home/mriveraceron/glv-research/data_tensors/gnn_proof_full'
sub_dir = '/home/mriveraceron/glv-research/data_tensors/gnn_proof_sub'
# data_dir = '/home/mriveraceron/glv-research/data_tensors/Boosted_filtered'
# Communities were pre-filtered so that extinctions and network topology were derived from survival nodes only

def generate_data(data_dir, dummy=False):
    paths = glob.glob(f'{data_dir}/*.pt')
    data_list = []
    for path in paths:
        data_list.extend(torch.load(path, weights_only=False))
    if dummy:
        for d in data_list:
            d.x = torch.ones(d.x.shape[0], 1, dtype=torch.float32)
            del d.x_dummy
    split = int(len(data_list) * .8)
    log.info(f"Loaded {'dummy ' if dummy else ''}data from: {data_dir} | Total samples: {len(data_list)}")
    return data_list[:split], data_list[split:]


train_full,       eval_full       = generate_data(full_dir)
train_sub,        eval_sub        = generate_data(sub_dir)
train_full_dummy, eval_full_dummy = generate_data(full_dir, dummy=True)
train_sub_dummy,  eval_sub_dummy  = generate_data(sub_dir,  dummy=True)
#-------------------------------
# Section: Run model
#-------------------------------

# Model hyperparameters
hyperparams = {'channels': 64,'layers': 5,'lr': 0.0001,'epochs': 700}
log.info(f"Selected hyperparameters for new variant: {hyperparams}")

# Device and loss
device  = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
loss_fn = nn.MSELoss()

# Hyperparameters
channels   = int(hyperparams['channels'])
layers     = int(hyperparams['layers'])
lr         = hyperparams['lr']
epochs     = hyperparams['epochs']   # removed the hardcoded override
batch_size = 50

# Model setup
n_seed     = 42
seed_fn(seed=n_seed)

tr_loop = partial(training_fn,
    device         = device,
    loss_fn        = loss_fn,
    epochs         = epochs,
    batch_size     = batch_size,
)

#------------------------
# Generate dataframe
#------------------------
df = pd.DataFrame({
    'model':            ['full_feats', 'full_dummy', 'sub_feats', 'sub_dummy'],
    'train_size':       None,
    'eval_size':        None,
    'ppv_idx':  None,
    'pearson_corr':  None,
    'spearman_corr': None,
    'elapsed_seconds':  None,
})

run_configs = [
    ('full_feats', 'weights_full_feats.pth', 'preds_full_feats.npz', train_full, eval_full, False),
    ('full_dummy', 'weights_full_dummy.pth', 'preds_full_dummy.npz', train_full_dummy, eval_full_dummy, True),
    ('sub_feats',  'weights_sub_feats.pth',  'preds_sub_feats.npz',  train_sub, eval_sub, False),
    ('sub_dummy',  'weights_sub_dummy.pth',  'preds_sub_dummy.npz',  train_sub_dummy, eval_sub_dummy, True)
]

start_line = f"""
gnn proof of concept summary
=========================
Hyperparameters for all models
Hidden channels: {hyperparams['channels']}
Model layers: {hyperparams['layers']},
Optimizer lr:   {hyperparams['lr']}
Number of epochs: {hyperparams['epochs']}
Seed used: {n_seed}\n
"""
with open(f'{results_dir}/training_summary.txt', 'w') as f:
    f.write(start_line)


def summarize(model_name, model_declared, tr_size, ev_size, perf_ev, saving_path, elapsed):
    summary = f"""
    \n-----------------------------------------------
    Model name: {model_name}
    Model declared: {model_declared}
    Training samples: {tr_size}
    Validation samples: {ev_size} \n
    Pearson Correlation validation:  {perf_ev.corrP}    
    Spearman Correlation validation: {perf_ev.corrS}   
    Maximum node ppv validation: {perf_ev.ppv} \n
    Running time seconds: {elapsed}
    Results saved at: {saving_path} \n
    """
    with open(f'{results_dir}/training_summary.txt', 'a') as f:
        f.write(summary)
    log.info(summary)


def run_single(model_name, weights_file, npz_file, data_train, data_ev, dummy):
    # Model setup
    in_channels  = 1 if dummy else 13
    model        = GraphConvModel(in_channels=in_channels, hidden_channels=channels, num_layers=layers).to(device)
    optimizer    = optim.Adam(model.parameters(), lr=lr)
    save_path    = f'{results_dir}/{model_name}'
    weights_path = os.path.join(save_path, weights_file)
    tr_size, ev_size = len(data_train), len(data_ev)
    os.makedirs(save_path, exist_ok=True)
    log.info(f'>> Starting {model_name} | train={tr_size} | lr={lr} | epochs={epochs}')
    log.info(f'>> Results → {save_path}')
    loss_history, total_elapsed = tr_loop(
        model=model, weights_path=weights_path,
        data_train=data_train, optimizer=optimizer
    )
    m_ev, p_ev = evaluate_split(data_ev, model, device, batch_size)
    # Update dataframe
    df.loc[df['model'] == model_name, ['train_size', 'eval_size', 'ppv_idx','pearson_corr', 
                                       'spearman_corr', 'elapsed_seconds']
                                       ] = [tr_size, ev_size, p_ev.ppv, p_ev.corrP, p_ev.corrS, total_elapsed]
    # Save arrays
    np.savez(os.path.join(save_path, npz_file), **{
        'idxt': m_ev.idxt, 'idxp': m_ev.idxp,
        'mt'  : m_ev.mt,   'mp'  : m_ev.mp,
        'loss': loss_history, 'nodes': m_ev.nodes,
    })
    summarize(model_name=model_name, model_declared=model, tr_size=tr_size, 
              ev_size=ev_size, perf_ev=p_ev, saving_path=save_path, elapsed=total_elapsed)
    df.to_csv(f'{results_dir}/result_table.csv', index=False)


# --- Loop ---
for config in run_configs:
    model_name, weights_file, npz_file, data_train, data_ev, dummy = config
    try:
        run_single(model_name, weights_file, npz_file, data_train, data_ev, dummy)
        print(df)
    except ValueError as e:
        log.error(f"Training aborted for '{model_name}': {e}")


pd.read_csv(f'{results_dir}/result_table.csv')


