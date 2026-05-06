"""
# Random Network Feature Analysis
# =================================
# Evaluates the contribution of node features and adjacency matrix structure
# to model performance using a previously generated random interaction network 
# dataset.
#
# Experiments:
#   - Node feature importance: dataset is re-run with dummy/features.
#   - Adjacency matrix importance: matrix will be filtered for surviving 
#     nodes or kept the same.
"""

import os
import sys
from glob import glob
import time
import random
import logging
import pickle
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
results_dir = '/home/mriveraceron/glv-research/Results/feats_analysis'
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

# Generate results directory
log.info(f'>> Results will be saved at: {results_dir}')

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
data_dir = '/home/mriveraceron/glv-research/data_null'
dirs = glob(f'{data_dir}/*_train/') + glob(f'{data_dir}/*_eval/')

def generate_data(data_dir, split_ratio=0.8):
    log.info(f'>> Generating data from: {data_dir}')
    paths = glob(f'{data_dir}/*.pt')
    data_list = [d for path in paths for d in torch.load(path, weights_only=False)]
    # Build dummy version from the same loaded data (no re-read)
    dummy_list = []
    for d in data_list:
        if d.edge_weights.shape[0] > d.x.shape[0] ** 2:
            log.info(f'Error in data, more edges than maximum theoretical edges.\n')
            continue  
        d_dummy = d.clone()
        d_dummy.x = torch.ones(d.x.shape[0], 1, dtype=torch.float32)
        dummy_list.append(d_dummy)
    split = int(len(data_list) * split_ratio)
    log.info(f"Loaded data from: {data_dir} | Total samples: {len(data_list)}")
    log.info(f"Train: {split} | Val: {len(data_list) - split}")
    return data_list[:split], data_list[split:], dummy_list[:split], dummy_list[split:]

train_feats, eval_feats, train_dummy, eval_dummy = generate_data(dirs[0])


#-------------------------------
# Section: Run model
#-------------------------------

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
run_configs = [
    ('features_run', train_feats, eval_feats, False),
    ('dummy_run', train_dummy, eval_dummy, True),
]


df = pd.DataFrame({
    'model':            [name for name, *_ in run_configs],
    'train_size':       None,
    'eval_size':        None,
    'tr_accuracy_idx':  None,
    'tr_pearson_corr':  None,
    'tr_spearman_corr': None,
    'ev_accuracy_idx':  None,
    'ev_pearson_corr':  None,
    'ev_spearman_corr': None,
    'elapsed_seconds':  None,
})

start_line = f"""
Adjacency matrix analysis
=========================
Hyperparameters for all models
Hidden channels: {channels}
Model layers: {layers},
Optimizer lr:   {lr}
Number of epochs: {epochs}
Seed used: {n_seed}
"""
with open(f'{results_dir}/training_summary.txt', 'w') as f:
    f.write(start_line)

def summarize(model_name, model_declared, tr_size, ev_size, perf_tr, perf_ev, saving_path, elapsed):
    summary = f"""
    \n-----------------------------------------------
    Model name: {model_name}
    Model declared: {model_declared}
    Training samples: {tr_size}
    Validation samples: {ev_size} \n
    Pearson Correlation training:  {perf_tr.corrP}    
    Spearman Correlation training: {perf_tr.corrS}   
    Maximum node ppv training: {perf_tr.ppv} \n
    Pearson Correlation validation:  {perf_ev.corrP}    
    Spearman Correlation validation: {perf_ev.corrS}   
    Maximum node ppv validation: {perf_ev.ppv} \n
    Running time seconds: {elapsed}
    Results saved at: {saving_path} \n
    """
    with open(f'{results_dir}/training_summary.txt', 'a') as f:
        f.write(summary)
    log.info(summary)

# Testing line
model_name, data_train, data_ev, dummy = run_configs[0]

for model_name, data_train, data_ev, dummy in run_configs:
    # Setup
    save_path    = os.path.join(results_dir, model_name)
    weights_path = os.path.join(save_path, 'model_weights.pth')
    os.makedirs(save_path, exist_ok=True)
    # Select model
    in_channels = 1 if dummy else 13
    data_label  = 'dummy' if dummy else 'features'
    log.info(f'Starting training with {data_label} data.')
    model = GraphConvModel(in_channels=in_channels, hidden_channels=channels, num_layers=layers).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    tr_size, ev_size = len(data_train), len(data_ev)
    log.info(f'>> [{model_name}] train_size={tr_size} | eval_size={ev_size} | lr={lr} | epochs={epochs} | saving to {save_path}')
    try:
        random.shuffle(data_train)
        loss_history, total_elapsed = tr_loop(
            model=model,
            weights_path=weights_path,
            data_train=data_train,
            optimizer=optimizer,
        )
        # Generate metrics and performance
        m_tr, p_tr = evaluate_split(data_train, model, device, batch_size)
        m_ev, p_ev = evaluate_split(data_ev,    model, device, batch_size)
        # Update results table
        mask = df['model'] == model_name  # ← point 4: cache mask
        df.loc[mask, [
            'train_size', 'tr_accuracy_idx', 'tr_pearson_corr', 'tr_spearman_corr',
            'eval_size',  'ev_accuracy_idx', 'ev_pearson_corr', 'ev_spearman_corr',
            'elapsed_seconds',
        ]] = [
            tr_size, p_tr.ppv, p_tr.corrP, p_tr.corrS,
            ev_size, p_ev.ppv, p_ev.corrP, p_ev.corrS,
            total_elapsed,
        ]
       # Save predictions
        values_path = os.path.join(save_path, 'model_preds.npz')
        np.savez(values_path,
            # Trainin data
            idxt_train  = m_tr.idxt,
            idxp_train  = m_tr.idxp,
            mt_train    = m_tr.mt,
            mp_train    = m_tr.mp,
            nodes_train = m_tr.nodes,
            loss        = loss_history,
            # Eval data
            idxt_eval   = m_ev.idxt,
            idxp_eval   = m_ev.idxp,
            mt_eval     = m_ev.mt,
            mp_eval     = m_ev.mp,
            nodes_eval  = m_ev.nodes,
        )
        # Add model to summary
        summarize(
            model_name=model_name, model_declared=model,
            tr_size=tr_size,       ev_size=ev_size,
            perf_tr=p_tr,          perf_ev=p_ev,
            saving_path=save_path, elapsed=total_elapsed,
        )
    except ValueError as e:
        log.error(f"[{model_name}] Training aborted: {e}")

df.to_csv(f'{results_dir}/result_table.csv', index=False)
pd.read_csv(f'{results_dir}/result_table.csv')

