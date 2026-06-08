
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
import pickle
from itertools import product

# Create result directory
results_dir = '/home/mriveraceron/glv-research/updated_results/GraphConv_hpo'
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
    

#-------------------------
# Section: Generate data
#-------------------------
# Data for training
def data_generator(data_path, split=5000):
    data_list = []
    paths = glob.glob(f'{data_path}/*.pt')
    if not paths:
        raise FileNotFoundError(f"No .pt files found under {data_path}/")
    for path in paths:
        data = torch.load(path, weights_only=False)
        data_list.extend(data)
    random.seed(42)  # For reproducibility
    random.shuffle(data_list)
    # log.info(f"Total samples generated: {len(data_list)}")
    return data_list[:split]

# Generate training data
data_path = '/home/mriveraceron/glv-research/data_null/d2f93775a813_train'
train_data = data_generator(data_path = data_path, split = 5000)
log.info(f"Loaded data from: {data_path} | Training samples: {len(train_data)}")

# Generate evaluation data
data_path = '/home/mriveraceron/glv-research/data_null/efa83c9fafa0_eval'
eval_data = data_generator(data_path, split = 1000)
log.info(f"Loaded data from: {data_path} | Evaluation samples: {len(eval_data)}")


#-------------------------------
# Section: Generate grid
#-------------------------------
layers = [5,10]
LearnR = [1e-1, 1e-3, 1e-5, 1e-7]
channels = [64,128]
pairs = list(product(layers, LearnR, channels))
names = [f'Variant_{i}' for i in range(1, len(pairs)+1)]

# Create a datafrane
tuning_df = pd.DataFrame({
    'model_id': names,
    'tr_size':None,
    'channels': [l[2] for l in pairs],
    'layers': [l[0] for l in pairs],
    'learning_rate': [l[1] for l in pairs],
    'epochs': 700,
    'elapsed(s)': None,
    'eval_size': None,
    'ppv_idx': None,
    'pearson_corr': None,
    'spearman_corr': None,
})

#-------------------------------
# Section: Declare model
#-------------------------------
class GraphConv_Model(nn.Module):
    def __init__(self, hidden_channels=64, num_layers=5):
        super().__init__()
        self.convs = nn.ModuleList()
        # First layer: 1 -> hidden_channels
        self.convs.append(GraphConv(13, hidden_channels))
        # Middle layers: hidden_channels -> hidden_channels
        for _ in range(num_layers - 2):
            self.convs.append(GraphConv(hidden_channels, hidden_channels))
        # Last layer: hidden_channels -> 1
        self.convs.append(GraphConv(hidden_channels, 1))
    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        # Apply all layers except the last
        for conv in self.convs[:-1]:
            x = conv(x, edge_index)
            x = F.relu(x)
        # Apply last layer with sigmoid
        x = self.convs[-1](x, edge_index)
        x = torch.sigmoid(x)
        return x  # [num_nodes]


#-------------------------------
# Section: Evaluation functions
#-------------------------------

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
                for y, o in zip(unbatch(batch.y, batch.batch), unbatch(out, batch.batch)):
                    idxt.append(torch.argmax(y, dim=0))   # 0-dim tensor
                    idxp.append(torch.argmax(o, dim=0))   # 0-dim tensor
                    mt.append(y)
                    mp.append(o)
                    nodes.append(y.shape[0])              # plain int
    finally:
        model_declared.train()
    return MetricsResult(
            torch.cat(idxt).cpu().numpy(),   # stack, not cat
            torch.cat(idxp).cpu().numpy(),   # stack, not cat
            torch.cat(mt).cpu().numpy(),
            torch.cat(mp).cpu().numpy(),
            np.array(nodes),                   # np.array, not torch.cat
        )

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

def evaluate_split(eval_loader, model, device):
    metrics     = collect_metrics(eval_loader, model, device)
    performance = compute_metrics(metrics)
    return metrics, performance

#-------------------------------
# Section: Seeding function
#-------------------------------
def seed_fn(seed=42):
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)        
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True, warn_only=True) 


#-------------------------------
# Section: Summary function
#-------------------------------
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

#-------------------------------
# Section: Training function 
#-------------------------------
def training_fn(model, model_name, device, data_train, weights_path, loss_fn, optimizer, epochs, batch_size=30):
    model.train()
    log.info(f'Starting training of model {model_name} \n')
    data = DataLoader(data_train, batch_size, shuffle = True)
    loss_history = []
    total_elapsed = 0
    for epoch in tqdm(range(epochs), desc="Training"):
        start = time.time()
        epoch_loss = 0
        for d in data:
            d = d.to(device)
            # optimizer.zero_grad()
            optimizer.zero_grad(set_to_none=True)
            out = model(d)
            loss = loss_fn(out, d.y)
            if torch.isnan(loss):
                raise ValueError(f"NaN loss detected at epoch {epoch}")
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            # epoch_loss += loss.item()
            epoch_loss += loss.detach()
        epoch_elapsed = time.time() - start		# Epoch elapsed time
        total_elapsed += epoch_elapsed			# Add it to model elapsed time
        loss_history.append(epoch_loss.item())			# Append loss
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
    return np.array(loss_history), total_elapsed


#-------------------------------
# Section: Run model
#-------------------------------

# Constant model parameters
device  = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
loss_fn = nn.MSELoss()
epochs = 700 
batch_size = 40
tr_size = len(train_data)
ev_size = len(eval_data)
eval_loader = DataLoader(eval_data, batch_size=batch_size, shuffle=False)

# Shuffle training data and declare function
tr_loop = partial(training_fn,
    data_train     = train_data,
    device         = device,
    loss_fn        = loss_fn,
    epochs         = epochs,
    batch_size     = batch_size,
)

# Pre-define result columns to update (avoids repeated string literals)
RESULT_COLS = ['ppv_idx', 'pearson_corr', 'spearman_corr', 'eval_size', 'tr_size', 'elapsed(s)']
RESULTS_CSV = f'{results_dir}/tuning_results.csv'

for i, row in tuning_df.iterrows():
    # ── Unpack config ──────────────────────────────────────────────
    lr         = row['learning_rate']
    model_name = row['model_id']
    epochs     = row['epochs']
    channels   = int(row['channels'])
    layers     = int(row['layers'])
    save_dir   = f'{results_dir}/{model_name}'
    log.info(f'>> Starting model {model_name} | Train size: {tr_size}')
    # ── Setup ──────────────────────────────────────────────────────
    os.makedirs(save_dir, exist_ok=True)
    log.info(f'>> Results will be saved at: {save_dir}')
    seed_fn(seed=42)
    model = GraphConv_Model(hidden_channels=channels, num_layers=layers).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    # ── Train ──────────────────────────────────────────────────────
    loss_history, total_elapsed = tr_loop(
        model=model,
        model_name = model_name,
        weights_path=f'{save_dir}/model_weights.pth',
        optimizer=optimizer,
    )
    log.info(f'Finished training in {total_elapsed:.2f}s')
    # ── Evaluate ───────────────────────────────────────────────────
    model_metrics, model_performance = evaluate_split(eval_loader, model, device)
    summarize(
        model_name=model_name, model=model,
        tr_size=tr_size, ev_size=ev_size,
        perf_tr=model_performance,
        saving_path=save_dir, elapsed=total_elapsed,
    )
    # ── Save predictions ───────────────────────────────────────────
    preds_path = f'{save_dir}/prediction_values.npz'
    np.savez(
        preds_path,
        idxt=model_metrics.idxt,
        idxp=model_metrics.idxp,
        mt=model_metrics.mt,
        mp=model_metrics.mp,
        nodes=model_metrics.nodes,
        loss_history=loss_history,
    )
    log.info(f'>> Predictions saved at: {preds_path}')
    # ── Update results dataframe ───────────────────────────────────
    tuning_df.loc[i, RESULT_COLS] = [
        model_performance.ppv,
        model_performance.corrP,
        model_performance.corrS,
        ev_size,
        tr_size,
        total_elapsed,
    ]
    tuning_df.to_csv(RESULTS_CSV, index=False)