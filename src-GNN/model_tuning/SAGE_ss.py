"""
Sample size and its effect on performance of SAGE model.

Previously hyperparameters were already tested and the best one was selected.
After determining the hyperparameters, sample size was tested with the expected tendency that more samples increased the performance of the model.

Therefore, same model will be trained at variable sample sizes.

16-March-2026
"""

# Imports
import os
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"  # Must be before torch import
import logging
import importlib.metadata as metadata
import sys
from collections import namedtuple
import torch
from torch_geometric.nn import SAGEConv
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


# Create result directory
results_dir = '/home/mriveraceron/glv-research/tuning_results/SAGE_ss_V2'
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
    return data_list

train_data = data_generator(data_dir, split='train')
eval_data = data_generator(data_dir, split='eval')

log.info(f"Loaded data from: {data_dir} | Training samples: {len(train_data)} | Evaluation samples: {len(eval_data)}")

#-------------------------------
# Section: Generate grid
#-------------------------------
size = [1000,10000,25000,40000, len(train_data)] 
names = [f'ss_{i}' for i in range(1, len(size)+1)]

# Create a datafrane
tuning_df = pd.DataFrame({
    'model_id': names,
    'train_size': size,
    'elapsed(s)':None,
    'mem_usage(mb)': None,
    'channels': 64,
    'layers': 5,
    'learning_rate': 1e-03,
    'epochs': 700,
    'eval_size': None,
    'ppv_idx': None,
    'pearson_corr': None,
    'spearman_corr': None
})


#-------------------------------
# Section: Declare model
#-------------------------------
class SageModel(nn.Module):
    def __init__(self, hidden_channels=64, num_layers=5):
        super().__init__()
        self.convs = nn.ModuleList()
        # First layer: 1 -> hidden_channels
        self.convs.append(SAGEConv(13, hidden_channels))
        # Middle layers: hidden_channels -> hidden_channels
        for _ in range(num_layers - 2):
            #self.convs.append(GATConv(hidden_channels*heads, hidden_channels, heads=heads))
            self.convs.append(SAGEConv(hidden_channels, hidden_channels))
        # Last layer: hidden_channels -> 1
        self.convs.append(SAGEConv(hidden_channels, 1))
    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        # Apply all layers except the last
        for i, conv in enumerate(self.convs[:-1]):
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
def training_fn(model, device, data_train, weights_path, loss_fn, optimizer, epochs, batch_size=30):
    model.train()
    log.info(f'Training with DataLoader method \n')
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
channels = 64
layers = 5
lr = 1e-05
epochs = 700 
batch_size = 40
eval_loader = DataLoader(eval_data, batch_size=batch_size, shuffle=False)

tr_loop = partial(training_fn,
    device         = device,
    loss_fn        = loss_fn,
    epochs         = epochs,
    batch_size     = batch_size,
)


for i, row in tuning_df.iterrows():
    #-------------------------
    # Declare model hyperparameters
    #-------------------------
    # row = tuning_df.iloc[0]
    size = row['train_size']
    model_name = row['model_id']
    log.info(f'>> Starting model {model_name} with train size {size}.')
    # Create result directory
    variant_dir = f'{results_dir}/{model_name}'
    log.info(f'>> Variant results will be saved at: {variant_dir}')
    os.makedirs(variant_dir, exist_ok=True)
    #-------------------------
    # Section: Run model
    #-------------------------
    # Seeding function
    n_seed = 42
    seed_fn(seed=n_seed)
    # Model parameters
    data = random.sample(train_data, size) 
    model = SageModel(hidden_channels=channels, num_layers=layers).to(device)
    loss_history, total_elapsed = tr_loop(
        model=model,
        weights_path=f'{variant_dir}/model_weights.pth',
        data_train=data,
        optimizer=optim.Adam(model.parameters(), lr=lr),
    )
    log.info(f'Finished training in {total_elapsed:.2f}s')
    # Evaluate model
    model_metrics, model_performance = evaluate_split(eval_loader, model, device)
    # Generate summary
    summary = summarize(model_name = model_name, model = model, 
                        tr_size = size, ev_size = len(eval_data), 
                        perf_tr = model_performance, saving_path = variant_dir, elapsed = total_elapsed)
    #------------------------
    # Section: Save metrics result_exp_dir
    #------------------------
    np.savez(f'{variant_dir}/prediction_values.npz',
        idxt  = model_metrics.idxt,
        idxp  = model_metrics.idxp,
        mt   = model_metrics.mt,
        mp   = model_metrics.mp,
        nodes = model_metrics.nodes,
        loss_history  = loss_history
    )
    #------------------------
    # Add results to dataframe
    #------------------------
    total_bytes = len(pickle.dumps(data))
    tuning_df.loc[i, :] = {
        **tuning_df.loc[i].to_dict(),  # preserve existing values
        'ppv_idx':        model_performance.ppv,
        'pearson_corr':   model_performance.corrP,
        'spearman_corr':  model_performance.corrS,
        'mem_usage(mb)':  total_bytes / 1e6,
        'eval_size':      len(eval_data),
        'elapsed(s)':     total_elapsed,
    }
    # Save table every row
    tuning_df.to_csv(f'{results_dir}/tuning_results.csv', index=False)
    
