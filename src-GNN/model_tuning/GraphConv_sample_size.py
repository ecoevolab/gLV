# Training script for random networks.
# 

import logging
import sys
import os 

# Create result directory
results_dir = '/home/mriveraceron/glv-research/tuning_results/GraphConv_sample_size'
os.makedirs(results_dir, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
    handlers=[
        logging.FileHandler(f'{results_dir}/run_log.txt', mode='w'),
        logging.StreamHandler(sys.stdout)
    ]
)
log = logging.getLogger(__name__)

#-------------------------------
# Section: Generate grid
#-------------------------------
"""
    We will be training model with random networks.
    Model to train will be Variant 3 hyperparameters at different sizes.
"""
import pandas as pd

size = [1000] + list(range(5000, 40000, 5000))
names = [f'Variant_{i}' for i in range(1, len(size)+1)]

# Create a datafrane
tuning_df = pd.DataFrame({
    'model_id': names,
    'train_size': size,
    'mem_usage(mb)': None,
    'channels': 64,
    'layers': 5,
    'learning_rate': 1e-03,
    'epochs': 700,
    'eval_size': None,
    'accuracy_idx': None,
    'pearson_corr': None,
    'spearman_corr': None
})

#-------------------------------
# Section: Declare model
#-------------------------------
import torch 
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv

class GCNModel(nn.Module):
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
# Section: Evaluation function at last epoch
#-------------------------------
from collections import namedtuple
from scipy.stats import pearsonr
from scipy.stats import spearmanr

MetricsResult     = namedtuple('MetricsResult', ['idxt', 'idxp', 'mt', 'mp'])
PerformanceResult = namedtuple('PerformanceResult', ['acc', 'corrP', 'corrS'])

def collect_metrics(loader, model_declared, device):
    idxt, idxp, mt, mp = [], [], [], []
    try:
        model_declared.eval()
        with torch.no_grad():
            for batch in loader:
                batch = batch.to(device)
                out = model_declared(batch)
                y_list = unbatch(batch.y, batch.batch)
                out_list = unbatch(out, batch.batch)
                for y, o in zip(y_list, out_list):
                    idxt.append(torch.argmax(y, dim=0))
                    idxp.append(torch.argmax(o, dim=0))
                    mt.append(y)
                    mp.append(o)
        # Convert to arrays
        idxt = torch.cat(idxt).cpu().numpy()
        idxp = torch.cat(idxp).cpu().numpy()
        mt = torch.cat(mt).cpu().numpy()
        mp = torch.cat(mp).cpu().numpy()
        return MetricsResult(idxt, idxp, mt, mp)
    finally:
        model_declared.train()

def compute_metrics(metrics_list):
    idxt, idxp = metrics_list.idxt, metrics_list.idxp
    mt, mp = metrics_list.mt, metrics_list.mp
    accuracy = np.mean(np.array(idxt) == np.array(idxp))
    if np.std(mt) == 0 or np.std(mp) == 0:
        log.warning("Cannot compute correlation: one input is constant.")
        correlationP = correlationS = float('nan')
    else:
        correlationP, _ = pearsonr(mt.flatten(), mp.flatten())
        correlationS, _ = spearmanr(mt.flatten(), mp.flatten())
    return PerformanceResult(accuracy, correlationP, correlationS)

#-------------------------------
# Section: Seeding function
#-------------------------------
import torch.optim as optim
import numpy as np
import random

def seed_fn(seed=42):
    # Set ALL seeds for full reproducibility
    torch.manual_seed(seed)                 # Seed CPU 
    torch.cuda.manual_seed(seed)            # Seed GPU
    np.random.seed(seed)                    # Seed numpy
    random.seed(seed)                       # Seed python random
    torch.backends.cudnn.deterministic = True   # Ensure deterministic behavior
    torch.backends.cudnn.benchmark = False 

#-------------------------------
# Section: Training function 
#-------------------------------
import glob
import time
from tqdm import tqdm
from torch_geometric.utils import unbatch
from torch_geometric.loader import DataLoader

def training_fn(model_declared, device, data_train, data_eval, weights_path, loss_fn, optimizer, epochs, eval_every=50, patience=2, batch_size=30):
    #------------------------------------------
    # Section: Declare training variables
    model_declared.train()
    loss_history  =  []          # Loss at epoch
    total_elapsed = 0            # Running time
    stop_early = False           # Early stopping flag
    metrics_list, performance_list = None, None            
    # Early stopping variables
    best_model = None          # Best model weights
    best_loss = float('inf')   # Best loss
    no_improve = 0             # counter for no improvement
    #---------------------
    # Section: Create batches of data
    loader_train = DataLoader(data_train, batch_size, shuffle=True)
    loader_eval = DataLoader(data_eval, batch_size, shuffle=False)  # Important: no shuffling for evaluation, to keep order for metrics
    for epoch in tqdm(range(epochs), desc="Training"):
        start = time.time()
        epoch_loss = 0
        #--------------------------
        for batch in loader_train:
            #----------------------
            # Move it to device and run model
            # data = data_list[0]
            batch = batch.to(device)
            optimizer.zero_grad()
            out = model_declared(batch)
            loss = loss_fn(out, batch.y)
            # Verify if loss is finite
            if torch.isnan(loss):
                log.error(f"NaN loss at epoch {epoch}, aborting.")
                raise ValueError(f"NaN loss detected at epoch {epoch}")
            loss.backward()         # Backpropagation, compute gradients
            torch.nn.utils.clip_grad_norm_(model_declared.parameters(), max_norm=1.0) # Gradient clipping
            optimizer.step()            # Update weights
            epoch_loss += loss.item()   # Accumulate loss
        #----------------------
        # Append epoch loss to history
        #----------------------
        loss_history.append(epoch_loss)
        # loss_history = np.append(loss_history, epoch_loss)
        elapsed = time.time() - start
        total_elapsed += elapsed
        # Print every n epochs
        if epoch % 10 == 0:
            log.info(f"Epoch {epoch}: Loss = {epoch_loss}, Elapsed time: {elapsed:.2f}")
        #----------------------
        # Section: Early stopping if evaluation loss does not improve for patience epochs
        #----------------------
        # Evalaution every 50 epochs
        if epoch % eval_every == 0:
            log.info(f"Checking early stopping at epoch {epoch} ")
            eval_loss = 0
            try:        
                model_declared.eval()
                with torch.no_grad():
                    for batch in loader_eval:
                        batch = batch.to(device)
                        out = model_declared(batch)
                        loss = loss_fn(out, batch.y)
                        eval_loss += loss.item()
                log.info(f"Evaluation loss at epoch {epoch}: {eval_loss}")
            finally:
                model_declared.train()
            # If evaluation loss does not improve, increase counter
            if eval_loss <= best_loss:
                best_loss = eval_loss
                no_improve = 0
                best_model = model_declared.state_dict()  # Save best model weights
            else :
                no_improve += 1
            if no_improve >= patience:
                log.info(f"Early stopping at epoch {epoch}, no improvement for {eval_every * patience} epochs.")
                stop_early = True
        #----------------------
        # Section: Evaluate model
        #----------------------
        # Model stops or is last epoch
        if stop_early or epoch == epochs - 1:
            loss_history = np.array(loss_history)
            metrics_list = collect_metrics(loader_eval, model_declared, device)
            performance_list = compute_metrics(metrics_list)
            # Save weights
            torch.save({
                'epoch': epoch,
                'model_state_dict': best_model if best_model is not None else model_declared.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'loss': epoch_loss
            }, weights_path)
            break
    # Summary
    log.info(f'>> the total elapsed time with {epochs} epochs is {total_elapsed:.2f} seconds ( {total_elapsed/60:.2f} minutes)')   
    return  loss_history, metrics_list, performance_list, total_elapsed

#-------------------------------
# Section: Function to save Summary
#-------------------------------
def summarize(model_declared, optimizer, row, train_dirs, eval_dirs, performance_list, result_exp_dir, extra_info):
    summary = f"""
    Model Training Summary
    =========================
    Model variant: {row['model_id']}
    Model: {model_declared}
    Samples for training {row['train_size']}
    Optimizer LR:   {optimizer.param_groups[0]['lr']}
    Number of epochs: {row['epochs']}
    Model layers: {row['layers']}
    Model hidden channels: {row['channels']}
    -----------------------------------------------
    Seed: {extra_info.n_seed}
    Evaluation interval: {extra_info.eval_interval} epochs
    Early stopping patience: {extra_info.patience} evaluations
    DataLoaders batch size: {extra_info.batch_size}
    Training data paths: \n\t{train_dirs}
    -----------------------------------------------
    Validation data path: {eval_dirs}
    Validation samples: {extra_info.validation_samples}
    Pearson Correlation:  {performance_list.corrP}    
    Spearman Correlation: {performance_list.corrS}   
    Maximum node accuracy: {performance_list.acc}  
    Running time seconds: {extra_info.total_elapsed}
    Epochs performed: {extra_info.epochs_runned}
    """
    with open(f'{result_exp_dir}/training_summary.txt', 'w') as f:
        f.write(summary)
    print(summary)
    log.info(summary)

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
    return data_list, paths

train_data, train_paths = data_generator(data_dir, split='train')
eval_data, eval_paths = data_generator(data_dir, split='eval')

#-------------------------------
# Section: Run model
#-------------------------------
import glob
import pickle


# Define namedtuple for results
ExtraInfo         = namedtuple('ExtraInfo', ['epochs_runned', 'n_seed', 'total_elapsed', 'validation_samples', 'eval_interval', 'patience', 'batch_size'])

# Get directories of data for training
train_dirs = '\n'.join(f'  {p}' for p in set(os.path.dirname(p) for p in train_paths))
eval_dirs  = '\n'.join(f'  {p}' for p in set(os.path.dirname(p) for p in eval_paths))

# Constant model parameters
device  = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
loss_fn = nn.MSELoss()

for i, row in tuning_df.iterrows():
    #-------------------------
    # Declare model hyperparameters
    #-------------------------
    # row = tuning_df.iloc[0]
    size = row['train_size']
    lr = row['learning_rate']
    model_name = row['model_id']
    epochs = row['epochs']
    channels = int(row['channels'])
    layers = int(row['layers'])
    eval_interval = 50
    patience = 2
    batch_size = 30
    log.info(f'>> Starting model {model_name} with train size {size} and learning rate {lr}')
    #-------------------------
    # Section: Run model
    #-------------------------
    # Seeding function
    n_seed = 42
    seed_fn(seed=n_seed)
    # Model parameters
    model_declared = GCNModel(hidden_channels=channels,num_layers=layers).to(device)
    optimizer = optim.Adam(model_declared.parameters(), lr=lr)
    # Create result directory
    result_exp_dir = f'{results_dir}/{model_name}'
    log.info(f'>> Variant results will be saved at: {result_exp_dir}')
    os.makedirs(result_exp_dir, exist_ok=True)
    # Slice data for training
    random.shuffle(train_data)
    data_to_train = train_data[:size]
    # Run training function
    weights_path = f'{result_exp_dir}/model_weights.pth'
    try:
        loss_history, metrics_list, performance_list, total_elapsed = training_fn(model_declared, device, data_to_train, eval_data, weights_path, loss_fn, optimizer, epochs, eval_interval, patience, batch_size)
    except ValueError as e:
        log.info(f"Error occurred: {e}")
        continue
    #------------------------
    # Section: Generate summary
    #------------------------
    # ['epochs_runned', 'n_seed', 'total_elapsed', 'validation_samples', 'eval_interval', 'patience', 'batch_size']
    extra_info = ExtraInfo(len(loss_history), n_seed, total_elapsed, len(eval_data), eval_interval, patience, batch_size)
    summarize(model_declared, optimizer, row, train_dirs, eval_dirs, performance_list, result_exp_dir, extra_info)
    #------------------------
    # Section: Save metrics result_exp_dir
    #------------------------
    np.savez(f'{result_exp_dir}/metric-values.npz',
        max_idx_true  = metrics_list.idxt,
        max_idx_pred  = metrics_list.idxp,
        values_true   = metrics_list.mt,
        values_pred   = metrics_list.mp,
        loss_history  = loss_history
    )
    #------------------------
    # Add results to dataframe
    #------------------------
    total_bytes = len(pickle.dumps(data_to_train))
    tuning_df.loc[i, 'accuracy_idx'] = performance_list.acc
    tuning_df.loc[i, 'pearson_corr'] = performance_list.corrP
    tuning_df.loc[i, 'spearman_corr'] = performance_list.corrS
    tuning_df.loc[i, 'mem_usage(mb)'] = total_bytes / 1e6
    tuning_df.loc[i, 'eval_size'] = len(eval_data)
    # Save table every row
    tuning_df.to_csv(f'{results_dir}/tuning_results.csv', index=False)
    