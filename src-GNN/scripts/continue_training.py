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
import sys
import torch
import numpy as np
import torch.optim as optim
import random
import glob
import textwrap
import shutil


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

# Reroute prints to logger file
class StreamToLogger:
    def __init__(self, logger, level=logging.INFO):
        self.logger = logger
        self.level = level
    def write(self, message):
        message = message.strip()
        if message:
            self.logger.log(self.level, message)
    def flush(self):
        pass

# Generate log file and redirect stdout to it. All prints will be in the log file.
log     = make_logger('run_log',  f'{results_dir}/run_log.txt')
sys.stdout = StreamToLogger(log, logging.INFO)
 

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

#-------------------------
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
def summarize( model, tr_size, ev_size, perf_tr, saving_path, elapsed):
    summary = textwrap.dedent(f"""
        \n-----------------------------------------------
        Model name:          {model.__class__.__name__ }
        Model declared:\n      {model}
        Training samples:    {tr_size}
        Validation samples:  {ev_size} \n
        Pearson Correlation:   {perf_tr.corrP}
        Spearman Correlation:  {perf_tr.corrS}
        Maximum node PPV:      {perf_tr.ppv} \n
        Running time:   {elapsed:.2f}s ({elapsed/60:.2f} min)
        Results saved:  \n{saving_path}\n
    """).strip()
    with open(f'{results_dir}/training_summary.txt', 'w') as f:
        f.write(summary)
    print('>> Training summary generated!')


#------------------------------------------------------------
# Section: Source training and evaluation functions.
import sys
sys.path.append('/home/mriveraceron/glv-research/gLV/src-GNN/FUN')
from training_fun import training_fn
from evaluation_fun import collect_metrics, compute_metrics, evaluate_split

#------------------------------------------------------------
# Section: Train model and evaluate it
# Declare model
device  = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = GraphConv_model(in_channels=13, hidden_channels=64, num_layers=5)
model.to(device)
loss_fn = nn.MSELoss()
extra_epochs = 10 
batch_size = 40
optimizer = optim.Adam(model.parameters())

# Load previously trained model weights and continue training. 
pretrained_weights_path = '/home/mriveraceron/glv-research/updated_results/GraphConv_ss/ss_5/model_weights.pth'
checkpoint = torch.load(pretrained_weights_path, weights_only=False)

# Current model weights, loss and optimizer
model.load_state_dict(checkpoint['model_state_dict'])
optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
epoch = checkpoint['epoch']
total_epochs = epoch + extra_epochs
loss = checkpoint['loss']
print(f"Loaded model weights from {pretrained_weights_path} | Current epoch: {epoch} | Current epoch loss: {loss:.6f}")

# Continue training
new_weights_path = f'{results_dir}/model_weights.pth'
loss_history, train_elapsed = training_fn(model, device, train_data[1:1000], new_weights_path, loss_fn, optimizer, extra_epochs, batch_size=30)

# Save the new parameters
torch.save({
    'epoch': total_epochs,
    'model_state_dict': model.state_dict(),
    'optimizer_state_dict': optimizer.state_dict(),
    'loss': loss_history[-1],
}, new_weights_path)


#------------------------------------------------------------
# Evaluate the model 
checkpoint = torch.load(new_weights_path, weights_only=False)
metrics, performance = evaluate_split(eval_data, model, device, batch_size=30)

# Generate summary
summarize(model, tr_size = len(train_data), ev_size = len(eval_data), 
          perf_tr = performance, 
          saving_path = f'{results_dir}/model_summary.txt', 
          elapsed = train_elapsed)

"""
Once we retrained our model we save the perdictions for later plotting or comparing purposes.
npz file will be composed of:
    - idxt: Node index with the maximum observed target value (e.g. 'keystoneness')
    - idxp: Node index with the maximum predicted target value (e.g. 'keystoneness')
    - mt: Observed true target values 
    - mp: Predicted true target values

"""
np.savez(f'{results_dir}/metrics.npz', **metrics._asdict(), **performance._asdict())
