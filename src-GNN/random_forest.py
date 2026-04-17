import torch
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from cuml.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from torch_geometric.data import Data

# For package versions
from importlib.metadata import packages_distributions, version
import importlib.metadata as metadata

# Logging
import os
import logging
import sys

# data generation
import glob

#--------------------------------------------
# Section: Setup logging 
#--------------------------------------------
results_dir = '/home/mriveraceron/glv-research/Results/RForest'
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

#--------------------------------------------
# Generate data
#--------------------------------------------
data_dir = '/home/mriveraceron/glv-research/data_null'

def data_generator(data_dir, split='train'):
    data_list = []
    paths = glob.glob(f'{data_dir}/*_{split}/*.pt')
    if not paths:
        raise FileNotFoundError(f"No .pt files found under {data_dir}/*_{split}/")
    for path in paths:
        all_data = torch.load(path, weights_only=False) # Data for GNN
        for data in all_data:
            new_data = Data(x=data.x, y=data.y.squeeze()) # Convert to Random Forest data
            data_list.extend([new_data])
    log.info(f"Total samples for {split}: {len(data_list)}")
    # Get directories of data for training
    used_dirs = '\n'.join(f'  {p}' for p in set(os.path.dirname(p) for p in paths))
    return data_list, used_dirs


train_data, train_dirs = data_generator(data_dir, split='train')
eval_data, eval_dirs = data_generator(data_dir, split='eval')

# Concatenate all training data
x_tr = torch.cat([d.x for d in train_data], dim=0).numpy()  # (169, 13)
y_tr = torch.cat([d.y for d in train_data], dim=0).numpy()  # (169,)
x_test = torch.cat([d.x for d in eval_data], dim=0).numpy()  # (169, 13)
y_test = torch.cat([d.y for d in eval_data], dim=0).numpy()  # (169,)

# Generate graph ids
tr_ids = np.concatenate([
    np.full(d.x.shape[0], i) for i, d in enumerate(train_data)
]) 
eval_ids = np.concatenate([
    np.full(d.x.shape[0], i) for i, d in enumerate(eval_data)
])
#------------------------
# Section: Generate random forest

# --- Train the Random Forest ---
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X, y)
preds = rf.predict(X_test)
rf.fit(X_train_np, y_train_np)

# --- Predict & evaluate ---
y_pred = rf.predict(X_test_np)
print(f"Accuracy: {accuracy_score(y_test_np, y_pred):.4f}")

# --- Convert predictions back to tensor if needed ---
y_pred_tensor = torch.tensor(y_pred)