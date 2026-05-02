import torch
import numpy as np
from sklearn.ensemble import  RandomForestRegressor
from sklearn.metrics import  mean_squared_error
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
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import time
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
to_np = lambda data, attr: torch.cat([getattr(d, attr) for d in data], dim=0).numpy()
#train_data = train_data[1:100]
#eval_data = eval_data[1:100]
x_tr, y_tr     = to_np(train_data, "x"), to_np(train_data, "y")
x_test, y_test = to_np(eval_data, "x"),  to_np(eval_data, "y")

# Generate graph ids
make_ids = lambda data: np.concatenate([np.full(d.x.shape[0], i) for i, d in enumerate(data)])
tr_ids   = make_ids(train_data)
eval_ids = make_ids(eval_data)

#------------------------
# Section: Generate random forest with different number of trees
estimator_range = range(100, 1000, 100)  # From 100-1000

# Generate data for results
rf_df = pd.DataFrame({
    "name": [f"rf_{n}" for n in estimator_range],
    "n_estimators": estimator_range,
    "rmse": 0.0,
    "pearson_corr": 0.0,
    "spearman_corr": 0.0,
    "node_acc": 0.0, 
    'elapsed_time': 0.0,
    'train_size': len(train_data),
    'eval_size': len(eval_data),
})

# Generate directories for results
weights_dir = f'{results_dir}/tree_params'
os.makedirs(weights_dir, exist_ok=True)

for i, row in rf_df.iterrows():
    start = time.time()
    n_estimators = row["n_estimators"]
    rf = RandomForestRegressor(n_estimators=n_estimators, n_jobs=-1, random_state=42)
    rf.fit(x_tr, y_tr)          # training
    #---------------------------
    # prediction with eval data
    y_pred = rf.predict(x_test) 
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    log.info(f"n_estimators: {n_estimators}, Test RMSE: {rmse:.4f}")
    #---------------------------
    # Save model
    # model_path = f'{weights_dir}/{row["name"]}.joblib'
    # joblib.dump(rf, model_path)
    # log.info(f"Random Forest model saved to {model_path}")
    #---------------------------
    # Predicted maximum node value per graph 
    results = pd.DataFrame({
        "y_pred": y_pred,
        "y_true": y_test,
        "graph_id": eval_ids
    })
    pred_max_per_graph = results.groupby("graph_id")["y_pred"].idxmax()
    true_max_per_graph = results.groupby("graph_id")["y_true"].idxmax()
    # Accuracy of predicting the graph with the maximum node value
    acc = (pred_max_per_graph == true_max_per_graph).mean()
    log.info(f"Graph-level accuracy (max node prediction): {acc:.4f}")
    #---------------------------
    # Generate correlations 
    pearson_corr, _ = pearsonr(y_test, y_pred)
    spearman_corr, _ = spearmanr(y_test, y_pred)
    log.info(f"Pearson correlation: {pearson_corr:.4f}")
    log.info(f"Spearman correlation: {spearman_corr:.4f}")
    #---------------------------
    # Save metrics
    np.savez(f'{weights_dir}/{row["name"]}_metrics.npz', 
        max_idx_true  = true_max_per_graph,
        max_idx_pred  = pred_max_per_graph,
        values_true   = y_test,
        values_pred   = y_pred,
    )
    #---------------------------
    # Add results to dataframe
    elapsed = time.time() - start
    rf_df.loc[i, ["rmse", "pearson_corr", "spearman_corr", "node_acc", "elapsed_time"]] = [rmse, pearson_corr, spearman_corr, acc, elapsed]
    rf_df.to_csv(f'{results_dir}/rf_results.csv', index=False)
    log.info(f"Tree completed. Number of estimators: {n_estimators}. Elapsed time: {elapsed:.2f}s\n")