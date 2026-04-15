"""
    Functions for GNN training and evaluation.

Description: 
    This script coontains all functions required for training, evaluating and summarizing training of GNNs.
    Functions can be imported and used in other scripts.

15-April-2026

How to import functions:
    import sys      
    sys.path.append("/path/to/your/folder")     # Specify functions directory path
    import my_functions      # Import functions
"""

#-----------------------
# Imports
import logging
import time
import random
import torch
import numpy as np
import sys
from collections import namedtuple
from scipy.stats import pearsonr, spearmanr
from tqdm import tqdm
from torch_geometric.utils import unbatch
from torch_geometric.loader import DataLoader

# Logging
log = logging.getLogger(__name__)
#-------------------------------
# Section: Evaluation function at last epoch
#-------------------------------

# Namedtuple
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
import torch
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