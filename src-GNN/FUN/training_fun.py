"""
Traning function for Graph Neural Network (GNN) 
=================================================
Purpose:
    Declare training function for gnn training. It employs DataLoader for faster training time.

Input Data:
    - model (torch.nn.Module): Previously declared model architecture (e.g., GCN, GraphConv, SAGEConv).
    - device (str | torch.device): Device to run the training on (e.g., 'cuda' or 'cpu').
    - data_train (List[torch_geometric.data.Data]): Training dataset, expected to be a list of Data() objects with x being the node features and y being the target variable.
    - weights_path (str | Path): Path to save the trained model weights.
    - loss_fn (Callable): Loss function to optimize (e.g., MSELoss, CrossEntropyLoss).
    - optimizer (torch.optim.Optimizer): Optimizer for training (e.g., Adam, SGD).
    - epochs (int): Number of training epochs.
    - batch_size (int, default=30): Number of samples per batch (default: 30).

Dependencies:
    torch==2.4.0+cu121, pandas==2.0.1, tqdm==4.66.5

Returns:
    - loss_history (np.ndarray): Array of loss values for each epoch.
    - total_elapsed (float): Total training time in seconds.

Author: Manuel Rivera
Date:   June 18, 2026
"""

# Import libraries
import os
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"  # Must be before torch import
import torch
import numpy as np
import time
from torch_geometric.loader import DataLoader
from tqdm import tqdm


# Declare function
def training_fn(model, device, data_train, weights_path, loss_fn, optimizer, epochs, batch_size=30):
    model.train()
    print(f'Starting training of model {model.__class__.__name__ } \n')
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
            print(f"Epoch {epoch}: Loss={epoch_loss:.6f} | Time={epoch_elapsed:.2f}s")
    # Save after all epochs
    torch.save({
        'epoch': epochs - 1,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'loss': epoch_loss
    }, weights_path)
    print(f"Total elapsed: {total_elapsed:.2f}s ({total_elapsed/60:.2f} min)")
    return np.array(loss_history), total_elapsed