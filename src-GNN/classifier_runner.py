# 24-February-2024
# This code is for running the model using all node features.

#-------------------------------
# Section: Generate model
import torch 
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv

class model(nn.Module):
    def __init__(self, hidden_channels=64, num_layers=5):
        super().__init__()
        self.convs = nn.ModuleList()
        # First layer: 1 -> hidden_channels
        self.convs.append(GraphConv(13, hidden_channels))
        # Middle layers: hidden_channels -> hidden_channels
        for _ in range(num_layers - 2):
            #self.convs.append(GATConv(hidden_channels*heads, hidden_channels, heads=heads))
            self.convs.append(GraphConv(hidden_channels, hidden_channels))
        # As we have two classes (keystone or non). The last layer will be: 2
        self.convs.append(GraphConv(hidden_channels, 2))
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_weights
        # Apply all layers except the last
        for i, conv in enumerate(self.convs[:-1]):
            x = conv(x, edge_index, edge_weight)
            x = F.relu(x)
        # Last layer with softmax for probabilities
        x = self.convs[-1](x, edge_index, edge_weight)
        # x = F.softmax(x, dim=1)
        return x  # [num_nodes]
    

#-------------------------------
# Section: Traning loop
import glob
import time
import random
import pandas as pd
from tqdm import tqdm
import os

def training_loop(model_declared, device, batched_paths, weights_dir, loss_fn, optimizer, epochs=100):
    model_declared.train()
    loss_history  =  []             # Loss at epoch
    total_elapsed = 0               # Running time
    class_pred, class_true = [], []   # Declare list for metrics
    empty_count = 0
    for i, iter in enumerate(tqdm(range(epochs), desc="Training")):
        start = time.time()
        epoch_loss = 0
        for path in batched_paths:
            # path = batched_paths[1]
            data_list = torch.load(path, weights_only=False)          
            for data in data_list:
                if data.y.shape == 0:
                    empty_count += 1
                    continue
                #----------------------
                # Move it to device and run model
                # data = data_list[0]
                data = data.to(device)
                optimizer.zero_grad()
                out = model_declared(data)
                # Loss function
                loss = loss_fn(out, data.y.view(-1).long())
                loss.backward()
                #----------------------
                # Section: Last layeer information
                if iter == epochs - 1:  # Only for the last epoch
                    #----------------------
                    # Append true class
                    class_true.append(data.y.cpu().numpy())
                    # Append most probable class
                    class_pred.append(torch.argmax(out, dim = 1).cpu().numpy())
                    #----------------------
                    # Save weights
                    path_save = f'{weights_dir}/{os.path.basename(data_dir)}.pth'
                    torch.save(model_declared.state_dict(), path_save)
                #torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()
                epoch_loss += loss.item()   # Accumulate loss
        # Append epoch loss to history
        loss_history.append(epoch_loss)
        elapsed = time.time() - start
        total_elapsed += elapsed
        # Print every 25 epochs
        if iter % 10 == 0:
            tqdm.write(f"Epoch {iter}: Loss = {epoch_loss:.4f}, Elapsed time: {elapsed:.2f}")
    # Summary
    print(f"Found {empty_count} graphs with empty labels")
    print(f'>> the total elapsed time with {epochs} epochs is {total_elapsed:.2f} seconds ( {total_elapsed/60:.2f} minutes)')   
    return  loss_history, class_pred, class_true


#-------------------------------
# Section: Run model
import torch.optim as optim
import glob
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


# As cross entropy expects raw values:
import torch.nn as nn
loss_fn = nn.CrossEntropyLoss()
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
seed_fn(42)
model_declared = model(hidden_channels=64, num_layers=5).to(device)
optimizer = optim.Adam(model.parameters(), lr=0.001)
epochs = 100
print('The number of epochs will be:', epochs, '\n')

# Load file paths
data_dir = '/home/mriveraceron/glv-research/data_tensors/Boosted_filter-classes'
weights_dir = '/home/mriveraceron/glv-research/model_weights'
batched_paths = glob.glob(f"{data_dir}/*.pt")
loss_history, class_pred, class_true = training_loop(model_declared, device, batched_paths, weights_dir, loss_fn, optimizer, epochs)

#-------------------------------
# Check data

for path in batched_paths:
    data_list = torch.load(path, weights_only=False)          
    for data in data_list:
        if data.y.numel() == 0:
            print(f"Empty labels found in: {path}")
        elif data.y.dim() == 0:  # scalar tensor
            print(f"Scalar label found: {data.y} in {path}")
        else:
            print(data.y.shape)  # check what shapes you actually have
#-------------------------------
# Flatten lists of indexes
# Force each element to be 1D before concatenating
class_true = np.concatenate([np.atleast_1d(x) for x in class_true])
class_pred = np.concatenate([np.atleast_1d(x) for x in class_pred])

#-------------------------------
# Section: Save loss and max data/out tensor
result_path = '/home/mriveraceron/glv-research/Results/Boosted_keystone/Filtered_Classes'
os.makedirs(result_path, exist_ok=True)
print('The results directory will be:', result_path, '\n')
#-------------------------------
# Save max indexes
np.save(f'{result_path}/label_true.npy', class_true)
np.save(f'{result_path}/label_pred.npy', class_pred)
#-------------------------------
# Save loss history
np.save(f'{result_path}/loss_history.npy', np.array(loss_history))
# np.load(f'{result_path}/Dummy_loss.npy')
# x= np.load(f'{result_path}/Dummy_max_true.npy')


#-------------------------------
# Section: Plot loss over time
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 5))
plt.plot(loss_history)
plt.title("Loss over time")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.grid(True)
plt.savefig(f'{result_path}/DataLoss_plot.png')

#-------------------------------
# Section: Plot indexes predicted vs true
import matplotlib.pyplot as plt

# Calculate accuracy
correct_counts   = [np.sum((class_true == l) & (class_pred == l)) for l in [0, 1]]
incorrect_counts = [np.sum((class_true == l) & (class_pred != l)) for l in [0, 1]]

# Sanity check
[np.sum((class_true == 1) & (class_pred == 1))]
np.sum((class_true == 1))

fig, ax = plt.subplots(figsize=(8, 5))

x = np.array([0, 1])
width = 0.35

bars1 = ax.bar(x - width/2, correct_counts,   width, label='Correct',   color='steelblue')
bars2 = ax.bar(x + width/2, incorrect_counts, width, label='Incorrect', color='tomato')

# Value labels on bars
for bar in bars1 + bars2:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
            str(int(bar.get_height())), ha='center', va='bottom')

ax.set_xticks(x)
ax.set_xticklabels(['Class 0', 'Class 1'])
ax.set_ylabel('Count')
ax.set_title('Predictions per Class')
ax.legend()
plt.tight_layout()
plt.show()
plt.tight_layout()
plt.savefig(f'{result_path}/classes_plot.png',dpi=150)

# Testing 
out = torch.randn(3, 2)
tgt = torch.randint(0, 2, (3, 1)).squeeze()  # random, unrelated to out
loss_fn = torch.nn.CrossEntropyLoss()
loss_fn(out, tgt)
# torch.argmax(out, dim = 1)