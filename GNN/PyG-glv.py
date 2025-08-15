import pandas as pd
import os
import torch 
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader  
from torch_geometric.nn import GCNConv
from sklearn.metrics import mean_squared_error, r2_score
import torch.nn.functional as F
import torch.optim as optim
import torch.nn as nn
import rpy2.robjects as ro
import numpy as np

#=======================  Data loading =======================
dir = "/home/mriveraceron/glv-research/Data/Exp06-D29-Apr"      # Parent-directory
ints_path = os.path.join(dir, "Interacts")                      # Interactions-directory
tsv_file = os.path.join(dir, "Exp06-D29-Apr.tsv")               # Parameters-table
mets_path = os.path.join(dir, "CP-info")                        # Metrics

# Get IDs from simulation from TSV table
ids =  pd.read_csv(tsv_file, sep="\t")["id"].tolist()           # Extract simulation IDs     

# List to store the Data objects
data_list = []
nodes =[]

id = "29fe5a3825"
for id in ids[:10]:
    #================= Declare paths =======================
    path_out = os.path.join(mets_path, "E_"+ id + "-Info.feather")
    path_in = os.path.join(ints_path, "Intrs_"+ id + ".feather")
    #================= Adjacency matrix (interactions) =======================                                                                         
    Inters_tensor = torch.tensor(pd.read_feather(path_in).round(4).to_numpy(), dtype=torch.float64)        # Interactions matrix conversion -> Tensor
    edges = torch.nonzero(Inters_tensor).t().contiguous()                                                  # Adjacency from->to
    weights = Inters_tensor[edges[0], edges[1]]                                                            # weights
    #================= Predictions (metrics) =======================
    input_nodes = Inters_tensor.shape[0]                                                                                    # Input-nodes
    pred_tensor = torch.tensor(pd.read_feather(path_out).drop(columns="id").round(4).to_numpy(), dtype=torch.float64)       # Prediction-features
    y = torch.zeros(input_nodes, pred_tensor.shape[1])                                                                      # shape [input_nodes, pred_features]
    y[:pred_tensor.shape[0], :] = pred_tensor                                                                               # fill first nodes with targets, rest zero
    #================= Input features =================
    n_species = pd.read_csv(tsv_file, sep="\t").query("id == @id")["n_species"].iloc[0]         # Extract number of species
    # print("# The number of species are:", n_species, " for simulation ", id, "\n")
    pop_seed = pd.read_csv(tsv_file, sep="\t").query("id == @id")["x0_seed"].iloc[0]            # Extract population seed of simulation
    ro.r(f"set.seed({pop_seed})")                                                               # Set seed (R)
    ro.r(f"x0 <- stats::runif({n_species}, min = 0.1, max = 1)")                                # Generate X0 (R)
    x0 = torch.tensor(np.array(ro.r("x0")).round(4), dtype=torch.float32)                       # Extract it to python
    #================= Save Data =================
    data = Data(x=x0.view(-1, 1), edge_index=edges, edge_attr=weights, y=y)
    data_list.append(data)

#=======================  GCN =======================
class GCNModel(torch.nn.Module):
    def __init__(self, in_features, OutIn_features, out_pred, pred_features):
        super(GCNModel, self).__init__()
        self.conv1 = GCNConv(in_features, OutIn_features)                       # Input layer --> GCNConv(in_channels: Size of each input sample, out_channels: Size of each out sample)
        self.conv2 = GCNConv(OutIn_features, out_pred)                          # Hidden layer
        self.output_layer = torch.nn.Linear(out_pred, pred_features)            # Regression
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_attr
        x = F.relu(self.conv1(x, data.edge_index, data.edge_weight))                   # Input layer                         
        x = F.relu(self.conv2(x, data.edge_index, data.edge_weight))                   # Layer 2
        out = self.output_layer(x)                                                     # Output layer
        return out
    

#======================= Move model and data to the device =======================
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# in_features -> Input features
# OutIn_features -> Number of features that exit on first hidden layer
# out_pred -> # of prediction features done (prediction)
# pred_features -> # of features to predict (true)
model = GCNModel(in_features=1, OutIn_features=5, out_pred=5, pred_features=5).to(device)


criterion = nn.MSELoss()                                                # Loss function for regression
optimizer = optim.Adam(model.parameters(), lr=0.01)                     # Optimizer
data_train = DataLoader(data_list[:8], batch_size=1, shuffle=True)         # Use batch_size > 1 if all graphs are compatible

model.train()                                               # Set model to training mode
for epoch in range(100):                                    # number of epochs
    total_loss = 0
    for data in data_train:
        data = data.to(device)                              # move data to gpu
        optimizer.zero_grad()                               # Clear gradients
        out = model(data)                                   # evaluation
        loss = criterion(out, data.y)                       # apply loss function
        loss.backward()                                     
        optimizer.step()                                    # Fix weights
        total_loss += loss.item()                           
    if epoch % 10 == 0:
        print(f'Epoch {epoch}, Loss: {loss.item():.6f}')

first_batch = next(iter(data_train))



#========================================== Set model to evaluation mode ==========================================
data_eval = DataLoader(data_list[8:], batch_size=1, shuffle=True)     # Use batch_size > 1 if all graphs are compatible
model.eval()  
all_preds = []
all_targets = []

with torch.no_grad():
    for data in data_eval:
        data = data.to(device)
        out = model(data)
        print("== this is the shape of predictions: ", out.shape, "\n")
        all_preds.append(out.cpu())
        print("== this is the shape of Y: ", data.y.shape, "\n")
        all_targets.append(data.y.cpu())
        print("Graph has y rows:", data.y.size(0))
        print("Model output rows:", out.size(0))


# Concatenate predictions and targets
print("== The shape of preditions are: ",  torch.cat(all_preds, dim=0).numpy().shape, "\n")
print("== The shape of Y is: ",  torch.cat(all_targets, dim=0).numpy().shape, "\n")

preds = torch.cat(all_preds, dim=0).numpy()
targets = torch.cat(all_targets, dim=0).numpy()

all_preds.numpy()
# Compute metrics
mse = mean_squared_error(targets, preds)
rmse = mse**0.5
r2 = r2_score(targets, preds)  # Closer to 1 is better

print(f"\n Evaluation:")
print(f"RMSE: {rmse:.4f}")
print(f"R² Score: {r2:.4f}, Accuracy: {r2*100:.2f}%")

#======================= ??????????????????? =========================
# Testing
model.train()                                           # Training mode
for epoch in range(10000):                              # # of epochs
    epoch = 1
    data = data.to(device)                              # move data to gpu
    optimizer.zero_grad()                               # Clear gradients
    out = model(data)                                   # Evaluation
    loss = criterion(out, data.y)                       # Loss-function
    loss.backward()                                     
    optimizer.step()                                    # Fix weights                           
    if epoch % 10 == 0:
        print(f'Epoch {epoch}, Loss: {loss.item():.6f}')




#=======================  Verification =======================

# Check if PyTorch is installed and if CUDA (GPU support) is available
print("PyTorch version:", torch.__version__)            # Prints the PyTorch version
print("CUDA available:", torch.cuda.is_available())     # Checks if CUDA is available (for GPU usage)
print("Number of GPUs:", torch.cuda.device_count())     # Number of available GPUs

# Check device properties (including number of cores and more)
gpu_properties = torch.cuda.get_device_properties(0)
print("GPU Properties:", gpu_properties)
