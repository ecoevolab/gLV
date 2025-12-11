# Generate grid
from sklearn.model_selection import ParameterGrid
import numpy as np 

# Define your parameter grid
param_grid = {
    'hidden_channels': [64],
    'num_layers': [2,5,10],
    'dropout':[0, 0.5],
    'learning_rate': [.0001]
}

# Create all combinations
grid = ParameterGrid(param_grid)

# Iterate through all combinations
for params in grid:
    print(params)
    # Train your model with these params
    # model = GNN(hidden_channels=params['hidden_channels'], ...)


# ==============================
# Import my files
import glob

# List files
mods_files = glob.glob(f'/home/mriveraceron/glv-research/gLV/GNN/architectures/Model-*.py')
for m in mods_files:
    with open(m) as f:
        exec(f.read())


for params in grid:
    print(params)
    ch = params['hidden_channels']
    l = params['num_layers']
    drop = params['dropout']
    model_list = [ModelGATConv(hidden_channels=ch, num_layers=l, dropout=drop, h=5), 
                  ModelGCNConv(hidden_channels=ch, num_layers=l, dropout=drop),
                  ModelGraphConv(hidden_channels=ch, num_layers=l, dropout=drop),
    ]
    # Your training/evaluation code here
    for model in model_list:
        # train and evaluate model
        print(model)
             
# ==============================

