"""
Learning batching for GNN
=================================================
Purpose:
    - Generate data to learn batches.
    - I am interested on learning how to predict maximum data.y value per graph.

Author: Manuel Rivera
Date:   March 17, 2026
"""
#---------------------
# Generate data
#---------------------
import torch
from torch_geometric.data import Data

data_list = []
n = 5

for i in range(n):
    data = Data(y=torch.rand(5, 1), out=torch.rand(5, 1),  num_nodes=5)
    data_list.append(data)
    node = torch.argmax(data.y, dim=0).detach().cpu().numpy().item()
    print('>> Maximum node expected:', node )
    expected_value =  f'{data.y.max().cpu().numpy().item():.4f}'
    predicted_value = f"{data.y[node].cpu().numpy().item():.4f}"
    print(f'>> Maximum Y value expected ({expected_value}) vs predicted ({predicted_value}) \n')

#---------------------
# Using batches
#---------------------
from torch_geometric.utils import unbatch
from torch_geometric.loader import DataLoader

loader = DataLoader(data_list, batch_size=3, shuffle=True)

for batch in loader:
    y_list = unbatch(batch.y, batch.batch) 
    out_list = unbatch(batch.out, batch.batch)
    # out = model(batch)  # forward pass → shape [total_nodes, 1]
    # out_list = unbatch(out, batch.batch)
    for y, o in zip(y_list, out_list): 
        node = torch.argmax(y, dim=0).item()
        expected_value =  f'{y.max().cpu().numpy().item():.4f}'
        predicted_value = f'{y[node].cpu().numpy().item():.4f}'
        print(f'>> Maximum Y value expected ({expected_value}) vs predicted ({predicted_value}) \n')

