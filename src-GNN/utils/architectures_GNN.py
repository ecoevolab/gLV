"""
Architectures for GNN
=================================================
Description:
   This script contains different GNN architectures used for training.
   It is intended to be sourced in future scenarios where training needs
   to be resumed or extended from a previous checkpoint.

Author: Manuel Rivera
Date:   May 1, 2026
"""

# Imports
import torch 
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv, SAGEConv


# Models
class GraphConv_model(nn.Module):
    def __init__(self, hidden_channels=64, num_layers=5):
        super().__init__()
        self.convs = nn.ModuleList()
        # First layer: 1 -> hidden_channels
        self.convs.append(GraphConv(13, hidden_channels))
        # Middle layers: hidden_channels -> hidden_channels
        for _ in range(num_layers - 2):
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
    
# SAGE architecture
class SAGE_model(nn.Module):
    def __init__(self, hidden_channels=64, num_layers=5):
        super().__init__()
        self.convs = nn.ModuleList()
        # First layer: 1 -> hidden_channels
        self.convs.append(SAGEConv(13, hidden_channels))
        # Middle layers: hidden_channels -> hidden_channels
        for _ in range(num_layers - 2):
            #self.convs.append(GATConv(hidden_channels*heads, hidden_channels, heads=heads))
            self.convs.append(SAGEConv(hidden_channels, hidden_channels))
        # Last layer: hidden_channels -> 1
        self.convs.append(SAGEConv(hidden_channels, 1))
    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        # Apply all layers except the last
        for i, conv in enumerate(self.convs[:-1]):
            x = conv(x, edge_index)
            x = F.relu(x)
        # Apply last layer with sigmoid
        x = self.convs[-1](x, edge_index)
        x = torch.sigmoid(x)
        return x  # [num_nodes]