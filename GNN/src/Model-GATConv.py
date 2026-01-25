import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATConv

class ModelGATConv(nn.Module):
    def __init__(self, hidden_channels=64, num_layers=5, dropout=0.5, h = 1, concat=False):
        super().__init__()
        self.convs = nn.ModuleList()
        self.dropout = dropout
        # First layer: 1 -> hidden_channels
        self.convs.append(GATConv(in_channels = 1, out_channels = hidden_channels, heads = h, concat=concat))
        # Middle layers: hidden_channels -> hidden_channels
        for _ in range(num_layers - 2):
            self.convs.append(GATConv(in_channels = hidden_channels, out_channels = hidden_channels, heads = h, concat=concat))
        # Last layer: hidden_channels -> 1
        self.convs.append(GATConv(in_channels = hidden_channels, out_channels = 1, heads = h, concat=concat))
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_weights
        # Apply all layers except the last
        for i, conv in enumerate(self.convs[:-1]):
            x = conv(x, edge_index, edge_weight)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout)
        # Apply last layer with sigmoid
        x = self.convs[-1](x, edge_index, edge_weight)
        x = torch.sigmoid(x)
        return x  # [num_nodes]
