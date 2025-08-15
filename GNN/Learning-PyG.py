#======================== Load data ========================
from torch_geometric.data import Data
from torch_geometric.nn import GATConv
from torch_geometric.datasets import Planetoid
import torch_geometric.transforms as T

import matplotlib.pyplot as plt

name_data = 'Cora'
dataset = Planetoid(root= '/tmp/' + name_data, name = name_data)
dataset.transform = T.NormalizeFeatures()
data = dataset[0]

print(f"Number of Classes in {name_data}:", dataset.num_classes)
print(f"Number of Node Features in {name_data}:", dataset.num_node_features)
print("Feature matrix shape:", data.x.shape)
print("First node's features:", data.x[0])

#============================= GAT =====================================
import numpy as np
import torch
from torch import nn as nn 
import torch.nn as nn
import torch.nn.functional as F

class GAT(torch.nn.Module):
    def __init__(self):
        super(GAT, self).__init__()
        self.hid = 8 # Output feature vector for each node after it has passed through an attention head.
        self.in_head = 8 #  Attention heads in the first GAT layer.  Processes the input features
        self.out_head = 1 # Attention heads at the second layer. Transforms the hidden layer's output to the final output.
        
        # Dropout to the input features to prevent overfitting, x% of features set to 0.
        self.conv1 = GATConv(dataset.num_features, self.hid, heads=self.in_head, dropout=0.6)
        self.conv2 = GATConv(self.hid*self.in_head, dataset.num_classes, concat=False,
                             heads=self.out_head, dropout=0.6)
    def forward(self, data):
        x, edge_index = data.x, data.edge_index  
        # training=self.training -> dropout is only applied during training.    
        x = F.dropout(x, p=0.6, training=self.training)
        x = self.conv1(x, edge_index)
        x = F.elu(x) #  ELU (Exponential Linear Unit)
        x = F.dropout(x, p=0.6, training=self.training)
        x = self.conv2(x, edge_index)
        
        return F.log_softmax(x, dim=1)
    
    
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# device = "cpu"

# Move model and data to the device
model = GAT().to(device)
data = dataset[0].to(device)


optimizer = torch.optim.Adam(model.parameters(), lr=0.005, weight_decay=5e-4)

model.train() # Set model in training phase
for epoch in range(1000):
    model.train()
    optimizer.zero_grad() # Dont accumulate gradients
    out = model(data)
    # Measures the difference between the predicted and true labels
    loss = F.nll_loss(out[data.train_mask], # Predicted labels
                      data.y[data.train_mask] # True labels
                    )
    
    if epoch%200 == 0:
        print(loss)
    
    loss.backward() # Compute how much each parameter contributed to the error.
    optimizer.step() # Optimize parameters
    
model.eval() # Set model on evaluation mode

# model(data): Runs a forward pass on the whole graph data.
# The model outputs a tensor of shape [num_nodes, num_classes] with log-probabilities.
# .max(dim=1): Finds the index (class label) of the highest score (log-probability) for each node.
# pred contains the predicted class for each node.
_, pred = model(data).max(dim=1)

# data.test_mask: Boolean mask indicating which nodes belong to the test set.
# pred[data.test_mask]: Predictions for the test nodes.
# data.y[data.test_mask]: True labels for the test nodes.
# .eq(...): Compares predicted labels to true labels → returns True where they match.
# .sum(): Counts how many predictions were correct.
# .item(): Converts the result from a tensor to a float.
correct = float(pred[data.test_mask].eq(data.y[data.test_mask]).sum().item())
acc = correct / data.test_mask.sum().item()
print('Accuracy: {:.4f}'.format(acc)) 