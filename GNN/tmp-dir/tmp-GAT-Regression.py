# Developing Graph Attention Network (GAT)
import torch
import torch.nn as nn
from torch_geometric.nn import GATConv
from torch_geometric.data import Data
import torch.nn.functional as F
import torch.optim as optim


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Generate input data
x = torch.tensor([
    [1.0, 0.3, 0.2],  # Species 0: interacts with 0,1,2
    [0.4, 1.0, 0.1],  # Species 1: interacts with 0,1,2
    [0.2, 0.5, 1.0],  # Species 2: interacts with 0,1,2
], dtype=torch.float)

edge_index = torch.tensor([
    [0, 1, 0, 2, 1, 2], # who
    [1, 0, 2, 0, 2, 1]  # to whom
], dtype=torch.long)

# Target for each node (e.g., extinction metrics)
y = torch.tensor([[2.0, 1.0, 0.0067],
                  [1.0, 0.0, 0.0050],
                  [3.0, 2.0, 0.0071]], dtype=torch.float)


# Number of samples
n_samples = 100

# Generate samples with slight noise
x_samples = x.unsqueeze(0) + 0.05 * torch.randn(n_samples, *x.shape)
print(x_samples.shape)
y_samples = y.unsqueeze(0) + 0.001 * torch.randn(n_samples, *y.shape)
print(y_samples.shape)

# Create a batch of graphs
batch = torch.zeros(x_samples.size(1), dtype=torch.long)  # If you have 3 nodes per graph

# Create the Data object for the batch
data = Data(x=x_samples.view(-1, x_samples.size(-1)),  # Flatten node features for batch
                  edge_index=edge_index,  # Same edge_index for all samples
                  y=y_samples.view(-1, y_samples.size(-1))  # Flatten target values for batch
                  ).to(device)

# Generate GAT
# in_channles --> Size of each input sample,
# out_channels --> Size of each output sample.
# heads --> Number of multi-head-attentions. 
class GATModel(torch.nn.Module):
    def __init__(self, in_channels, out_channels, heads, heads_out):
        super(GATModel, self).__init__()
        self.conv1 = GATConv(in_channels, out_channels, heads=heads, dropout=0.2, concat=True)
        # Features * number of attention heads
        self.conv2 = GATConv(out_channels * heads, out_channels, heads=heads_out, concat=False)
    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.elu(x)                                # Exponential Linear-unit
        x = self.conv2(x, edge_index)
        return x
    
  
# device = "cpu"

# Move model and data to the device
model = GATModel(in_channels=3, out_channels=3, heads=4, heads_out=1).to(device)

criterion = nn.MSELoss()                            # Loss function for regression
optimizer = optim.Adam(model.parameters(), lr=0.01) # Optimizer

model.train()  # Set model to training mode
for epoch in range(10000):  # number of epochs
    optimizer.zero_grad()           # Clear gradients
    out = model(data)  # Forward pass
    loss = criterion(out, data.y)         # Compute loss
    loss.backward()                # Backpropagation
    optimizer.step()               # Update weights
    if epoch % 100 == 0:
        print(f'Epoch {epoch}, Loss: {loss.item():.6f}')



# Set the model to evaluation mode (this disables dropout, etc.)
model.eval()

# Perform a forward pass (get the predictions)
with torch.no_grad():  # We don't need gradients for the test phase
    y_pred = model(data)

# Calculate the Mean Squared Error (MSE) for regression task
# Assuming y_pred and test_data.y are both of shape [num_nodes, num_features] 
from sklearn.metrics import mean_absolute_error
y_true = data.y
mse = mean_absolute_error(y_true.cpu().numpy(), y_pred.cpu().numpy())
percentage_mse = (mse / torch.mean(y_true).item()) * 100  # Scaling relative to the mean of the true values
print(f"Mean Squared Error (Percentage): {percentage_mse:.2f}%")
print(f"Mean Squared Error: {mse}")