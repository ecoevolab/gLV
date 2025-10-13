## 1. Create batches
### 1.1 Create function
To create batches, we first define a function `load_single_data` that generates:

* The input tensor `x_tensor`
* The target tensor `y_tensor`
* The graph weights (`edge_weight`)  and adjacency matrix (`edge_index`) 
* Next, we store this data in a `Data` object from `torch_geometric.data`

```python
from torch_geometric.data import Data
import pyarrow.feather as feather

def load_single_data(id, A_dir, tgt_path):
    # Load adjacency matrix 
    A_path = os.path.join(A_dir, f"A_{id}.feather")
    A = pd.read_feather(A_path).to_numpy(dtype=np.float32)
    # Vector of edge weights
    row_idx, col_idx = np.nonzero(A)
    edge_weights = A[row_idx, col_idx]
    # Convert to torch tensors efficiently
    edge_index = torch.from_numpy(np.vstack([row_idx, col_idx]).astype(np.int64))
    edge_weights = torch.from_numpy(edge_weights)
    # Load target features 
    tgt_path = os.path.join(tgt_dir, f"tgt_{id}.feather")
    tgt_table =  feather.read_table(tgt_path, columns=['K_s'])
    y_tensor = torch.from_numpy(tgt_table.to_pandas().to_numpy(dtype=np.float32))   
    # Node features - simple ones vector
    n = A.shape[0]
    x_tensor = torch.ones(n, 1, dtype=torch.float32)
    # Clean up large intermediate
    del A, tgt_table
    # Create Data object
    data = Data(
        x=x_tensor,
        edge_weights=edge_weights,
        edge_index=edge_index,
        y=y_tensor
    )
    return data
```
### 1.2 Divide into training and validation set
We then split the data into training and validation sets by first shuffling the indices, and then assigning the first 80% to the training set and the remaining 20% to the validation set.

```python
import random

# Load all data samples (for demo, we use only first 100 samples)
indices = list(range(1, len(data_ids)))  # Indices 1-100
random.shuffle(indices)  # Uses Python's random module (already seeded)

# Now select first 80 for training, rest for validation
indx = round(len(indices) * .8)
train_indices = indices[:indx]            # First 80 shuffled indices
val_indices = indices[indx:]              # Last 20 shuffled indices
```
### 1.3 Parallelize function
For faster data generation, we define the function `generate_data_parallel`, which parallelizes `load_single_data` across 4 cores (`num_workers=4`).
```python
from concurrent.futures import ThreadPoolExecutor
import time

def generate_data_parallel(idx, A_dir, tgt_dir, num_workers=4):  # idx is a list
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        data_list = list(executor.map(
            load_single_data,           
            idx,                       # List of IDs to iterate over
            [A_dir]*len(idx),          # Repeat A_dir for each ID
            [tgt_dir]*len(idx)         # Repeat tgt_dir for each ID
        ))
    return data_list
```
To batch the data, we first create the directory where the batches will be saved (`path`).
```python
# Path were batches will be saved at
path = '/mnt/data/sur/users/mrivera/Train-sims/4379fd40-9f0a/batching'
os.makedirs(path, exist_ok=True)
```
Next, we calculate the number of batches needed to cover all training indices, generate the data, and save each batch to the specified `path`.
```python
def datgen_fn(num_batches, path, prefix, batch_ids_subset):
    for i in range(num_batches):
        # Batch ids
        start_idx = i * batch_size
        end_idx = min(start_idx + batch_size, len(batch_ids_subset))
        batch_ids = batch_ids_subset[start_idx:end_idx]
        # Generate data
        data = generate_data_parallel(batch_ids, A_dir, tgt_dir, num_workers=6)
        # Save data
        filename = os.path.join(path, f'{prefix}Batch_{i}.pt')
        torch.save(data, filename)
        # Free memory immediately
        del data
        print(f"Saved {prefix.lower()} batch {i}/{num_batches-1}: {filename}")

# Generate the TRAINING batches
batch_size = 1000
train_ids = data_ids[train_indices] 
num_batches = (len(batching_ids) + batch_size - 1) // batch_size
datgen_fn(num_batches, path, prefix = 'Train', batch_ids_subset = train_ids)

# Generate the VALIDATION batches
# Using the same batch size
# We recalculate the number of batches
val_ids = data_ids[val_indices] 
num_batches = (len(val_ids) + batch_size - 1) // batch_size
datgen_fn(num_batches, path, prefix = 'Val', batch_ids_subset = val_ids)
```

## 2. Copy data from fenix->GPU
We use secure shell protocol to copy files from the cluster to `data` directory.
```bash
# Copy files with secure shell protocol
# In this case we are interested in experiment 4379fd40-9f0a
scp ssh mrivera@fenix.lavis.unam.mx:/mnt/data/sur/users/mrivera/Train-sims/4379fd40-9f0a/batching /home/mriveraceron/data/4379fd40-9f0a
```